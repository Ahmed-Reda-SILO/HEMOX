
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D


# =============================================================================
# USER SETTINGS
# =============================================================================

DATA_ROOT = Path(__file__).resolve().parent

PARTICIPANTS = [f"P{i:02d}" for i in range(1, 11)]
RUNS = ["S01", "S02", "NS01", "NS02"]

# Preferred source behavior:
# AUTO      -> use exported threshold history when present, otherwise reconstruct
# EXPORTED  -> use only *_MOTION_THRESHOLDS.csv; skip missing runs
# RECONSTRUCT -> reconstruct every available run from RAW
SOURCE_MODE = "AUTO"

# Set False if you prefer AUTO to skip runs that lack *_MOTION_THRESHOLDS.csv
# instead of reconstructing them.
RECONSTRUCT_MISSING = True

OUT_DIR = DATA_ROOT / "Adaptive_Motion_Thresholds_AndroidAligned_ALL"
INDIVIDUAL_DIR = OUT_DIR / "Individual_Runs"

FIGSIZE_INDIVIDUAL = (7.2, 4.5)
FIGSIZE_COMBINED = (7.2, 4.8)

PNG_DPI = 600
SAVE_SVG = True
SAVE_CSV = True

# Combined plot:
# All run trajectories are faint; the across-run median is drawn thick.
SHOW_COMBINED_MEDIAN = True
COMBINED_RUN_ALPHA = 0.20
COMBINED_RUN_LINEWIDTH = 0.75
COMBINED_MEDIAN_LINEWIDTH = 2.3

# Android chart colors from AdaptiveThresholdPlotView.
CH1_COLOR = "#0058A0"    # Color.rgb(0, 88, 160)
CH2_COLOR = "#008755"    # Color.rgb(0, 135, 85)
GRID_COLOR = "#E2E9E6"
GREY = "#777F87"

# Individual plots use small points as Android does.
SHOW_POINTS = True


# =============================================================================
# EXACT ANDROID CONSTANTS
# =============================================================================

FS = 50
WINDOW_SECONDS = 8.0
STEP_SECONDS = 2.0

ACCEL_RMS_FLOOR_G = 0.0020
ACCEL_RMS_CEILING_G = 0.0100
GYRO_RMS_FLOOR_DPS = 0.15
GYRO_RMS_CEILING_DPS = 0.75

ADAPTIVE_BASELINE_FRACTION = 0.40
ADAPTIVE_SIGMA_MULTIPLIER = 5.0

CALIBRATION_SETTLING_SECONDS = 2.0
CALIBRATION_MIN_SECONDS = 5.0
CALIBRATION_SEARCH_SECONDS = 30.0
CALIBRATION_SIGMA_MULTIPLIER = 6.0
CALIBRATION_RELEASE_SIGMA_MULTIPLIER = 3.0
CALIBRATION_ACCEL_MAX_ON_G = 0.030
CALIBRATION_GYRO_MAX_ON_DPS = 3.5

MOTION_RELEASE_RATIO = 0.80

MOTION_TRACKING_LOOKBACK_SECONDS = 30.0
MOTION_TRACKING_UPDATE_SECONDS = 10.0
MOTION_TRACKING_MIN_SECONDS = 5.0
MOTION_TRACKING_CANDIDATE_MULTIPLIER = 1.50
MOTION_TRACKING_MIN_CANDIDATE_FRACTION = 0.60
MOTION_TRACKING_MAX_RISE_FRACTION = 0.20
MOTION_TRACKING_ON_QUANTILE = 0.95
MOTION_TRACKING_OFF_QUANTILE = 0.75

IMU_RMS_SECONDS = 0.5

MIN_IR_DC = 10000.0
MIN_RED_DC = 3000.0

# MotionAwareMetrics:
# 2 s contact + 2 s IMU settling + 5 s IMU calibration + 5 s stabilization
BIOMETRIC_START_SECONDS = 15.0


# =============================================================================
# STYLE
# =============================================================================

def setup_style():
    installed = {f.name for f in font_manager.fontManager.ttflist}
    preferred = ["Arial Narrow", "Arial", "Liberation Sans", "DejaVu Sans"]
    font = next((f for f in preferred if f in installed), "DejaVu Sans")

    plt.rcParams.update({
        "font.family": font,
        "font.size": 9,
        "font.weight": "bold",
        "axes.labelsize": 7.5,
        "axes.labelweight": "bold",
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


# =============================================================================
# FILE DISCOVERY
# =============================================================================

def compact_token(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def infer_participant_run(path: Path):
    """
    Robustly infer Pxx and run. NS01/NS02 are matched before S01/S02 so
    stationary runs are not accidentally confused with non-stationary runs.
    """
    full = compact_token(str(path))
    name = compact_token(path.name)

    pm = re.search(r"P(0[1-9]|10)", full)
    participant = f"P{pm.group(1)}" if pm else None

    # Prefer the basename because it normally carries the explicit run ID.
    for token in (name, full):
        rm = re.search(r"(NS01|NS02|S01|S02)", token)
        if rm:
            return participant, rm.group(1)

    return participant, None


def choose_best(paths: list[Path], suffix: str):
    if not paths:
        return None

    def key(p: Path):
        exact = 0 if p.name.upper().endswith(suffix.upper()) else 1
        duplicate = 1 if "(" in p.name else 0
        return (exact, duplicate, len(str(p)))

    return sorted(paths, key=key)[0]


def discover_all_sources(root: Path):
    slots = {
        (p, r): {"threshold": [], "raw": []}
        for p in PARTICIPANTS
        for r in RUNS
    }

    for path in root.rglob("*"):
        if not path.is_file():
            continue

        upper = path.name.upper()
        if not (
            upper.endswith("_MOTION_THRESHOLDS.CSV")
            or upper.endswith("_RAW.CSV")
        ):
            continue

        participant, run = infer_participant_run(path)
        if participant not in PARTICIPANTS or run not in RUNS:
            continue

        if upper.endswith("_MOTION_THRESHOLDS.CSV"):
            slots[(participant, run)]["threshold"].append(path)
        elif upper.endswith("_RAW.CSV"):
            slots[(participant, run)]["raw"].append(path)

    out = {}
    for key, item in slots.items():
        out[key] = {
            "threshold": choose_best(
                item["threshold"],
                "_MOTION_THRESHOLDS.CSV",
            ),
            "raw": choose_best(
                item["raw"],
                "_RAW.CSV",
            ),
        }
    return out


# =============================================================================
# ANDROID NUMERICAL HELPERS
# =============================================================================

def java_median(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    x = np.sort(x)
    middle = x.size // 2
    if x.size % 2 == 1:
        return float(x[middle])
    return float((x[middle - 1] + x[middle]) / 2.0)


def java_quantile(values, q):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    x = np.sort(x)
    position = (x.size - 1) * float(q)
    index = int(np.floor(position))
    fraction = position - index
    j = min(index + 1, x.size - 1)
    return float(x[index] + (x[j] - x[index]) * fraction)


def centered_standard_deviation(values, width):
    x = np.asarray(values, dtype=float)
    n = x.size
    radius = int(width // 2)

    valid = np.isfinite(x)
    xv = np.where(valid, x, 0.0)

    c_count = np.concatenate([[0], np.cumsum(valid.astype(np.int64))])
    c_sum = np.concatenate([[0.0], np.cumsum(xv)])
    c_sq = np.concatenate([[0.0], np.cumsum(xv * xv)])

    idx = np.arange(n)
    start = np.maximum(0, idx - radius)
    end = np.minimum(n, idx + radius + 1)

    count = c_count[end] - c_count[start]
    total = c_sum[end] - c_sum[start]
    total_sq = c_sq[end] - c_sq[start]

    out = np.full(n, np.nan, dtype=float)
    ok = count > 0

    mean = np.zeros(n, dtype=float)
    mean[ok] = total[ok] / count[ok]

    variance = np.zeros(n, dtype=float)
    variance[ok] = total_sq[ok] / count[ok] - mean[ok] ** 2
    variance = np.maximum(0.0, variance)

    out[ok] = np.sqrt(variance[ok])
    return out


def centered_rms(values, width):
    x = np.asarray(values, dtype=float)
    n = x.size
    radius = int(width // 2)

    valid = np.isfinite(x)
    xv = np.where(valid, x, 0.0)

    c_count = np.concatenate([[0], np.cumsum(valid.astype(np.int64))])
    c_sq = np.concatenate([[0.0], np.cumsum(xv * xv)])

    idx = np.arange(n)
    start = np.maximum(0, idx - radius)
    end = np.minimum(n, idx + radius + 1)

    count = c_count[end] - c_count[start]
    total_sq = c_sq[end] - c_sq[start]

    out = np.full(n, np.nan, dtype=float)
    ok = count > 0
    out[ok] = np.sqrt(total_sq[ok] / count[ok])
    return out


def adaptive_limit(values, floor, ceiling):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0:
        return np.nan, float(floor)

    x = np.sort(x)

    baseline_count = max(
        1,
        int(np.ceil(
            ADAPTIVE_BASELINE_FRACTION * x.size
        )),
    )
    baseline_values = x[:baseline_count]

    baseline = java_median(baseline_values)
    mad = java_median(np.abs(baseline_values - baseline))
    robust_sigma = 1.4826 * mad if np.isfinite(mad) else 0.0

    threshold = baseline + ADAPTIVE_SIGMA_MULTIPLIER * robust_sigma
    threshold = min(float(ceiling), max(float(floor), float(threshold)))

    return float(baseline), float(threshold)


def calibration_stats(values, start, end):
    start = max(0, int(start))
    end = min(len(values), int(end))

    x = np.asarray(values[start:end], dtype=float)
    x = x[np.isfinite(x)]

    required = max(5, (end - start) // 2)
    if x.size < required:
        return None

    x = np.sort(x)
    middle = java_median(x)
    robust_sigma = 1.4826 * java_median(np.abs(x - middle))

    return {
        "median": float(middle),
        "robust_sigma": float(robust_sigma),
        "p75": java_quantile(x, 0.75),
        "p95": java_quantile(x, 0.95),
        "p99": java_quantile(x, 0.99),
    }


def calibration_motion_limits(
    accel_rms,
    gyro_rms,
    accel_adaptive,
    gyro_adaptive,
):
    settling = max(
        0,
        int(round(CALIBRATION_SETTLING_SECONDS * FS)),
    )
    segment = max(
        1,
        int(round(CALIBRATION_MIN_SECONDS * FS)),
    )
    search_end = min(
        len(accel_rms),
        int(round(CALIBRATION_SEARCH_SECONDS * FS)),
    )
    step = max(1, FS)

    best_score = np.inf
    best_accel = None
    best_gyro = None
    best_start = -1

    for start in range(
        settling,
        search_end - segment + 1,
        step,
    ):
        accel = calibration_stats(accel_rms, start, start + segment)
        gyro = calibration_stats(gyro_rms, start, start + segment)

        if accel is None or gyro is None:
            continue

        accel_guard = max(
            accel["p99"],
            accel["median"]
            + CALIBRATION_SIGMA_MULTIPLIER * accel["robust_sigma"],
        )
        gyro_guard = max(
            gyro["p99"],
            gyro["median"]
            + CALIBRATION_SIGMA_MULTIPLIER * gyro["robust_sigma"],
        )

        if (
            accel_guard > CALIBRATION_ACCEL_MAX_ON_G
            or gyro_guard > CALIBRATION_GYRO_MAX_ON_DPS
        ):
            continue

        score = (
            accel_guard / ACCEL_RMS_FLOOR_G
            + gyro_guard / GYRO_RMS_FLOOR_DPS
        )

        if score < best_score:
            best_score = score
            best_accel = accel
            best_gyro = gyro
            best_start = start

    accel_baseline, accel_adaptive_threshold = accel_adaptive
    gyro_baseline, gyro_adaptive_threshold = gyro_adaptive

    if best_start < 0:
        return {
            "accel_on": accel_adaptive_threshold,
            "gyro_on": gyro_adaptive_threshold,
            "accel_off": max(
                MOTION_RELEASE_RATIO * accel_adaptive_threshold,
                accel_baseline,
            ),
            "gyro_off": max(
                MOTION_RELEASE_RATIO * gyro_adaptive_threshold,
                gyro_baseline,
            ),
        }

    accel_cal_on = max(
        best_accel["p99"],
        best_accel["median"]
        + CALIBRATION_SIGMA_MULTIPLIER * best_accel["robust_sigma"],
    )
    gyro_cal_on = max(
        best_gyro["p99"],
        best_gyro["median"]
        + CALIBRATION_SIGMA_MULTIPLIER * best_gyro["robust_sigma"],
    )

    accel_on = max(accel_adaptive_threshold, accel_cal_on)
    gyro_on = max(gyro_adaptive_threshold, gyro_cal_on)

    accel_off = max(
        MOTION_RELEASE_RATIO * accel_on,
        max(
            best_accel["p95"],
            best_accel["median"]
            + CALIBRATION_RELEASE_SIGMA_MULTIPLIER
            * best_accel["robust_sigma"],
        ),
    )
    gyro_off = max(
        MOTION_RELEASE_RATIO * gyro_on,
        max(
            best_gyro["p95"],
            best_gyro["median"]
            + CALIBRATION_RELEASE_SIGMA_MULTIPLIER
            * best_gyro["robust_sigma"],
        ),
    )

    return {
        "accel_on": float(accel_on),
        "gyro_on": float(gyro_on),
        "accel_off": float(accel_off),
        "gyro_off": float(gyro_off),
    }


def dynamic_motion_limits(
    ir,
    red,
    accel_rms,
    gyro_rms,
    initial,
    scientific_start,
):
    n = len(accel_rms)

    current_accel_on = float(initial["accel_on"])
    current_gyro_on = float(initial["gyro_on"])
    current_accel_off = float(initial["accel_off"])
    current_gyro_off = float(initial["gyro_off"])

    lookback = max(
        1,
        int(round(MOTION_TRACKING_LOOKBACK_SECONDS * FS)),
    )
    update_step = max(
        1,
        int(round(MOTION_TRACKING_UPDATE_SECONDS * FS)),
    )
    minimum = max(
        1,
        int(round(MOTION_TRACKING_MIN_SECONDS * FS)),
    )

    next_update = min(
        n,
        int(scientific_start) + minimum,
    )

    for i in range(n):
        if i >= next_update:
            start = max(
                int(scientific_start),
                i - lookback,
            )
            sl = slice(start, i)

            contact = (
                np.isfinite(ir[sl])
                & np.isfinite(red[sl])
                & (ir[sl] >= MIN_IR_DC)
                & (red[sl] >= MIN_RED_DC)
                & np.isfinite(accel_rms[sl])
                & np.isfinite(gyro_rms[sl])
            )

            contact_count = int(contact.sum())

            low_amplitude = (
                contact
                & (
                    accel_rms[sl]
                    <= MOTION_TRACKING_CANDIDATE_MULTIPLIER
                    * current_accel_on
                )
                & (
                    gyro_rms[sl]
                    <= MOTION_TRACKING_CANDIDATE_MULTIPLIER
                    * current_gyro_on
                )
            )

            candidate_count = int(low_amplitude.sum())

            required = max(
                minimum,
                int(np.ceil(
                    MOTION_TRACKING_MIN_CANDIDATE_FRACTION
                    * contact_count
                )),
            )

            if (
                contact_count >= minimum
                and candidate_count >= required
            ):
                accel_candidates = np.sort(accel_rms[sl][low_amplitude])
                gyro_candidates = np.sort(gyro_rms[sl][low_amplitude])

                accel = calibration_stats(
                    accel_candidates, 0, len(accel_candidates)
                )
                gyro = calibration_stats(
                    gyro_candidates, 0, len(gyro_candidates)
                )

                if accel is not None and gyro is not None:
                    desired_accel_on = max(
                        initial["accel_on"],
                        max(
                            java_quantile(
                                accel_candidates,
                                MOTION_TRACKING_ON_QUANTILE,
                            ),
                            accel["median"]
                            + CALIBRATION_SIGMA_MULTIPLIER
                            * accel["robust_sigma"],
                        ),
                    )
                    desired_gyro_on = max(
                        initial["gyro_on"],
                        max(
                            java_quantile(
                                gyro_candidates,
                                MOTION_TRACKING_ON_QUANTILE,
                            ),
                            gyro["median"]
                            + CALIBRATION_SIGMA_MULTIPLIER
                            * gyro["robust_sigma"],
                        ),
                    )

                    desired_accel_on = min(
                        CALIBRATION_ACCEL_MAX_ON_G,
                        desired_accel_on,
                    )
                    desired_gyro_on = min(
                        CALIBRATION_GYRO_MAX_ON_DPS,
                        desired_gyro_on,
                    )

                    next_accel_on = min(
                        desired_accel_on,
                        current_accel_on
                        * (1.0 + MOTION_TRACKING_MAX_RISE_FRACTION),
                    )
                    next_gyro_on = min(
                        desired_gyro_on,
                        current_gyro_on
                        * (1.0 + MOTION_TRACKING_MAX_RISE_FRACTION),
                    )

                    current_accel_on = max(
                        current_accel_on,
                        next_accel_on,
                    )
                    current_gyro_on = max(
                        current_gyro_on,
                        next_gyro_on,
                    )

                    current_accel_off = min(
                        current_accel_on,
                        max(
                            MOTION_RELEASE_RATIO * current_accel_on,
                            max(
                                java_quantile(
                                    accel_candidates,
                                    MOTION_TRACKING_OFF_QUANTILE,
                                ),
                                accel["median"]
                                + CALIBRATION_RELEASE_SIGMA_MULTIPLIER
                                * accel["robust_sigma"],
                            ),
                        ),
                    )
                    current_gyro_off = min(
                        current_gyro_on,
                        max(
                            MOTION_RELEASE_RATIO * current_gyro_on,
                            max(
                                java_quantile(
                                    gyro_candidates,
                                    MOTION_TRACKING_OFF_QUANTILE,
                                ),
                                gyro["median"]
                                + CALIBRATION_RELEASE_SIGMA_MULTIPLIER
                                * gyro["robust_sigma"],
                            ),
                        ),
                    )

            next_update += update_step

    return {
        "final_accel_on": float(current_accel_on),
        "final_gyro_on": float(current_gyro_on),
    }


def motion_profile(
    ir,
    red,
    ax,
    ay,
    az,
    gx,
    gy,
    gz,
):
    n = len(ir)

    accel_magnitude = np.sqrt(
        ax * ax + ay * ay + az * az
    )

    gyro_bias_x = java_median(gx)
    gyro_bias_y = java_median(gy)
    gyro_bias_z = java_median(gz)

    gyro_magnitude = np.sqrt(
        (gx - gyro_bias_x) ** 2
        + (gy - gyro_bias_y) ** 2
        + (gz - gyro_bias_z) ** 2
    )

    rms_samples = max(
        5,
        int(round(IMU_RMS_SECONDS * FS)),
    )

    accel_rms = centered_standard_deviation(
        accel_magnitude,
        rms_samples,
    )
    gyro_rms = centered_rms(
        gyro_magnitude,
        rms_samples,
    )

    accel_limit = adaptive_limit(
        accel_rms,
        ACCEL_RMS_FLOOR_G,
        ACCEL_RMS_CEILING_G,
    )
    gyro_limit = adaptive_limit(
        gyro_rms,
        GYRO_RMS_FLOOR_DPS,
        GYRO_RMS_CEILING_DPS,
    )

    calibrated = calibration_motion_limits(
        accel_rms,
        gyro_rms,
        accel_limit,
        gyro_limit,
    )

    scientific_start = min(
        n,
        int(round(BIOMETRIC_START_SECONDS * FS)),
    )

    dynamic = dynamic_motion_limits(
        ir,
        red,
        accel_rms,
        gyro_rms,
        calibrated,
        scientific_start,
    )

    return {
        "accel_threshold_g": dynamic["final_accel_on"],
        "gyro_threshold_dps": dynamic["final_gyro_on"],
    }


# =============================================================================
# RAW RECONSTRUCTION OF Android motionThresholdHistory()
# =============================================================================

def read_raw(raw_path: Path):
    needed = [
        "ch1_ir_raw",
        "ch1_red_raw",
        "ch1_acc_x_g",
        "ch1_acc_y_g",
        "ch1_acc_z_g",
        "ch1_gyro_x_dps",
        "ch1_gyro_y_dps",
        "ch1_gyro_z_dps",
        "ch2_ir_raw",
        "ch2_red_raw",
        "ch2_acc_x_g",
        "ch2_acc_y_g",
        "ch2_acc_z_g",
        "ch2_gyro_x_dps",
        "ch2_gyro_y_dps",
        "ch2_gyro_z_dps",
    ]

    df = pd.read_csv(raw_path, low_memory=False)

    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(
            "RAW file is missing required columns:\n"
            + "\n".join(missing)
        )

    out = {}
    for c in needed:
        out[c] = pd.to_numeric(
            df[c],
            errors="coerce",
        ).to_numpy(float)

    return out


def channel_prefix(raw, ch, end):
    return {
        "ir": raw[f"ch{ch}_ir_raw"][:end],
        "red": raw[f"ch{ch}_red_raw"][:end],
        "ax": raw[f"ch{ch}_acc_x_g"][:end],
        "ay": raw[f"ch{ch}_acc_y_g"][:end],
        "az": raw[f"ch{ch}_acc_z_g"][:end],
        "gx": raw[f"ch{ch}_gyro_x_dps"][:end],
        "gy": raw[f"ch{ch}_gyro_y_dps"][:end],
        "gz": raw[f"ch{ch}_gyro_z_dps"][:end],
    }


def reconstruct_android_history(raw_path: Path):
    raw = read_raw(raw_path)
    n = len(raw["ch1_ir_raw"])

    first_samples = max(
        int(round(WINDOW_SECONDS * FS)),
        FS,
    )
    step_samples = max(
        1,
        int(round(STEP_SECONDS * FS)),
    )

    rows = []
    end = min(first_samples, n)

    while end <= n:
        p1 = motion_profile(**channel_prefix(raw, 1, end))
        p2 = motion_profile(**channel_prefix(raw, 2, end))

        rows.append({
            "elapsed_s": (end - 1.0) / FS,
            "ch1_accel_threshold_g": p1["accel_threshold_g"],
            "ch2_accel_threshold_g": p2["accel_threshold_g"],
            "ch1_gyro_threshold_dps": p1["gyro_threshold_dps"],
            "ch2_gyro_threshold_dps": p2["gyro_threshold_dps"],
            "mode": "bounded_adaptive_reconstructed",
        })

        if end == n:
            break

        end = min(n, end + step_samples)

    return pd.DataFrame(rows)


# =============================================================================
# LOAD EXPORTED HISTORY
# =============================================================================

def load_exported_history(path: Path):
    df = pd.read_csv(path, low_memory=False)

    required = [
        "elapsed_s",
        "ch1_accel_threshold_g",
        "ch2_accel_threshold_g",
        "ch1_gyro_threshold_dps",
        "ch2_gyro_threshold_dps",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "MOTION_THRESHOLDS.csv is missing columns:\n"
            + "\n".join(missing)
        )

    for c in required:
        df[c] = pd.to_numeric(
            df[c],
            errors="coerce",
        )

    keep = required + (["mode"] if "mode" in df.columns else [])

    return (
        df[keep]
        .dropna(subset=["elapsed_s"])
        .sort_values("elapsed_s")
        .reset_index(drop=True)
    )


# =============================================================================
# RUN LOADING
# =============================================================================

def load_one_run(participant, run, files):
    mode = SOURCE_MODE.strip().upper()

    if mode not in {"AUTO", "EXPORTED", "RECONSTRUCT"}:
        raise ValueError(
            "SOURCE_MODE must be AUTO, EXPORTED, or RECONSTRUCT"
        )

    threshold_path = files.get("threshold")
    raw_path = files.get("raw")

    if mode in {"AUTO", "EXPORTED"} and threshold_path is not None:
        df = load_exported_history(threshold_path)
        source = f"Android exported: {threshold_path.name}"
        return df, source

    if mode == "EXPORTED":
        raise FileNotFoundError(
            "No *_MOTION_THRESHOLDS.csv"
        )

    if mode == "AUTO" and not RECONSTRUCT_MISSING:
        raise FileNotFoundError(
            "No *_MOTION_THRESHOLDS.csv and RECONSTRUCT_MISSING=False"
        )

    if raw_path is None:
        raise FileNotFoundError(
            "No matching *_RAW.csv"
        )

    df = reconstruct_android_history(raw_path)
    source = f"Android-algorithm reconstruction: {raw_path.name}"
    return df, source


# =============================================================================
# INDIVIDUAL PLOTS
# =============================================================================

def plot_individual(
    df: pd.DataFrame,
    participant: str,
    run: str,
):
    fig, axes = plt.subplots(
        2,
        1,
        figsize=FIGSIZE_INDIVIDUAL,
        sharex=True,
    )

    t = df["elapsed_s"].to_numpy(float)

    panels = [
        (
            axes[0],
            "ch1_accel_threshold_g",
            "ch2_accel_threshold_g",
            "Acceleration threshold (g RMS)",
            0.0025,
        ),
        (
            axes[1],
            "ch1_gyro_threshold_dps",
            "ch2_gyro_threshold_dps",
            "Gyroscope threshold (deg/s RMS)",
            0.2,
        ),
    ]

    for ax, c1, c2, ylabel, min_ymax in panels:
        y1 = df[c1].to_numpy(float)
        y2 = df[c2].to_numpy(float)

        kwargs1 = {
            "color": CH1_COLOR,
            "linewidth": 1.8,
            "label": "CH1",
        }
        kwargs2 = {
            "color": CH2_COLOR,
            "linewidth": 1.8,
            "label": "CH2",
        }

        if SHOW_POINTS:
            kwargs1.update(marker="o", markersize=2.7)
            kwargs2.update(marker="o", markersize=2.7)

        ax.plot(t, y1, **kwargs1)
        ax.plot(t, y2, **kwargs2)

        finite_values = np.concatenate([
            y1[np.isfinite(y1)],
            y2[np.isfinite(y2)],
        ])
        if finite_values.size:
            ax.set_ylim(
                0.0,
                max(
                    min_ymax,
                    1.15 * float(np.max(finite_values)),
                ),
            )

        ax.set_ylabel(ylabel)
        ax.grid(
            True,
            color=GRID_COLOR,
            linewidth=0.65,
            alpha=0.85,
        )
        ax.legend(
            frameon=False,
            loc="upper left",
            ncol=2,
        )

    axes[0].text(
        -0.075, 1.01, "(a)",
        transform=axes[0].transAxes,
        fontsize=10,
        fontweight="bold",
    )
    axes[1].text(
        -0.075, 1.01, "(b)",
        transform=axes[1].transAxes,
        fontsize=10,
        fontweight="bold",
    )

    axes[1].set_xlabel("Elapsed time (s)")

    fig.text(
        0.985,
        0.985,
        f"{participant}/{run}",
        ha="right",
        va="top",
        fontsize=8.5,
        fontweight="bold",
    )

    fig.tight_layout(
        rect=(0.035, 0.025, 0.995, 0.965)
    )

    INDIVIDUAL_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    stem = (
        INDIVIDUAL_DIR
        / f"Adaptive_Motion_Thresholds_{participant}_{run}"
    )

    png = stem.with_suffix(".png")
    svg = stem.with_suffix(".svg")

    fig.savefig(
        png,
        dpi=PNG_DPI,
        bbox_inches="tight",
    )

    if SAVE_SVG:
        fig.savefig(
            svg,
            bbox_inches="tight",
        )

    plt.close(fig)

    return png, svg


# =============================================================================
# COMBINED PLOT
# =============================================================================

def common_time_grid(run_data: list[pd.DataFrame]):
    """
    Build the union of available elapsed times rounded to milliseconds.
    Android histories normally use the same 2-s grid, but union+interpolation
    keeps the combined summary robust to shorter or incomplete runs.
    """
    times = []
    for df in run_data:
        t = pd.to_numeric(
            df["elapsed_s"],
            errors="coerce",
        ).to_numpy(float)
        times.extend(
            np.round(
                t[np.isfinite(t)],
                3,
            ).tolist()
        )

    if not times:
        return np.array([], dtype=float)

    return np.array(
        sorted(set(times)),
        dtype=float,
    )


def interpolate_without_extrapolation(t, y, grid):
    t = np.asarray(t, float)
    y = np.asarray(y, float)

    ok = np.isfinite(t) & np.isfinite(y)
    t = t[ok]
    y = y[ok]

    out = np.full(
        grid.shape,
        np.nan,
        dtype=float,
    )

    if t.size < 2:
        return out

    inside = (
        (grid >= np.min(t))
        & (grid <= np.max(t))
    )

    out[inside] = np.interp(
        grid[inside],
        t,
        y,
    )
    return out


def plot_combined(all_runs: list[dict]):
    fig, axes = plt.subplots(
        2,
        1,
        figsize=FIGSIZE_COMBINED,
        sharex=True,
    )

    panels = [
        (
            axes[0],
            "ch1_accel_threshold_g",
            "ch2_accel_threshold_g",
            "Acceleration threshold (g RMS)",
            0.0025,
        ),
        (
            axes[1],
            "ch1_gyro_threshold_dps",
            "ch2_gyro_threshold_dps",
            "Gyroscope threshold (deg/s RMS)",
            0.2,
        ),
    ]

    for ax, c1, c2, ylabel, min_ymax in panels:
        all_values = []

        # All individual trajectories.
        for item in all_runs:
            df = item["df"]
            t = df["elapsed_s"].to_numpy(float)
            y1 = df[c1].to_numpy(float)
            y2 = df[c2].to_numpy(float)

            ax.plot(
                t,
                y1,
                color=CH1_COLOR,
                linewidth=COMBINED_RUN_LINEWIDTH,
                alpha=COMBINED_RUN_ALPHA,
            )
            ax.plot(
                t,
                y2,
                color=CH2_COLOR,
                linewidth=COMBINED_RUN_LINEWIDTH,
                alpha=COMBINED_RUN_ALPHA,
            )

            all_values.extend(
                y1[np.isfinite(y1)].tolist()
            )
            all_values.extend(
                y2[np.isfinite(y2)].tolist()
            )

        # Median trajectory across all available runs.
        if SHOW_COMBINED_MEDIAN:
            grid = common_time_grid(
                [item["df"] for item in all_runs]
            )

            if grid.size:
                ch1_matrix = []
                ch2_matrix = []

                for item in all_runs:
                    df = item["df"]
                    t = df["elapsed_s"].to_numpy(float)

                    ch1_matrix.append(
                        interpolate_without_extrapolation(
                            t,
                            df[c1].to_numpy(float),
                            grid,
                        )
                    )
                    ch2_matrix.append(
                        interpolate_without_extrapolation(
                            t,
                            df[c2].to_numpy(float),
                            grid,
                        )
                    )

                ch1_matrix = np.vstack(ch1_matrix)
                ch2_matrix = np.vstack(ch2_matrix)

                with np.errstate(all="ignore"):
                    med1 = np.nanmedian(
                        ch1_matrix,
                        axis=0,
                    )
                    med2 = np.nanmedian(
                        ch2_matrix,
                        axis=0,
                    )

                ax.plot(
                    grid,
                    med1,
                    color=CH1_COLOR,
                    linewidth=COMBINED_MEDIAN_LINEWIDTH,
                    label="CH1 median",
                    zorder=5,
                )
                ax.plot(
                    grid,
                    med2,
                    color=CH2_COLOR,
                    linewidth=COMBINED_MEDIAN_LINEWIDTH,
                    label="CH2 median",
                    zorder=5,
                )

        if all_values:
            ax.set_ylim(
                0.0,
                max(
                    min_ymax,
                    1.10 * float(np.max(all_values)),
                ),
            )

        ax.set_ylabel(ylabel)
        ax.grid(
            True,
            color=GRID_COLOR,
            linewidth=0.65,
            alpha=0.85,
        )

    axes[0].text(
        -0.075, 1.01, "(a)",
        transform=axes[0].transAxes,
        fontsize=10,
        fontweight="bold",
    )
    axes[1].text(
        -0.075, 1.01, "(b)",
        transform=axes[1].transAxes,
        fontsize=10,
        fontweight="bold",
    )

    axes[1].set_xlabel("Elapsed time (s)")

    legend_handles = [
        Line2D(
            [0], [0],
            color=CH1_COLOR,
            linewidth=COMBINED_MEDIAN_LINEWIDTH,
            label="CH1 median",
        ),
        Line2D(
            [0], [0],
            color=CH2_COLOR,
            linewidth=COMBINED_MEDIAN_LINEWIDTH,
            label="CH2 median",
        ),
        Line2D(
            [0], [0],
            color=GREY,
            linewidth=COMBINED_RUN_LINEWIDTH,
            alpha=COMBINED_RUN_ALPHA,
            label="Individual runs",
        ),
    ]

    axes[0].legend(
        handles=legend_handles,
        frameon=False,
        loc="upper left",
        ncol=3,
    )

    fig.text(
        0.985,
        0.985,
        f"All available runs (n={len(all_runs)})",
        ha="right",
        va="top",
        fontsize=8.5,
        fontweight="bold",
    )

    fig.tight_layout(
        rect=(0.035, 0.025, 0.995, 0.965)
    )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    stem = (
        OUT_DIR
        / "Adaptive_Motion_Thresholds_ALL_RUNS_COMBINED"
    )

    png = stem.with_suffix(".png")
    svg = stem.with_suffix(".svg")

    fig.savefig(
        png,
        dpi=PNG_DPI,
        bbox_inches="tight",
    )

    if SAVE_SVG:
        fig.savefig(
            svg,
            bbox_inches="tight",
        )

    plt.close(fig)

    return png, svg


# =============================================================================
# MAIN
# =============================================================================

def main():
    setup_style()

    if not DATA_ROOT.exists():
        raise FileNotFoundError(
            f"DATA_ROOT does not exist:\n{DATA_ROOT}"
        )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    INDIVIDUAL_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    sources = discover_all_sources(
        DATA_ROOT
    )

    all_runs = []
    status_rows = []
    combined_rows = []

    for participant in PARTICIPANTS:
        for run in RUNS:
            files = sources.get(
                (participant, run),
                {},
            )

            try:
                df, source = load_one_run(
                    participant,
                    run,
                    files,
                )

                if df.empty:
                    raise RuntimeError(
                        "Threshold history is empty"
                    )

                # Individual figure
                png, svg = plot_individual(
                    df,
                    participant,
                    run,
                )

                all_runs.append({
                    "participant": participant,
                    "run": run,
                    "source": source,
                    "df": df,
                })

                export = df.copy()
                export.insert(
                    0,
                    "run",
                    run,
                )
                export.insert(
                    0,
                    "participant",
                    participant,
                )
                export.insert(
                    2,
                    "source",
                    source,
                )
                combined_rows.append(export)

                if SAVE_CSV:
                    out_csv = (
                        INDIVIDUAL_DIR
                        / f"Adaptive_Motion_Thresholds_{participant}_{run}_DATA.csv"
                    )
                    export.to_csv(
                        out_csv,
                        index=False,
                    )

                status_rows.append({
                    "participant": participant,
                    "run": run,
                    "status": "OK",
                    "source": source,
                    "points": len(df),
                    "figure_png": str(png),
                })

                print(
                    f"[OK] {participant}/{run} | "
                    f"{len(df)} points | {source}"
                )

            except Exception as exc:
                status_rows.append({
                    "participant": participant,
                    "run": run,
                    "status": "SKIPPED",
                    "source": "",
                    "points": 0,
                    "figure_png": "",
                    "reason": str(exc),
                })

                print(
                    f"[SKIP] {participant}/{run}: {exc}"
                )

    if not all_runs:
        raise RuntimeError(
            "No threshold histories were available. "
            "Check DATA_ROOT and source files."
        )

    combined_png, combined_svg = plot_combined(
        all_runs
    )

    if SAVE_CSV:
        all_data = pd.concat(
            combined_rows,
            ignore_index=True,
        )
        all_data.to_csv(
            OUT_DIR
            / "Adaptive_Motion_Thresholds_ALL_RUNS_DATA.csv",
            index=False,
        )

    pd.DataFrame(
        status_rows
    ).to_csv(
        OUT_DIR
        / "Adaptive_Motion_Thresholds_Run_Status.csv",
        index=False,
    )

    print()
    print("=" * 68)
    print("HEMOX adaptive threshold analysis completed")
    print("=" * 68)
    print(
        f"Runs plotted individually: {len(all_runs)}/"
        f"{len(PARTICIPANTS) * len(RUNS)}"
    )
    print(
        f"Individual folder: {INDIVIDUAL_DIR}"
    )
    print(
        f"Combined PNG: {combined_png}"
    )
    if SAVE_SVG:
        print(
            f"Combined SVG: {combined_svg}"
        )


if __name__ == "__main__":
    main()
