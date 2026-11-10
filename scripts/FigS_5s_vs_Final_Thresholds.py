from __future__ import annotations

import re
from pathlib import Path

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

FS = 50.0
IMU_WINDOW_S = 0.5
ROBUST_K_CAL = 6.0

OUT_DIR = DATA_ROOT / "HEMOX_5s_vs_Final_Thresholds"
PNG_DPI = 300
SAVE_SVG = True
SAVE_CSV = True

# HEMOX / IEEE colors
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREY = "#777F87"
GRID = "#CFD5DC"
BLACK = "#222222"

# Four runs are horizontally separated within each participant.
RUN_OFFSETS = {
    "S01": -0.30,
    "S02": -0.10,
    "NS01": +0.10,
    "NS02": +0.30,
}


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
        "axes.labelsize": 8,
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

def compact_token(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", s.upper())


def infer_participant_run(path: Path):
    token = compact_token(str(path))

    pm = re.search(r"P(0[1-9]|10)", token)
    participant = f"P{pm.group(1)}" if pm else None

    # NS first so S01 is not accidentally matched inside NS01.
    rm = re.search(r"(NS01|NS02|S01|S02)", token)
    run = rm.group(1) if rm else None

    return participant, run


def choose_best(files: list[Path], suffix: str):
    if not files:
        return None

    def key(p: Path):
        exact = 0 if p.name.upper().endswith(suffix.upper()) else 1
        duplicate = 1 if "(" in p.name else 0
        return (exact, duplicate, len(str(p)))

    return sorted(files, key=key)[0]


def discover_files(root: Path):
    slots = {}

    for p in root.rglob("*"):
        if not p.is_file():
            continue

        name = p.name.upper()
        if name.endswith("_REPORT.TXT"):
            kind, suffix = "REPORT", "_REPORT.TXT"
        elif name.endswith("_RAW.CSV"):
            kind, suffix = "RAW", "_RAW.CSV"
        else:
            continue

        participant, run = infer_participant_run(p)
        if participant not in PARTICIPANTS or run not in RUNS:
            continue

        slot = slots.setdefault(
            (participant, run),
            {"REPORT": [], "RAW": []},
        )
        slot[kind].append(p)

    out = {}
    for key, slot in slots.items():
        report = choose_best(slot["REPORT"], "_REPORT.TXT")
        raw = choose_best(slot["RAW"], "_RAW.CSV")
        if report is not None and raw is not None:
            out[key] = {"REPORT": report, "RAW": raw}

    return out


# =============================================================================
# REPORT PARSING
# =============================================================================

def channel_block(text: str, ch: int) -> str:
    start = re.search(rf"^Channel\s+{ch}\s*$", text, flags=re.I | re.M)
    if not start:
        return ""

    tail = text[start.end():]
    end = re.search(
        r"^(Channel\s+[12]|Dual-channel fusion|Decision rule|Processing specification)\s*$",
        tail,
        flags=re.I | re.M,
    )
    return tail[:end.start()] if end else tail


def find_float(text: str, label: str) -> float:
    m = re.search(
        rf"{re.escape(label)}\s*:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
        text,
        flags=re.I,
    )
    return float(m.group(1)) if m else np.nan


def parse_calibration_interval(block: str):
    m = re.search(
        r"Calibration noise guard\s*:\s*.*?interval\s+"
        r"([-+]?\d+(?:\.\d+)?)\s*-\s*([-+]?\d+(?:\.\d+)?)\s*s",
        block,
        flags=re.I,
    )
    if not m:
        return np.nan, np.nan
    return float(m.group(1)), float(m.group(2))


def parse_report(report_path: Path):
    text = report_path.read_text(encoding="utf-8", errors="replace")

    vm = re.search(r"Algorithm version\s*:\s*([^\n\r]+)", text, flags=re.I)
    algorithm = vm.group(1).strip() if vm else "N/A"

    result = {"algorithm_version": algorithm}

    for ch in (1, 2):
        block = channel_block(text, ch)
        cal_start, cal_end = parse_calibration_interval(block)

        result[ch] = {
            "cal_start_s": cal_start,
            "cal_end_s": cal_end,
            "acc_final": find_float(
                block, "Effective acceleration threshold"
            ),
            "gyro_final": find_float(
                block, "Effective gyroscope threshold"
            ),
        }

    return result


# =============================================================================
# IMU FEATURES + 5-s CALIBRATION RULE
# =============================================================================

def robust_sigma(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan

    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return 1.4826 * mad


def rolling_std(x: np.ndarray, window: int):
    return (
        pd.Series(np.asarray(x, float))
        .rolling(
            window=window,
            center=True,
            min_periods=max(3, window // 3),
        )
        .std(ddof=0)
        .to_numpy()
    )


def rolling_rms(x: np.ndarray, window: int):
    return np.sqrt(
        pd.Series(np.asarray(x, float) ** 2)
        .rolling(
            window=window,
            center=True,
            min_periods=max(3, window // 3),
        )
        .mean()
        .to_numpy()
    )


def calibration_guard(x: np.ndarray):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]

    if x.size < 5:
        return np.nan, np.nan, np.nan, np.nan

    p99 = float(np.percentile(x, 99))
    med = float(np.median(x))
    sig = robust_sigma(x)

    threshold = max(
        p99,
        med + ROBUST_K_CAL * sig,
    )

    return float(threshold), p99, med, sig


def raw_time_seconds(df: pd.DataFrame):
    if "common_time_us" in df.columns:
        t = pd.to_numeric(
            df["common_time_us"],
            errors="coerce",
        ).to_numpy(float)

        if np.isfinite(t).sum() >= 2:
            return (t - np.nanmin(t)) / 1e6

    return np.arange(len(df), dtype=float) / FS


def read_imu_activity(raw_path: Path, ch: int):
    required = [
        "common_time_us",
        f"ch{ch}_acc_x_g",
        f"ch{ch}_acc_y_g",
        f"ch{ch}_acc_z_g",
        f"ch{ch}_gyro_x_dps",
        f"ch{ch}_gyro_y_dps",
        f"ch{ch}_gyro_z_dps",
    ]

    header = pd.read_csv(raw_path, nrows=0)
    available = [c for c in required if c in header.columns]

    df = pd.read_csv(
        raw_path,
        usecols=available,
        low_memory=False,
    )

    imu_required = required[1:]
    missing = [c for c in imu_required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{raw_path.name}: missing required IMU columns: {missing}"
        )

    t = raw_time_seconds(df)

    acc = np.column_stack([
        pd.to_numeric(
            df[f"ch{ch}_acc_{axis}_g"],
            errors="coerce",
        ).to_numpy(float)
        for axis in "xyz"
    ])

    gyro = np.column_stack([
        pd.to_numeric(
            df[f"ch{ch}_gyro_{axis}_dps"],
            errors="coerce",
        ).to_numpy(float)
        for axis in "xyz"
    ])

    # Acceleration magnitude includes gravity.
    acc_mag = np.sqrt(np.nansum(acc * acc, axis=1))

    # Gyroscope bias removal: subtract per-axis recording median.
    gyro_bias = np.nanmedian(gyro, axis=0)
    gyro_mag = np.sqrt(
        np.nansum((gyro - gyro_bias) ** 2, axis=1)
    )

    win = max(3, int(round(IMU_WINDOW_S * FS)))

    acc_activity = rolling_std(acc_mag, win)
    gyro_activity = rolling_rms(gyro_mag, win)

    return t, acc_activity, gyro_activity


# =============================================================================
# BUILD ALL-RUN TABLE
# =============================================================================

def build_table(root: Path):
    files = discover_files(root)
    rows = []

    for participant in PARTICIPANTS:
        for run in RUNS:
            pair = files.get((participant, run))
            if pair is None:
                print(f"[SKIP] {participant}/{run}: REPORT + RAW not found")
                continue

            try:
                report = parse_report(pair["REPORT"])

                for ch in (1, 2):
                    cfg = report[ch]

                    if not (
                        np.isfinite(cfg["cal_start_s"])
                        and np.isfinite(cfg["cal_end_s"])
                        and cfg["cal_end_s"] > cfg["cal_start_s"]
                    ):
                        raise ValueError(
                            f"CH{ch}: calibration interval not reported"
                        )

                    t, acc_activity, gyro_activity = read_imu_activity(
                        pair["RAW"],
                        ch,
                    )

                    cal_mask = (
                        np.isfinite(t)
                        & (t >= cfg["cal_start_s"])
                        & (t < cfg["cal_end_s"])
                    )

                    acc_cal, acc_p99, acc_med, acc_sig = calibration_guard(
                        acc_activity[cal_mask]
                    )
                    gyro_cal, gyro_p99, gyro_med, gyro_sig = calibration_guard(
                        gyro_activity[cal_mask]
                    )

                    rows.append({
                        "participant": participant,
                        "run": run,
                        "channel": f"CH{ch}",
                        "algorithm_version": report["algorithm_version"],
                        "calibration_start_s": cfg["cal_start_s"],
                        "calibration_end_s": cfg["cal_end_s"],
                        "calibration_duration_s":
                            cfg["cal_end_s"] - cfg["cal_start_s"],

                        "acc_calibration_guard_g": acc_cal,
                        "acc_cal_P99_g": acc_p99,
                        "acc_cal_median_g": acc_med,
                        "acc_cal_robust_sigma_g": acc_sig,
                        "acc_final_reported_g": cfg["acc_final"],

                        "gyro_calibration_guard_dps": gyro_cal,
                        "gyro_cal_P99_dps": gyro_p99,
                        "gyro_cal_median_dps": gyro_med,
                        "gyro_cal_robust_sigma_dps": gyro_sig,
                        "gyro_final_reported_dps": cfg["gyro_final"],

                        "report_file": str(pair["REPORT"]),
                        "raw_file": str(pair["RAW"]),
                    })

                print(f"[OK] {participant}/{run}")

            except Exception as exc:
                print(f"[ERROR] {participant}/{run}: {exc}")

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError(
            "No runs were processed. Check DATA_ROOT and folder structure."
        )

    return df


# =============================================================================
# PLOT
# =============================================================================

def x_position(participant: str, run: str):
    p = int(participant[1:]) - 1
    return p + RUN_OFFSETS[run]


def plot_figure(df: pd.DataFrame):
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(7.2, 4.8),
        sharex=True,
    )

    panels = [
        (
            axes[0],
            "acc_calibration_guard_g",
            "acc_final_reported_g",
            "Acceleration threshold (g)",
        ),
        (
            axes[1],
            "gyro_calibration_guard_dps",
            "gyro_final_reported_dps",
            "Gyroscope threshold (deg/s)",
        ),
    ]

    for ax, cal_col, final_col, ylabel in panels:

        for _, row in df.iterrows():
            cal = row[cal_col]
            final = row[final_col]

            if not np.isfinite(cal):
                continue

            x = x_position(
                row["participant"],
                row["run"],
            )

            color = (
                BLUE
                if row["channel"] == "CH1"
                else ORANGE
            )

            # Thin connector from 5-s guard to final effective ON threshold.
            if np.isfinite(final):
                ax.plot(
                    [x, x],
                    [cal, final],
                    color=color,
                    linewidth=0.75,
                    alpha=0.65,
                    zorder=1,
                )

            # 5-s calibration guard
            ax.scatter(
                x,
                cal,
                marker="o",
                s=31,
                facecolor=color,
                edgecolor="white",
                linewidth=0.45,
                zorder=3,
            )

            # Final reported ON
            if np.isfinite(final):
                ax.scatter(
                    x,
                    final,
                    marker="*",
                    s=55,
                    facecolor=color,
                    edgecolor=BLACK,
                    linewidth=0.45,
                    zorder=4,
                )

        # Faint separators only; NO P02/P03 dashed cohort line.
        for sep in np.arange(0.5, 9.5, 1.0):
            ax.axvline(
                sep,
                color=GRID,
                linewidth=0.55,
                alpha=0.75,
                zorder=0,
            )

        ax.grid(
            True,
            axis="y",
            color=GRID,
            linewidth=0.65,
            alpha=0.75,
        )

        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)

    # Panel labels
    axes[0].text(
        -0.075,
        1.01,
        "(a)",
        transform=axes[0].transAxes,
        fontsize=10,
        fontweight="bold",
    )

    axes[1].text(
        -0.075,
        1.01,
        "(b)",
        transform=axes[1].transAxes,
        fontsize=10,
        fontweight="bold",
    )

    axes[1].set_xlabel("Participant / run")
    axes[1].set_xticks(np.arange(10))
    axes[1].set_xticklabels(PARTICIPANTS)

    # Compact legend matching the requested figure.
    handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=GREY,
            markeredgecolor="none",
            markersize=6.5,
            label="Initial 5-s calibration",
        ),
        Line2D(
            [0], [0],
            marker="*",
            linestyle="none",
            markerfacecolor=GREY,
            markeredgecolor=BLACK,
            markersize=8,
            label="Final adaptive threshold",
        ),
        Line2D(
            [0], [0],
            color=BLUE,
            linewidth=2.0,
            label="CH1",
        ),
        Line2D(
            [0], [0],
            color=ORANGE,
            linewidth=2.0,
            label="CH2",
        ),
    ]

    axes[0].legend(
        handles=handles,
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.21),
        columnspacing=1.8,
        handletextpad=0.7,
    )

    fig.tight_layout(
        rect=(0.035, 0.03, 0.995, 0.94)
    )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    png = OUT_DIR / "FigS_5s_vs_Final_Thresholds_AllRuns.png"
    svg = OUT_DIR / "FigS_5s_vs_Final_Thresholds_AllRuns.svg"

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

    df = build_table(DATA_ROOT)

    if SAVE_CSV:
        csv_path = OUT_DIR / "FigS_5s_vs_Final_Thresholds_Data.csv"
        df.to_csv(
            csv_path,
            index=False,
        )
        print(f"\nSaved data: {csv_path}")

    png, svg = plot_figure(df)

    print("\nDone.")
    print(f"PNG: {png}")
    if SAVE_SVG:
        print(f"SVG: {svg}")

    # Verification of calibration duration.
    d = df["calibration_duration_s"].dropna().to_numpy(float)
    if d.size:
        n5 = int(
            np.isclose(
                d,
                5.0,
                atol=1e-6,
            ).sum()
        )
        print(
            f"Exactly 5.0-s calibration intervals: "
            f"{n5}/{len(d)} channel-records"
        )


if __name__ == "__main__":
    main()
