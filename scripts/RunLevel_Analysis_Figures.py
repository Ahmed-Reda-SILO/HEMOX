from __future__ import annotations

import argparse
from io import BytesIO, StringIO
from pathlib import Path, PurePosixPath
import re
import warnings
from zipfile import ZipFile, BadZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

# ============================== USER SETTINGS ==============================
DATA_FOLDER = r""       # Blank = folder containing this script
OUTPUT_FOLDER = r""     # Blank = Run_Level_R2_Standalone_Output under data root
PAIR_TOLERANCE_S = 1.0
FROZEN_M1_OFFSET_PP = 0.906754   # native fused SpO2 minus reference, P01/P02
DPI = 600
IMAGE_FORMAT = "jpg"    # jpg or png
# ===========================================================================

PARTICIPANTS = tuple(f"P{i:02d}" for i in range(3, 11))
RUNS = ("S01", "S02", "NS01", "NS02")
CONDITION = {
    "S01": "Stationary",
    "S02": "Stationary",
    "NS01": "Controlled motion",
    "NS02": "Controlled motion",
}

BLUE = "#0072B2"
RED = "#C62828"
GREY = "#777F87"
ORANGE = "#D55E00"
GRID = "#CFD5DC"
INK = "#253348"

MEASURES = (
    dict(variable="SpO2", ref="sp20_spo2", pred="selected_spo2",
         valid="spo2_pair_valid", unit="pp", color=BLUE,
         ylabel=r"SpO$_2$ RMSE (pp)", band=3.0),
    dict(variable="HR", ref="sp20_hr", pred="selected_hr",
         valid="hr_pair_valid", unit="bpm", color=RED,
         ylabel="HR RMSE (bpm)", band=5.0),
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", type=Path,
                   default=Path(DATA_FOLDER) if DATA_FOLDER.strip() else None)
    p.add_argument("--out-dir", type=Path,
                   default=Path(OUTPUT_FOLDER) if OUTPUT_FOLDER.strip() else None)
    p.add_argument("--tolerance", type=float, default=PAIR_TOLERANCE_S)
    p.add_argument("--spo2-offset", type=float, default=FROZEN_M1_OFFSET_PP)
    p.add_argument("--dpi", type=int, default=DPI)
    p.add_argument("--format", choices=("jpg", "png"), default=IMAGE_FORMAT)
    return p.parse_args(argv)


def normal_subject(v):
    m = re.fullmatch(r"P\s*0*(\d{1,2})", str(v).strip(), re.I)
    return f"P{int(m.group(1)):02d}" if m and 1 <= int(m.group(1)) <= 10 else None


def run_token(name: str, participant: str):
    text = str(name).replace("\\", "/").upper()
    n = int(participant[1:])
    m = re.search(rf"P0*{n}[_-]?(NS0[12]|S0[12])(?=$|[^A-Z0-9])", text, re.I)
    return m.group(1).upper() if m else None


class Source:
    """Read either an extracted participant directory or one ZIP."""
    def __init__(self, path: Path):
        self.path = path.resolve()
        self.kind = "directory" if self.path.is_dir() else "zip"
        if self.kind == "directory":
            self.names = [str(p) for p in self.path.rglob("*") if p.is_file()]
        else:
            try:
                with ZipFile(self.path) as z:
                    self.names = [n for n in z.namelist() if not n.endswith("/")]
            except BadZipFile as exc:
                raise ValueError(f"Bad ZIP: {self.path}") from exc

    def read(self, member) -> bytes:
        if self.kind == "directory":
            return Path(member).read_bytes()
        with ZipFile(self.path) as z:
            return z.read(member)


def locate_participant_source(root: Path, participant: str) -> Path:
    """Prefer extracted Pxx folder. Otherwise accept one unambiguous Pxx ZIP."""
    folders = [p for p in root.iterdir()
               if p.is_dir() and p.name.strip().upper() == participant]
    if len(folders) == 1:
        return folders[0]
    if len(folders) > 1:
        raise ValueError(f"{participant}: multiple folders found: {folders}")

    zips = [p for p in root.iterdir()
            if p.is_file() and p.suffix.lower() == ".zip"
            and re.fullmatch(re.escape(participant) + r"(?:\(\d+\))?", p.stem, re.I)]
    exact = [p for p in zips if p.stem.upper() == participant]
    if len(exact) == 1:
        return exact[0]
    if len(zips) == 1:
        return zips[0]
    if not zips:
        raise FileNotFoundError(f"{participant}: no extracted folder or participant ZIP under {root}")
    raise ValueError(f"{participant}: ambiguous ZIP versions {zips}. Extract/select one version.")


def index_run_files(src: Source, participant: str):
    out = {r: {"hemo": [], "sp20": []} for r in RUNS}
    for name in src.names:
        base = PurePosixPath(str(name).replace("\\", "/")).name
        run = run_token(base, participant)
        if run not in RUNS:
            continue
        upper = base.upper()
        if upper.endswith("_FUSED_TRENDS.CSV") and "CALIBRATED" not in upper:
            out[run]["hemo"].append(name)
        elif upper.startswith("SP20") and upper.endswith(".CSV"):
            out[run]["sp20"].append(name)
    return out


def choose_unique(items, label, participant, run):
    if len(items) != 1:
        raise FileNotFoundError(
            f"{participant}/{run}: expected exactly one {label}; found {len(items)}: {items}")
    return items[0]


def decode_text(payload: bytes) -> str:
    if payload.startswith((b"\xff\xfe", b"\xfe\xff")) or b"\x00" in payload[:200]:
        return payload.decode("utf-16")
    try:
        return payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        return payload.decode("cp1252")


def read_hemox(payload: bytes):
    df = pd.read_csv(BytesIO(payload), low_memory=False)
    required = {"elapsed_s", "fused_spo2_percent", "fused_hr_bpm", "quality_valid"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"FUSED_TRENDS missing columns: {sorted(missing)}")

    elapsed = pd.to_numeric(df["elapsed_s"], errors="coerce")
    if elapsed.isna().any() or elapsed.duplicated().any() or not elapsed.is_monotonic_increasing:
        raise ValueError("FUSED_TRENDS elapsed_s must be numeric, unique, and increasing")

    quality = df["quality_valid"].astype(str).str.strip().str.lower().isin(
        ("true", "1", "yes", "valid"))

    out = pd.DataFrame({
        "hemox_elapsed_s": elapsed.astype(float),
        "quality_valid": quality.astype(bool),
        "hemox_spo2_m0": pd.to_numeric(df["fused_spo2_percent"], errors="coerce"),
        "hemox_hr_m0": pd.to_numeric(df["fused_hr_bpm"], errors="coerce"),
    })
    # Same R2 correlation convention: evaluate only quality-valid fused outputs.
    out = out.loc[out["hemox_elapsed_s"].between(0, 300) & out["quality_valid"]].reset_index(drop=True)
    return out


def parse_sp20(payload: bytes):
    """Robustly parse SP-20 export and express time relative to its first row."""
    text = decode_text(payload)
    lines = text.splitlines()

    # Find the actual table header rather than assuming a fixed line count.
    header_idx = None
    for i, line in enumerate(lines):
        clean = line.lstrip("#").strip().lower()
        if "time" in clean and "spo2" in clean and ("pr" in clean or "pulse" in clean):
            header_idx = i
            break
    if header_idx is None:
        # Fallback compatible with the standard one-line SP-20 metadata header.
        header_idx = 1

    csv_text = "\n".join(lines[header_idx:])
    df = pd.read_csv(StringIO(csv_text), sep=None, engine="python", dtype=str)
    df.columns = [str(c).strip().lstrip("#").strip() for c in df.columns]

    # Locate columns with tolerant naming.
    time_col = next((c for c in df.columns if c.strip().lower() == "time"), None)
    spo2_col = next((c for c in df.columns if "spo2" in c.lower()), None)
    hr_col = next((c for c in df.columns if re.search(r"\bpr\b|pulse|hr", c, re.I)), None)
    if not all((time_col, spo2_col, hr_col)):
        raise ValueError(f"SP-20 columns not recognized: {list(df.columns)}")

    raw_time = df[time_col].astype(str).str.strip()
    times = pd.to_datetime(raw_time, errors="coerce")
    if times.isna().all():
        formats = (
            "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %I:%M:%S %p",
            "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %I:%M:%S %p",
            "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %I:%M:%S %p",
        )
        for fmt in formats:
            candidate = pd.to_datetime(raw_time, format=fmt, errors="coerce")
            if candidate.notna().sum() > times.notna().sum():
                times = candidate
    keep = times.notna()
    if not keep.any():
        raise ValueError("No readable SP-20 timestamps")

    times = times.loc[keep].reset_index(drop=True)
    spo2 = pd.to_numeric(df.loc[keep, spo2_col], errors="coerce").reset_index(drop=True)
    hr = pd.to_numeric(df.loc[keep, hr_col], errors="coerce").reset_index(drop=True)

    if times.duplicated().any() or not times.is_monotonic_increasing:
        raise ValueError("SP-20 timestamps must be unique and increasing")

    elapsed = (times - times.iloc[0]).dt.total_seconds().astype(float)
    out = pd.DataFrame({
        "sp20_elapsed_s": elapsed,
        "sp20_timestamp": times.astype(str),
        "sp20_spo2": spo2,
        "sp20_hr": hr,
    })
    return out.loc[out["sp20_elapsed_s"].between(0, 300)].reset_index(drop=True)


def pair_one_to_one(sp20: pd.DataFrame, hemo: pd.DataFrame, tolerance: float):
    """SP-20-first, one-to-one nearest elapsed-time pairing used by the R2 correlation analysis."""
    if sp20.empty or hemo.empty:
        return pd.DataFrame()

    h = hemo.sort_values("hemox_elapsed_s").reset_index(drop=True)
    s = sp20.sort_values("sp20_elapsed_s").reset_index(drop=True)
    ht = h["hemox_elapsed_s"].to_numpy(float)
    used = np.zeros(len(h), dtype=bool)

    low = max(float(s.sp20_elapsed_s.iloc[0]), float(ht[0]))
    high = min(float(s.sp20_elapsed_s.iloc[-1]), float(ht[-1]))
    if high < low:
        return pd.DataFrame()

    rows = []
    for x in s.itertuples(index=False):
        rt = float(x.sp20_elapsed_s)
        if not (low <= rt <= high):
            continue
        at = int(np.searchsorted(ht, rt))
        left, right = max(0, at - 2), min(len(ht), at + 3)
        choices = [i for i in range(left, right)
                   if not used[i] and low <= ht[i] <= high]
        if not choices:
            continue
        i = min(choices, key=lambda j: (abs(ht[j] - rt), j))
        delta = float(ht[i] - rt)
        if abs(delta) > tolerance:
            continue
        used[i] = True
        y = h.iloc[i]
        rows.append({
            "sp20_elapsed_s": rt,
            "sp20_timestamp": x.sp20_timestamp,
            "hemox_elapsed_s": float(y.hemox_elapsed_s),
            "delta_t_s": delta,
            "quality_valid": bool(y.quality_valid),
            "sp20_spo2": x.sp20_spo2,
            "sp20_hr": x.sp20_hr,
            "hemox_spo2_m0": y.hemox_spo2_m0,
            "hemox_hr_m0": y.hemox_hr_m0,
        })
    return pd.DataFrame(rows)


def build_pairs(root: Path, tolerance: float, spo2_offset: float):
    all_pairs = []
    qa = []

    for participant in PARTICIPANTS:
        src_path = locate_participant_source(root, participant)
        src = Source(src_path)
        indexed = index_run_files(src, participant)
        print(f"{participant}: {src_path}")

        for run in RUNS:
            q = {"participant": participant, "run": run, "source": str(src_path)}
            try:
                he_name = choose_unique(indexed[run]["hemo"], "FUSED_TRENDS", participant, run)
                sp_name = choose_unique(indexed[run]["sp20"], "SP-20 CSV", participant, run)
                he = read_hemox(src.read(he_name))
                sp = parse_sp20(src.read(sp_name))
                paired = pair_one_to_one(sp, he, tolerance)
                if paired.empty:
                    raise ValueError("No HEMOX/SP-20 pairs within tolerance")

                paired["participant"] = participant
                paired["run"] = run
                paired["split"] = "Heldout"
                paired["selected_spo2"] = np.clip(
                    pd.to_numeric(paired["hemox_spo2_m0"], errors="coerce") - spo2_offset,
                    0, 100)
                paired["selected_hr"] = pd.to_numeric(paired["hemox_hr_m0"], errors="coerce")

                # Metric-specific validity, matching the exported R2 pair table logic.
                paired["sp20_spo2"] = pd.to_numeric(paired["sp20_spo2"], errors="coerce")
                paired["sp20_hr"] = pd.to_numeric(paired["sp20_hr"], errors="coerce")
                paired["spo2_pair_valid"] = (
                    paired["quality_valid"] &
                    paired["sp20_spo2"].between(0, 100) &
                    paired["selected_spo2"].between(0, 100) &
                    np.isfinite(paired["sp20_spo2"]) &
                    np.isfinite(paired["selected_spo2"])
                )
                paired["hr_pair_valid"] = (
                    paired["quality_valid"] &
                    paired["sp20_hr"].gt(0) & paired["selected_hr"].gt(0) &
                    np.isfinite(paired["sp20_hr"]) & np.isfinite(paired["selected_hr"])
                )

                all_pairs.append(paired)
                q.update(status="OK", hemox_quality_rows=len(he), sp20_rows=len(sp),
                         paired_rows=len(paired),
                         spo2_valid_pairs=int(paired.spo2_pair_valid.sum()),
                         hr_valid_pairs=int(paired.hr_pair_valid.sum()),
                         median_abs_delta_ms=1000 * float(np.median(np.abs(paired.delta_t_s))))
                print(f"  {run}: pairs={len(paired):3d}, "
                      f"SpO2={q['spo2_valid_pairs']:3d}, HR={q['hr_valid_pairs']:3d}")
            except Exception as exc:
                q.update(status=f"ERROR: {exc}", paired_rows=0,
                         spo2_valid_pairs=0, hr_valid_pairs=0)
                print(f"  {run}: ERROR: {exc}")
            qa.append(q)

    qa_df = pd.DataFrame(qa)
    bad = qa_df.loc[~qa_df.status.eq("OK")]
    if len(qa_df) != len(PARTICIPANTS) * len(RUNS) or not bad.empty:
        raise RuntimeError(
            "Incomplete P03--P10 source set. See source_run_QA.csv after output folder is created.\n" +
            (bad[["participant", "run", "status"]].to_string(index=False) if not bad.empty else ""))
    if not all_pairs:
        raise ValueError("No paired observations were generated")

    pairs = pd.concat(all_pairs, ignore_index=True)
    return pairs, qa_df


def errors(group, m):
    selected = group.loc[group[m["valid"]], [m["ref"], m["pred"]]]
    return selected[m["pred"]].to_numpy(float) - selected[m["ref"]].to_numpy(float)


def stats(e, band):
    if not len(e):
        return dict(n_valid_pairs=0, bias=np.nan, MAE=np.nan, RMSE=np.nan,
                    within_limit_pct=np.nan)
    return dict(
        n_valid_pairs=int(len(e)),
        bias=float(np.mean(e)),
        MAE=float(np.mean(np.abs(e))),
        RMSE=float(np.sqrt(np.mean(np.square(e)))),
        within_limit_pct=float(100 * np.mean(np.abs(e) <= band)),
    )


def calculate(df):
    run_rows, participant_rows, condition_rows, pooled_rows = [], [], [], []
    for m in MEASURES:
        for subject in PARTICIPANTS:
            person = df.loc[df.participant.eq(subject)]
            for run in RUNS:
                subset = person.loc[person.run.eq(run)]
                e = errors(subset, m)
                run_rows.append(dict(participant=subject, run=run,
                    condition=CONDITION[run], variable=m["variable"], unit=m["unit"],
                    status="OK" if len(e) else "N/A: no eligible pairs",
                    **stats(e, m["band"])))
            participant_rows.append(dict(participant=subject,
                variable=m["variable"], unit=m["unit"],
                **stats(errors(person, m), m["band"])))
            for condition in ("Stationary", "Controlled motion"):
                chosen = person.loc[person.run.map(CONDITION).eq(condition)]
                condition_rows.append(dict(participant=subject, condition=condition,
                    variable=m["variable"], unit=m["unit"],
                    contributing_runs=int(chosen.loc[chosen[m["valid"]], "run"].nunique()),
                    **stats(errors(chosen, m), m["band"])))
        pooled_rows.append(dict(variable=m["variable"], unit=m["unit"],
                                **stats(errors(df, m), m["band"])))

    run_df = pd.DataFrame(run_rows)
    participant_df = pd.DataFrame(participant_rows)
    condition_df = pd.DataFrame(condition_rows)
    pooled_df = pd.DataFrame(pooled_rows)

    contrasts = []
    for m in MEASURES:
        c = condition_df.loc[condition_df.variable.eq(m["variable"])]
        s = c.loc[c.condition.eq("Stationary")].set_index("participant")
        ns = c.loc[c.condition.eq("Controlled motion")].set_index("participant")
        for subject in PARTICIPANTS:
            a, b = s.loc[subject], ns.loc[subject]
            good = np.isfinite(a.RMSE) and np.isfinite(b.RMSE)
            contrasts.append(dict(
                participant=subject, variable=m["variable"], unit=m["unit"],
                stationary_RMSE=a.RMSE, motion_RMSE=b.RMSE,
                motion_minus_stationary=(b.RMSE-a.RMSE if good else np.nan),
                stationary_n=int(a.n_valid_pairs), motion_n=int(b.n_valid_pairs),
            ))
    return run_df, participant_df, condition_df, pooled_df, pd.DataFrame(contrasts)


def setup_style(dpi):
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((s for s in ("Arial Narrow", "Arial", "Liberation Sans", "DejaVu Sans")
                 if s in installed), "DejaVu Sans")
    plt.rcParams.update({
        "font.family": font, "font.weight": "bold", "axes.labelweight": "bold",
        "axes.labelsize": 12, "xtick.labelsize": 10.5, "ytick.labelsize": 11,
        "axes.linewidth": 1.3, "axes.spines.top": False, "axes.spines.right": False,
        "savefig.dpi": dpi, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    return font


def save(fig, path, dpi):
    opts = dict(dpi=dpi, facecolor="white", bbox_inches="tight")
    if path.suffix.lower() in (".jpg", ".jpeg"):
        opts["pil_kwargs"] = dict(quality=96, subsampling=0)
    fig.savefig(path, **opts)
    plt.close(fig)
    print("Saved:", path)


def heatmap(run_df, out, dpi):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.5))
    for ax, m, letter in zip(axes, MEASURES, "ab"):
        d = (run_df.loc[run_df.variable.eq(m["variable"])]
             .pivot(index="participant", columns="run", values="RMSE")
             .reindex(index=PARTICIPANTS, columns=RUNS))
        z = d.to_numpy(float)
        vmax = np.nanmax(z) if np.isfinite(z).any() else 1.0
        cmap = LinearSegmentedColormap.from_list(
            m["variable"], ["#F1F5F9", m["color"]]).copy()
        cmap.set_bad("#ECEFF2")
        img = ax.imshow(np.ma.masked_invalid(z), cmap=cmap, vmin=0,
                        vmax=max(float(vmax), 0.01), aspect="auto",
                        interpolation="nearest")
        ax.set_xticks(range(4), ["S01\nStationary", "S02\nStationary",
                                 "NS01\nMotion", "NS02\nMotion"])
        ax.set_yticks(range(8), PARTICIPANTS)
        ax.tick_params(length=0)
        ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 8, 1), minor=True)
        ax.grid(which="minor", color="white", lw=2.0)
        ax.tick_params(which="minor", length=0)
        ax.set_xlabel("Recording run")
        ax.set_ylabel("Participant")
        ax.text(-.11, 1.035, f"({letter})", transform=ax.transAxes,
                weight="bold", fontsize=13, ha="left", va="bottom")
        for i in range(8):
            for j in range(4):
                val = z[i, j]
                color = "white" if np.isfinite(val) and val > vmax * .57 else INK
                ax.text(j, i, f"{val:.3f}" if np.isfinite(val) else "N/A",
                        ha="center", va="center", fontsize=9.2,
                        color=color if np.isfinite(val) else GREY, weight="bold")
        cb = fig.colorbar(img, ax=ax, fraction=.046, pad=.027)
        cb.ax.set_ylabel(m["ylabel"], fontsize=10, weight="bold")
        cb.ax.tick_params(labelsize=9)
    fig.subplots_adjust(left=.095, right=.96, bottom=.135, top=.91, wspace=.35)
    save(fig, out, dpi)


def paired_figure(contrasts, out, dpi):
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.9))
    for ax, m, letter in zip(axes, MEASURES, "ab"):
        c = contrasts[contrasts.variable.eq(m["variable"])].set_index("participant").reindex(PARTICIPANTS)
        values = np.r_[c.stationary_RMSE.to_numpy(float), c.motion_RMSE.to_numpy(float)]
        peak = max(1., np.nanmax(values) if np.isfinite(values).any() else 1.)
        for _, item in c.iterrows():
            if np.isfinite(item.stationary_RMSE) and np.isfinite(item.motion_RMSE):
                ax.plot([0, 1], [item.stationary_RMSE, item.motion_RMSE],
                        lw=1.2, color=GREY, alpha=.55, zorder=2)
            if np.isfinite(item.stationary_RMSE):
                ax.scatter(0, item.stationary_RMSE, color=BLUE, s=55,
                           edgecolor="white", lw=.8, zorder=3)
            if np.isfinite(item.motion_RMSE):
                ax.scatter(1, item.motion_RMSE, color=ORANGE, s=55,
                           edgecolor="white", lw=.8, zorder=3)
        ax.set(xlim=(-.32, 1.33), ylim=(0, peak * 1.16))
        ax.set_xticks([0, 1], ["Stationary", "Controlled motion"])
        ax.set_ylabel(m["ylabel"])
        ax.text(.01, .97, f"({letter})", transform=ax.transAxes,
                fontsize=12, va="top", weight="bold")
        ax.yaxis.grid(True, color=GRID, lw=.6, alpha=.55)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=.105, right=.98, bottom=.14, top=.975, wspace=.31)
    save(fig, out, dpi)


def main(argv=None):
    a = parse_args(argv)
    if a.dpi < 72:
        raise ValueError("DPI must be >= 72")
    if not np.isfinite(a.tolerance) or a.tolerance <= 0:
        raise ValueError("Pairing tolerance must be positive")
    if not np.isfinite(a.spo2_offset):
        raise ValueError("SpO2 M1 offset must be finite")

    here = Path(__file__).resolve().parent
    root = (a.data_dir or here).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Data root does not exist: {root}")

    out = (a.out_dir or (root / "Run_Level_R2_Standalone_Output")).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    print("HEMOX R2 standalone run-level analysis")
    print("Data root:", root)
    print("SpO2 model: M1 = native fused SpO2 -", a.spo2_offset, "pp")
    print("HR model: M0 native fused HR")
    print("Pairing tolerance: +/-", a.tolerance, "s")
    print("Font:", setup_style(a.dpi))

    # Build directly from raw participant exports.
    try:
        pairs, qa = build_pairs(root, a.tolerance, a.spo2_offset)
    except Exception:
        # Write any possible clue about the intended folder structure.
        raise

    qa.to_csv(out / "source_run_QA.csv", index=False)
    pairs.to_csv(out / "paired_observations_self_calculated.csv", index=False,
                 float_format="%.8f")

    run_df, participant_df, condition_df, pooled_df, contrasts = calculate(pairs)
    run_df.to_csv(out / "run_level_bias_MAE_RMSE.csv", index=False, float_format="%.8f")
    participant_df.to_csv(out / "participant_level_bias_MAE_RMSE.csv", index=False, float_format="%.8f")
    condition_df.to_csv(out / "participant_stationary_vs_motion_metrics.csv", index=False, float_format="%.8f")
    pooled_df.to_csv(out / "pooled_heldout_descriptive_metrics.csv", index=False, float_format="%.8f")
    contrasts.to_csv(out / "stationary_vs_motion_paired_RMSE_differences.csv", index=False, float_format="%.8f")

    heatmap(run_df, out / f"Heldout_RunLevel_RMSE_Heatmap.{a.format}", a.dpi)
    paired_figure(contrasts, out / f"Heldout_Stationary_vs_Motion_Paired_RMSE.{a.format}", a.dpi)

    print("\nPooled P03--P10 metrics:")
    print(pooled_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\nExpected manuscript check (same source/alignment/model):")
    expected = {"SpO2": 0.968, "HR": 2.495}
    for row in pooled_df.itertuples(index=False):
        target = expected[row.variable]
        delta = float(row.RMSE) - target
        tag = "OK" if abs(delta) <= 0.01 else "CHECK SOURCE/PAIRING"
        print(f"  {row.variable}: RMSE={row.RMSE:.3f}; manuscript={target:.3f}; "
              f"difference={delta:+.3f} [{tag}]")

    readme = f"""HEMOX R2 standalone run-level analysis
Data root: {root}
Participants: P03-P10
Runs: S01, S02, NS01, NS02
SpO2 model: M1 = native fused SpO2 - {a.spo2_offset:.6f} pp, clipped 0-100
HR model: M0 native fused HR
HEMOX inclusion: quality_valid=true
Pairing: SP-20-first, nearest unique HEMOX elapsed-time observation, +/-{a.tolerance:g} s
RMSE: sqrt(mean(error^2)) over metric-valid paired observations
Pooled RMSE is observation-weighted and is NOT the arithmetic mean of heatmap cells.
No precomputed native_and_selected_reference_pairs.csv is required.
"""
    (out / "README_standalone.txt").write_text(readme, encoding="utf-8")
    print("\nOutput:", out)
    return out


if __name__ == "__main__":
    main()
