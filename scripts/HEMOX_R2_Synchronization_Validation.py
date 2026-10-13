from __future__ import annotations

import argparse
import io
import re
import sys
import zipfile
import textwrap
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd


def configure_console_encoding():
    """Use UTF-8 for Windows/VS Code console output when supported.

    The analysis itself already writes CSV/TXT files as Unicode.  This only
    prevents Windows cp1252 terminals from crashing on symbols such as Δ.
    """
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except (AttributeError, ValueError, OSError):
                pass


def console_safe(value):
    """Fallback-safe text for consoles whose encoding cannot represent Unicode."""
    value = str(value)
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    try:
        value.encode(enc)
        return value
    except UnicodeEncodeError:
        return (value.replace("Δ", "Delta")
                     .replace("±", "+/-")
                     .replace("–", "-"))


# ============================= USER SETTINGS =============================

DATA_ROOT = r""
PAIRS_CSV = r""
OFFSETS_CSV = r""
OUTPUT_FOLDER = r""

COHORT = "heldout"               # "heldout" -> P03-P10; "all" -> P01-P10
PAIRING_TOLERANCE_S = 1.0
NOMINAL_FS_HZ = 50.0
RUN_DURATION_S = 300.0

IMAGE_FORMAT = "png"             # png, jpg, svg
DPI = 300                         # 300 is publication-safe and Windows-memory-safe

# R2 plotting palette
BLUE = "#0072B2"
RED = "#C62828"
GREY = "#777F87"
ORANGE = "#D55E00"
GREEN = "#009E73"
INK = "#253348"
GRID = "#D7DCE2"
BORDER = "#AAB4BF"
LIGHT = "#F5F7F9"

RUNS = ("S01", "S02", "NS01", "NS02")
ALL_PARTICIPANTS = tuple(f"P{i:02d}" for i in range(1, 11))
HELDOUT_PARTICIPANTS = tuple(f"P{i:02d}" for i in range(3, 11))

RAW_REQUIRED = (
    "counter", "common_time_us", "sync_delta_us",
    "ch1_ir_raw", "ch1_red_raw", "ch2_ir_raw", "ch2_red_raw",
)

# =======================================================================


@dataclass(frozen=True)
class Source:
    kind: str                 # "directory" or "zip"
    path: Path


@dataclass(frozen=True)
class FileRef:
    source: Source
    item: str                 # absolute path for directory, member name for zip


# ------------------------------- CLI ------------------------------------

def parser(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--data-root", type=Path,
                   help="Folder containing P01...P10 folders or participant ZIPs")
    p.add_argument("--pairs-csv", type=Path,
                   help="Existing native_and_selected_reference_pairs.csv; used for offsets/provenance only")
    p.add_argument("--offsets-csv", type=Path,
                   help="Optional verified per-run offset CSV; overrides offsets found in pairs CSV")
    p.add_argument("--out-dir", type=Path)
    p.add_argument("--cohort", choices=("heldout", "all"), default=COHORT)
    p.add_argument("--pairing-tolerance", type=float, default=PAIRING_TOLERANCE_S)
    p.add_argument("--fs", type=float, default=NOMINAL_FS_HZ)
    p.add_argument("--duration", type=float, default=RUN_DURATION_S)
    p.add_argument("--format", choices=("png", "jpg", "svg"), default=IMAGE_FORMAT)
    p.add_argument("--dpi", type=int, default=DPI)
    p.add_argument("--allow-incomplete", action="store_true",
                   help="Allow missing participant/runs; outputs will be marked incomplete")
    return p.parse_args(argv)


# ------------------------------ Styling ---------------------------------

def setup_style():
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next(
        (x for x in ("Arial Narrow", "Arial", "Liberation Sans", "DejaVu Sans")
         if x in installed),
        "DejaVu Sans"
    )
    plt.rcParams.update({
        "font.family": font,
        "font.weight": "bold",
        "font.size": 11.5,
        "axes.labelweight": "bold",
        "axes.labelsize": 13.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 10.5,
        "axes.linewidth": 1.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })
    return font


def save_figure(fig, path: Path, dpi: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".svg":
        fig.savefig(path, facecolor="white")
    elif suffix in (".jpg", ".jpeg"):
        fig.savefig(path, dpi=dpi, facecolor="white",
                    pil_kwargs={"quality": 96, "subsampling": 0})
    else:
        fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


# -------------------------- Naming / discovery ---------------------------

def normalize_subject(value) -> Optional[str]:
    m = re.search(r"P\s*0*(\d{1,2})", str(value), re.I)
    if not m:
        return None
    n = int(m.group(1))
    return f"P{n:02d}" if 1 <= n <= 10 else None


def detect_run(text: str, participant: str) -> Optional[str]:
    u = str(text).replace("\\", "/").upper()
    n = int(participant[1:])

    # Most HEMOX names contain strings such as HEMOXP03S01 or P03NS02.
    m = re.search(rf"P0*{n}[_-]?(NS0[12]|S0[12])", u)
    if m:
        return m.group(1)

    # Fallback for files nested inside an already participant-specific folder.
    m = re.search(r"(?<![A-Z0-9])(NS0[12]|S0[12])(?=$|[^A-Z0-9])", u)
    return m.group(1) if m else None


def choose_root(arg: Optional[Path]) -> Path:
    if arg is not None:
        root = arg.expanduser().resolve()
    elif DATA_ROOT.strip():
        root = Path(DATA_ROOT).expanduser().resolve()
    else:
        root = Path(__file__).resolve().parent
    if not root.exists():
        raise FileNotFoundError(f"Data root does not exist: {root}")
    return root


def choose_participant_sources(root: Path, participants) -> dict[str, Source]:
    """Prefer extracted participant folders; otherwise use one participant ZIP."""
    result = {}
    for participant in participants:
        if root.is_dir():
            dirs = [p for p in root.iterdir()
                    if p.is_dir() and p.name.upper() == participant]
            if len(dirs) > 1:
                raise ValueError(f"Conflicting directories for {participant}: {dirs}")
            if dirs:
                result[participant] = Source("directory", dirs[0])
                continue

            zips = [p for p in root.iterdir()
                    if p.is_file() and p.suffix.lower() == ".zip"
                    and re.fullmatch(re.escape(participant) + r"(?:\(\d+\))?", p.stem, re.I)]
            exact = [p for p in zips if p.stem.upper() == participant]
            if len(exact) == 1:
                result[participant] = Source("zip", exact[0])
            elif len(zips) == 1:
                result[participant] = Source("zip", zips[0])
            elif len(zips) > 1:
                raise ValueError(f"Ambiguous ZIP archives for {participant}: {zips}")
        elif root.is_file() and root.suffix.lower() == ".zip":
            # A single archive may contain one participant only.
            result[participant] = Source("zip", root)
    return result


def list_csv_items(source: Source):
    if source.kind == "directory":
        return [str(p.resolve()) for p in source.path.rglob("*.csv") if p.is_file()]
    with zipfile.ZipFile(source.path) as z:
        return [n for n in z.namelist() if n.lower().endswith(".csv")]


def read_bytes(ref: FileRef) -> bytes:
    if ref.source.kind == "directory":
        return Path(ref.item).read_bytes()
    with zipfile.ZipFile(ref.source.path) as z:
        return z.read(ref.item)


def index_participant(source: Source, participant: str):
    """Return RAW/FUSED/SP20 file references for each run."""
    out = {run: {"raw": [], "fused": [], "sp20": []} for run in RUNS}
    for item in list_csv_items(source):
        run = detect_run(item, participant)
        if run not in RUNS:
            continue
        base = Path(item).name.upper()
        ref = FileRef(source, item)
        if base.endswith("_RAW.CSV"):
            out[run]["raw"].append(ref)
        if base.endswith("_FUSED_TRENDS.CSV") and "CALIBRATED" not in base:
            out[run]["fused"].append(ref)
        if re.match(r"SP[-_ ]?20", base):
            out[run]["sp20"].append(ref)
    return out


def require_unique(refs, label, participant, run, allow_incomplete):
    if len(refs) == 1:
        return refs[0]
    if len(refs) == 0 and allow_incomplete:
        return None
    if len(refs) == 0:
        raise FileNotFoundError(f"{participant}/{run}: missing {label} CSV")
    raise ValueError(
        f"{participant}/{run}: found {len(refs)} {label} CSVs. "
        "Remove duplicate exports or use a clean participant folder."
    )


# ----------------------------- Input parsers -----------------------------

def decode_text(payload: bytes) -> str:
    if payload.startswith((b"\xff\xfe", b"\xfe\xff")) or b"\x00" in payload[:200]:
        return payload.decode("utf-16")
    try:
        return payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        return payload.decode("cp1252")


def read_raw(ref: FileRef) -> pd.DataFrame:
    df = pd.read_csv(io.BytesIO(read_bytes(ref)), low_memory=False)
    missing = [c for c in RAW_REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"RAW CSV lacks required columns: {missing}")
    return df


def read_fused(ref: FileRef) -> pd.DataFrame:
    df = pd.read_csv(io.BytesIO(read_bytes(ref)), low_memory=False)
    if "elapsed_s" not in df.columns:
        raise ValueError("FUSED_TRENDS lacks elapsed_s")
    t = pd.to_numeric(df["elapsed_s"], errors="coerce")
    if t.isna().any() or t.duplicated().any() or not t.is_monotonic_increasing:
        raise ValueError("FUSED_TRENDS elapsed_s must be finite, unique and monotonic")
    return pd.DataFrame({"elapsed_s": t.to_numpy(float)})


def parse_sp20_time(value: str):
    value = " ".join(str(value).strip().split())
    formats = (
        "%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %H:%M:%S",
        "%d/%m/%Y %I:%M:%S %p", "%d/%m/%Y %H:%M:%S",
        "%Y-%m-%d %I:%M:%S %p", "%Y-%m-%d %H:%M:%S",
    )
    for fmt in formats:
        try:
            return pd.Timestamp(datetime.strptime(value, fmt))
        except ValueError:
            pass
    return pd.to_datetime(value, errors="coerce")


def read_sp20(ref: FileRef) -> pd.DataFrame:
    text = decode_text(read_bytes(ref))
    lines = text.splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        normalized = line.replace("#", "")
        if "Time" in normalized and "SpO2(%)" in normalized and "PR(bpm)" in normalized:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("SP-20 CSV header containing Time/SpO2(%)/PR(bpm) not found")

    df = pd.read_csv(io.StringIO("\n".join(lines[header_idx:])),
                     skipinitialspace=True, dtype=str)
    df.columns = [str(c).strip().lstrip("#").strip() for c in df.columns]
    needed = {"Time", "SpO2(%)", "PR(bpm)"}
    if not needed.issubset(df.columns):
        raise ValueError(f"SP-20 CSV missing {needed - set(df.columns)}")

    parsed = pd.Series([parse_sp20_time(v) for v in df["Time"]])
    good = parsed.notna()
    df = df.loc[good].reset_index(drop=True)
    parsed = parsed.loc[good].reset_index(drop=True)
    if len(parsed) < 2:
        raise ValueError("SP-20 CSV has fewer than two valid timestamp rows")
    if not parsed.is_monotonic_increasing or parsed.duplicated().any():
        raise ValueError("SP-20 timestamps are non-monotonic or duplicated")

    elapsed = (parsed - parsed.iloc[0]).dt.total_seconds().to_numpy(float)
    return pd.DataFrame({"elapsed_s": elapsed, "timestamp": parsed.astype(str)})


# ------------------------- Alignment/provenance --------------------------

def locate_pairs_csv(arg: Optional[Path], root: Path) -> Optional[Path]:
    if arg is not None:
        p = arg.expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(p)
        return p
    if PAIRS_CSV.strip():
        p = Path(PAIRS_CSV).expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(p)
        return p

    candidates = []
    if root.is_dir():
        direct = root / "Study_HEMOX_Correlation_outputs" / "native_and_selected_reference_pairs.csv"
        if direct.is_file():
            candidates.append(direct)
        candidates += [p for p in root.rglob("native_and_selected_reference_pairs.csv") if p.is_file()]
    unique = []
    seen = set()
    for p in candidates:
        rp = p.resolve()
        if rp not in seen:
            unique.append(rp)
            seen.add(rp)
    if not unique:
        return None
    if len(unique) == 1:
        return unique[0]
    # Prefer canonical output folder.
    canonical = [p for p in unique if p.parent.name == "Study_HEMOX_Correlation_outputs"]
    return canonical[0] if len(canonical) == 1 else unique[0]


def offsets_from_pairs(path: Optional[Path], participants):
    offsets = {}
    accepted_counts = {}
    methods = {}
    if path is None:
        return offsets, accepted_counts, methods

    df = pd.read_csv(path, low_memory=False)
    if not {"participant", "run"}.issubset(df.columns):
        return offsets, accepted_counts, methods
    df["participant"] = df["participant"].map(normalize_subject)
    df["run"] = df["run"].astype(str).str.upper().str.strip()
    df = df[df.participant.isin(participants) & df.run.isin(RUNS)].copy()

    for (p, r), g in df.groupby(["participant", "run"], dropna=False):
        accepted_counts[(p, r)] = len(g)
        if "offset_seconds" in g.columns:
            x = pd.to_numeric(g["offset_seconds"], errors="coerce").dropna().to_numpy(float)
            if len(x):
                if np.ptp(x) > 1e-9:
                    raise ValueError(f"{p}/{r}: pairs CSV contains inconsistent offset_seconds")
                offsets[(p, r)] = float(np.median(x))
        if "alignment_method" in g.columns:
            vals = [str(v).strip() for v in g["alignment_method"].dropna().unique()]
            methods[(p, r)] = vals[0] if len(vals) == 1 else ";".join(vals)
    return offsets, accepted_counts, methods


def offsets_from_csv(path: Optional[Path], participants):
    if path is None and OFFSETS_CSV.strip():
        path = Path(OFFSETS_CSV).expanduser().resolve()
    if path is None:
        return {}, {}, {}
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    needed = {"participant", "run", "offset_seconds"}
    if not needed.issubset(df.columns):
        raise ValueError(f"Offsets CSV must contain {needed}")
    offsets, methods, evidence = {}, {}, {}
    for row in df.itertuples(index=False):
        p = normalize_subject(row.participant)
        r = str(row.run).upper().strip()
        if p not in participants or r not in RUNS:
            continue
        key = (p, r)
        if key in offsets:
            raise ValueError(f"Duplicate offset row for {p}/{r}")
        off = float(row.offset_seconds)
        if not np.isfinite(off):
            raise ValueError(f"Nonfinite offset for {p}/{r}")
        offsets[key] = off
        methods[key] = str(getattr(row, "method", "verified_offset_csv")).strip()
        evidence[key] = str(getattr(row, "evidence", "")).strip()
    return offsets, methods, evidence


# --------------------------- RAW calculations ----------------------------

def finite_numeric(df, col):
    return pd.to_numeric(df[col], errors="coerce").to_numpy(float)


def analyze_raw(df: pd.DataFrame, fs: float, duration_s: float):
    nrows = len(df)

    # 1) Direct exported synchronization residual.
    sync_us = finite_numeric(df, "sync_delta_us")
    sync_us = sync_us[np.isfinite(sync_us)]
    abs_sync_us = np.abs(sync_us)

    # 2) Synchronized-record counter integrity.
    counter = finite_numeric(df, "counter")
    counter = counter[np.isfinite(counter)]
    counter_integer = len(counter) and np.max(np.abs(counter - np.rint(counter))) <= 1e-6

    if len(counter) >= 2 and counter_integer:
        c = np.rint(counter).astype(np.int64)
        dc = np.diff(c)
        missing = int(np.sum(np.clip(dc - 1, 0, None)))
        duplicates = int(np.sum(dc == 0))
        backwards = int(np.sum(dc < 0))
        counter_loss_pct = 100.0 * missing / (len(c) + missing) if len(c) + missing else np.nan
        counter_start = int(c[0])
        counter_end = int(c[-1])
    else:
        missing = duplicates = backwards = np.nan
        counter_loss_pct = np.nan
        counter_start = counter_end = np.nan

    # 3) Bilateral optical-row completeness.
    bilateral_mask = np.ones(nrows, dtype=bool)
    for col in ("ch1_ir_raw", "ch1_red_raw", "ch2_ir_raw", "ch2_red_raw"):
        bilateral_mask &= np.isfinite(finite_numeric(df, col))
    bilateral_complete = int(np.sum(bilateral_mask))
    bilateral_incomplete = int(nrows - bilateral_complete)
    bilateral_pct = 100.0 * bilateral_complete / nrows if nrows else np.nan

    # 4) Common-timeline integrity.
    t_us = finite_numeric(df, "common_time_us")
    t_us = t_us[np.isfinite(t_us)]
    nominal_us = 1e6 / fs
    if len(t_us) >= 2:
        dt_us = np.diff(t_us)
        jitter_us = np.abs(dt_us - nominal_us)
        positive = dt_us > 0
        gap_intervals = np.where(
            dt_us > 1.5 * nominal_us,
            np.maximum(np.rint(dt_us / nominal_us).astype(int) - 1, 0),
            0,
        )
        timeline_gap_estimate = int(np.sum(gap_intervals))
        nonpositive_dt = int(np.sum(~positive))
        common_duration_s = float((t_us[-1] - t_us[0]) / 1e6)
        median_dt_ms = float(np.median(dt_us) / 1000.0)
        p95_jitter_ms = float(np.percentile(jitter_us, 95) / 1000.0)
        max_jitter_ms = float(np.max(jitter_us) / 1000.0)
    else:
        timeline_gap_estimate = nonpositive_dt = np.nan
        common_duration_s = median_dt_ms = p95_jitter_ms = max_jitter_ms = np.nan

    expected_rows = int(round(fs * duration_s))
    row_coverage_pct = 100.0 * nrows / expected_rows if expected_rows else np.nan

    metrics = {
        "raw_rows": nrows,
        "expected_rows_nominal": expected_rows,
        "row_coverage_pct": row_coverage_pct,
        "sync_n": len(sync_us),
        "sync_median_signed_us": float(np.median(sync_us)) if len(sync_us) else np.nan,
        "sync_median_abs_us": float(np.median(abs_sync_us)) if len(abs_sync_us) else np.nan,
        "sync_mean_abs_us": float(np.mean(abs_sync_us)) if len(abs_sync_us) else np.nan,
        "sync_p95_abs_us": float(np.percentile(abs_sync_us, 95)) if len(abs_sync_us) else np.nan,
        "sync_max_abs_us": float(np.max(abs_sync_us)) if len(abs_sync_us) else np.nan,
        "counter_start": counter_start,
        "counter_end": counter_end,
        "counter_missing_increments": missing,
        "counter_duplicates": duplicates,
        "counter_backwards": backwards,
        "counter_loss_pct": counter_loss_pct,
        "bilateral_complete_rows": bilateral_complete,
        "bilateral_incomplete_rows": bilateral_incomplete,
        "bilateral_completeness_pct": bilateral_pct,
        "common_duration_s": common_duration_s,
        "median_dt_ms": median_dt_ms,
        "p95_abs_timeline_jitter_ms": p95_jitter_ms,
        "max_abs_timeline_jitter_ms": max_jitter_ms,
        "timeline_gap_estimate": timeline_gap_estimate,
        "nonpositive_common_dt": nonpositive_dt,
    }
    return metrics, sync_us


# ---------------------- SP20 pre-tolerance analysis ----------------------

def nearest_pre_tolerance(sp: pd.DataFrame, fused: pd.DataFrame,
                          offset_s: float, tolerance_s: float,
                          participant: str, run: str):
    """
    For every SP-20 timestamp inside the common overlap, find the nearest aligned
    HEMOX FUSED_TRENDS timestamp WITHOUT applying the tolerance and WITHOUT
    excluding a residual because it exceeds the tolerance.

    This distribution is therefore suitable for validating the tolerance itself.
    It is intentionally different from the final one-to-one accepted-pair table.
    """
    s = sp["elapsed_s"].to_numpy(float)
    h_native = fused["elapsed_s"].to_numpy(float)
    h = h_native + float(offset_s)
    if len(s) == 0 or len(h) == 0:
        return pd.DataFrame()

    low = max(float(s[0]), float(h[0]))
    high = min(float(s[-1]), float(h[-1]))
    if high < low:
        return pd.DataFrame()

    rows = []
    for ref_idx, x in enumerate(s):
        if x < low or x > high:
            continue
        at = int(np.searchsorted(h, x))
        candidates = []
        if at < len(h):
            candidates.append(at)
        if at > 0:
            candidates.append(at - 1)
        if not candidates:
            continue
        j = min(candidates, key=lambda k: (abs(h[k] - x), k))
        signed = float(h[j] - x)
        rows.append({
            "participant": participant,
            "run": run,
            "sp20_index": ref_idx,
            "sp20_elapsed_s": float(x),
            "nearest_hemox_index": j,
            "hemox_elapsed_native_s": float(h_native[j]),
            "hemox_elapsed_aligned_s": float(h[j]),
            "offset_seconds": float(offset_s),
            "residual_s": signed,
            "abs_residual_s": abs(signed),
            "within_tolerance": abs(signed) <= tolerance_s + 1e-12,
        })
    return pd.DataFrame(rows)


# ----------------------------- Summaries ---------------------------------

def pooled_summary(run_df: pd.DataFrame, residual_df: pd.DataFrame,
                   sync_pool: list[np.ndarray], tolerance_s: float,
                   accepted_counts: dict, fs: float):
    sync = np.concatenate([x for x in sync_pool if len(x)]) if any(len(x) for x in sync_pool) else np.array([])
    abs_sync = np.abs(sync)

    total_rows = int(pd.to_numeric(run_df["raw_rows"], errors="coerce").fillna(0).sum())
    total_expected = int(pd.to_numeric(run_df["expected_rows_nominal"], errors="coerce").fillna(0).sum())
    total_missing = int(pd.to_numeric(run_df["counter_missing_increments"], errors="coerce").fillna(0).sum())
    total_duplicates = int(pd.to_numeric(run_df["counter_duplicates"], errors="coerce").fillna(0).sum())
    total_backwards = int(pd.to_numeric(run_df["counter_backwards"], errors="coerce").fillna(0).sum())
    total_bilateral = int(pd.to_numeric(run_df["bilateral_complete_rows"], errors="coerce").fillna(0).sum())
    total_incomplete = int(pd.to_numeric(run_df["bilateral_incomplete_rows"], errors="coerce").fillna(0).sum())
    total_timeline_gaps = int(pd.to_numeric(run_df["timeline_gap_estimate"], errors="coerce").fillna(0).sum())

    residual = residual_df["abs_residual_s"].to_numpy(float) if len(residual_df) else np.array([])
    within = residual <= tolerance_s if len(residual) else np.array([], dtype=bool)

    rows = []
    rows.append({
        "metric": "CH1-CH2 synchronization residual",
        "symbol_or_definition": "sync_delta_us; |Δt_CH1-CH2|",
        "n": len(sync),
        "result": (
            f"median |Δt|={np.median(abs_sync)/1000:.3f} ms; "
            f"p95={np.percentile(abs_sync,95)/1000:.3f} ms; "
            f"max={np.max(abs_sync)/1000:.3f} ms"
            if len(abs_sync) else "N/A"
        ),
        "criterion_or_context": f"Nominal 50-Hz sample interval={1000/fs:.1f} ms",
    })
    rows.append({
        "metric": "Synchronized-record counter integrity",
        "symbol_or_definition": "counter increments",
        "n": total_rows,
        "result": (
            f"missing increments={total_missing}; duplicates={total_duplicates}; "
            f"backward jumps={total_backwards}; nominal rows={total_rows}/{total_expected}"
        ),
        "criterion_or_context": "Ideal: 0 missing, 0 duplicates, 0 backward jumps",
    })
    rows.append({
        "metric": "Bilateral optical-row completeness",
        "symbol_or_definition": "finite CH1 RED/IR and CH2 RED/IR in same synchronized row",
        "n": total_rows,
        "result": (
            f"{total_bilateral}/{total_rows} complete "
            f"({100*total_bilateral/total_rows:.4f}%); incomplete={total_incomplete}"
            if total_rows else "N/A"
        ),
        "criterion_or_context": "Directly verifies bilateral content in exported synchronized records",
    })
    rows.append({
        "metric": "Common-timeline integrity",
        "symbol_or_definition": "common_time_us spacing",
        "n": int(pd.to_numeric(run_df["raw_rows"], errors="coerce").notna().sum()),
        "result": (
            f"median run p95 jitter={pd.to_numeric(run_df['p95_abs_timeline_jitter_ms'],errors='coerce').median():.3f} ms; "
            f"worst run p95={pd.to_numeric(run_df['p95_abs_timeline_jitter_ms'],errors='coerce').max():.3f} ms; "
            f"inferred missing intervals={total_timeline_gaps}"
        ),
        "criterion_or_context": f"Nominal interval={1000/fs:.1f} ms; counter result determines whether a timing gap is packet loss",
    })
    rows.append({
        "metric": "HEMOX-SP20 pre-tolerance nearest-timestamp residual",
        "symbol_or_definition": "min_j |(t_HEMOX,j + offset) - t_SP20| before tolerance gate",
        "n": len(residual),
        "result": (
            f"median={np.median(residual):.3f} s; p95={np.percentile(residual,95):.3f} s; "
            f"max={np.max(residual):.3f} s; within ±{tolerance_s:g} s={100*np.mean(within):.2f}%"
            if len(residual) else "N/A"
        ),
        "criterion_or_context": f"Configured pairing tolerance=±{tolerance_s:g} s; computed BEFORE gating",
    })
    if accepted_counts:
        rows.append({
            "metric": "Previously accepted post-gate pairs",
            "symbol_or_definition": "native_and_selected_reference_pairs.csv rows",
            "n": int(sum(accepted_counts.values())),
            "result": f"accepted pair rows={int(sum(accepted_counts.values()))}",
            "criterion_or_context": "Provenance only; not used to validate the tolerance",
        })
    return pd.DataFrame(rows)


def latex_escape(value):
    return (str(value).replace("\\", r"\textbackslash{}")
            .replace("&", r"\&").replace("%", r"\%")
            .replace("_", r"\_").replace("#", r"\#"))


def write_latex_table(summary: pd.DataFrame, path: Path):
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Quantitative synchronization, stream-integrity, and reference-pairing validation.}",
        r"\label{tab:synchronization_validation}",
        r"\begin{tabular}{p{0.22\textwidth}p{0.24\textwidth}p{0.32\textwidth}p{0.17\textwidth}}",
        r"\hline",
        r"\textbf{Metric} & \textbf{Definition} & \textbf{Measured result} & \textbf{Criterion / context} \\",
        r"\hline",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"{latex_escape(row.metric)} & {latex_escape(row.symbol_or_definition)} & "
            f"{latex_escape(row.result)} & {latex_escape(row.criterion_or_context)} \\\\"
        )
    lines += [r"\hline", r"\end{tabular}", r"\end{table*}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


# ------------------------------- Plots -----------------------------------

def run_axis(run_df, participants):
    order = {(p, r): i for i, p in enumerate(participants) for r in RUNS}
    # correct order index manually
    order = {}
    k = 0
    for p in participants:
        for r in RUNS:
            order[(p, r)] = k
            k += 1
    out = run_df.copy()
    out["_x"] = [order.get((p, r), np.nan) for p, r in zip(out.participant, out.run)]
    out = out.sort_values("_x").reset_index(drop=True)
    centers = [np.mean([order[(p, r)] for r in RUNS]) for p in participants]
    return out, centers


def decorate_run_axis(ax, centers, participants, n_positions):
    ax.set_xticks(centers, participants, rotation=0)
    for i in range(1, len(participants)):
        ax.axvline(i * len(RUNS) - 0.5, color=GRID, lw=0.6, alpha=0.7, zorder=0)
    ax.set_xlim(-0.6, n_positions - 0.4)


def plot_three_panel(run_df: pd.DataFrame, residual_df: pd.DataFrame,
                     participants, tolerance_s: float, fs: float,
                     out_path: Path, dpi: int):
    rd, centers = run_axis(run_df, participants)
    x = rd["_x"].to_numpy(float)

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.7))

    # (a) Direct sync_delta_us
    ax = axes[0]
    p95_ms = pd.to_numeric(rd["sync_p95_abs_us"], errors="coerce").to_numpy(float) / 1000.0
    max_ms = pd.to_numeric(rd["sync_max_abs_us"], errors="coerce").to_numpy(float) / 1000.0
    good = np.isfinite(p95_ms)
    ax.scatter(x[good], p95_ms[good], marker="o", s=30, color=BLUE,
               edgecolor="white", linewidth=0.45, zorder=3,
               label=r"Run p95 $|\Delta t|$")
    good_max = np.isfinite(max_ms)
    ax.scatter(x[good_max], max_ms[good_max], marker="^", s=30, color=ORANGE,
               edgecolor="white", linewidth=0.4, zorder=4,
               label=r"Run max $|\Delta t|$")
    # lightly connect p95 to max as a whisker
    for xi, a, b in zip(x, p95_ms, max_ms):
        if np.isfinite(a) and np.isfinite(b):
            ax.vlines(xi, a, b, color=GREY, lw=0.65, alpha=0.45, zorder=1)

    all_p95 = p95_ms[np.isfinite(p95_ms)]
    if len(all_p95) and np.nanmax(max_ms) == 0:
        ax.set_ylim(-0.01, 0.10)
    else:
        ymax = np.nanmax(max_ms) if np.isfinite(max_ms).any() else 1.0
        ax.set_ylim(0, max(0.1, 1.12 * ymax))
    ax.set_ylabel(r"CH$_1$–CH$_2$ sync residual (ms)")
    ax.set_xlabel("Held-out participant")
    decorate_run_axis(ax, centers, participants, len(rd))
    ax.grid(True, axis="y", color=GRID, lw=0.55, alpha=0.65)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.6, loc="upper right")
    ax.text(-.075, 1.015, "(a)", transform=ax.transAxes, ha="left", va="top",
            fontsize=15, weight="bold", clip_on=False)

    # (b) Counter integrity + common timeline gap distinction
    ax = axes[1]
    missing = pd.to_numeric(rd["counter_missing_increments"], errors="coerce").fillna(0).to_numpy(float)
    dup = pd.to_numeric(rd["counter_duplicates"], errors="coerce").fillna(0).to_numpy(float)
    back = pd.to_numeric(rd["counter_backwards"], errors="coerce").fillna(0).to_numpy(float)
    time_gaps = pd.to_numeric(rd["timeline_gap_estimate"], errors="coerce").fillna(0).to_numpy(float)

    width = 0.72
    ax.bar(x, missing, width=width, color=BLUE, alpha=0.90, label="Missing counter increments")
    ax.bar(x, dup, bottom=missing, width=width, color=ORANGE, alpha=0.82, label="Duplicate counters")
    ax.bar(x, back, bottom=missing + dup, width=width, color=RED, alpha=0.82, label="Backward jumps")
    ax.scatter(x, time_gaps, marker="D", s=24, color=GREEN, zorder=4,
               label=r"Common-time inferred gaps")

    max_anom = float(np.nanmax(np.r_[missing + dup + back, time_gaps])) if len(rd) else 0
    ax.set_ylim(-0.08 if max_anom == 0 else 0,
                1.0 if max_anom == 0 else max(1.0, 1.25 * max_anom))
    ax.set_ylabel("Integrity events (count)")
    ax.set_xlabel("Held-out participant")
    decorate_run_axis(ax, centers, participants, len(rd))
    ax.grid(True, axis="y", color=GRID, lw=0.55, alpha=0.65)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7.6, loc="lower left")
    ax.text(-.075, 1.015, "(b)", transform=ax.transAxes, ha="left", va="top",
            fontsize=15, weight="bold", clip_on=False)

    bilateral = pd.to_numeric(rd["bilateral_completeness_pct"], errors="coerce")
    total_rows = int(pd.to_numeric(rd["raw_rows"], errors="coerce").fillna(0).sum())
    total_missing = int(missing.sum())
    total_dup = int(dup.sum())
    total_back = int(back.sum())
    bilateral_pooled = (
        100.0 * pd.to_numeric(rd["bilateral_complete_rows"], errors="coerce").fillna(0).sum() / total_rows
        if total_rows else np.nan
    )
    ax.text(.98, .96,
            f"Rows = {total_rows:,}\n"
            f"Counter anomalies = {total_missing + total_dup + total_back}\n"
            f"Bilateral completeness = {bilateral_pooled:.4f}%",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.8,
            bbox=dict(boxstyle="round,pad=.32", fc="white", ec=BORDER, alpha=.95))

    # (c) Pre-tolerance nearest residual ECDF
    ax = axes[2]
    if len(residual_df):
        abs_r = np.sort(residual_df["abs_residual_s"].to_numpy(float))
        y = 100.0 * np.arange(1, len(abs_r) + 1) / len(abs_r)
        ax.plot(abs_r, y, color=GREEN, lw=2.0,
                label="Pre-tolerance nearest residual")
        ax.axvline(tolerance_s, color=RED, lw=1.5, ls="--",
                   label=f"±{tolerance_s:g} s tolerance")
        med = float(np.median(abs_r))
        p95 = float(np.percentile(abs_r, 95))
        mx = float(np.max(abs_r))
        within = 100.0 * np.mean(abs_r <= tolerance_s)
        ax.text(.035, .14,
                f"Median = {med:.3f} s\n"
                f"95th = {p95:.3f} s\n"
                f"Max = {mx:.3f} s\n"
                f"Within tolerance = {within:.2f}%",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=9.0,
                bbox=dict(boxstyle="round,pad=.32", fc="white", ec=BORDER, alpha=.95))
        xmax = max(tolerance_s * 1.05, mx * 1.05, 0.05)
        ax.set_xlim(0, xmax)
    else:
        ax.text(.5, .5, "No pre-tolerance residuals available",
                transform=ax.transAxes, ha="center", va="center")
    ax.set_ylim(0, 101)
    ax.set_xlabel(r"Nearest HEMOX–SP-20 residual $|\Delta t|$ (s)")
    ax.set_ylabel("Cumulative reference observations (%)")
    ax.grid(True, color=GRID, lw=0.55, alpha=0.65)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.1, loc="lower right")
    ax.text(.012, .98, "(c)", transform=ax.transAxes, ha="left", va="top",
            fontsize=15, weight="bold")

    fig.subplots_adjust(left=.063, right=.992, bottom=.16, top=.965, wspace=.31)
    save_figure(fig, out_path, dpi)


def plot_summary_table(summary: pd.DataFrame, path: Path, dpi: int):
    display = summary[["metric", "result", "criterion_or_context"]].copy()
    display.columns = ["Metric", "Measured result", "Criterion / context"]

    # Wrap long strings so the rendered table is readable and not clipped.
    widths = [32, 58, 46]
    wrapped_rows = []
    row_line_counts = []
    for values in display.astype(str).itertuples(index=False, name=None):
        wrapped = [textwrap.fill(v, width=w, break_long_words=False)
                   for v, w in zip(values, widths)]
        wrapped_rows.append(wrapped)
        row_line_counts.append(max(x.count("\n") + 1 for x in wrapped))

    fig, ax = plt.subplots(figsize=(14.6, 5.2))
    ax.axis("off")
    table = ax.table(
        cellText=wrapped_rows,
        colLabels=display.columns,
        loc="center",
        cellLoc="left",
        colLoc="left",
        colWidths=[0.23, 0.43, 0.34],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.4)

    # Header + content row heights scaled to wrapped line count.
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(BORDER)
        cell.set_linewidth(0.7)
        cell.PAD = 0.10
        if row == 0:
            cell.set_facecolor("#E9EEF3")
            cell.set_text_props(weight="bold", color=INK, va="center")
            cell.set_height(0.095)
        else:
            cell.set_facecolor("white" if row % 2 else LIGHT)
            cell.set_text_props(color=INK, weight="bold" if col == 0 else "normal",
                                va="center")
            lines = row_line_counts[row - 1]
            cell.set_height(0.070 + 0.028 * max(0, lines - 1))

    fig.subplots_adjust(left=.015, right=.985, top=.98, bottom=.02)
    save_figure(fig, path, dpi)


# -------------------------------- Main -----------------------------------

def main(argv=None):
    a = parser(argv)
    if a.pairing_tolerance <= 0 or a.fs <= 0 or a.duration <= 0:
        raise ValueError("pairing tolerance, fs and duration must be positive")
    if a.dpi < 72:
        raise ValueError("--dpi must be >=72")

    font = setup_style()
    root = choose_root(a.data_root)
    participants = HELDOUT_PARTICIPANTS if a.cohort == "heldout" else ALL_PARTICIPANTS
    sources = choose_participant_sources(root, participants)

    missing_participants = [p for p in participants if p not in sources]
    if missing_participants and not a.allow_incomplete:
        raise FileNotFoundError(
            f"Missing participant sources: {missing_participants}. "
            "Use --allow-incomplete only intentionally."
        )

    pairs_path = locate_pairs_csv(a.pairs_csv, root)
    pair_offsets, accepted_counts, pair_methods = offsets_from_pairs(pairs_path, participants)
    explicit_offsets, explicit_methods, explicit_evidence = offsets_from_csv(a.offsets_csv, participants)

    # Explicit offsets override offsets recovered from the pairs CSV.
    offsets = dict(pair_offsets)
    offsets.update(explicit_offsets)
    methods = dict(pair_methods)
    methods.update(explicit_methods)

    out = (
        a.out_dir
        or (Path(OUTPUT_FOLDER).expanduser() if OUTPUT_FOLDER.strip()
            else ((pairs_path.parent if pairs_path else root) / "Synchronization_Validation_v2_Output"))
    ).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    run_rows = []
    sync_pool = []
    residual_frames = []

    for participant in participants:
        source = sources.get(participant)
        if source is None:
            continue
        indexed = index_participant(source, participant)

        for run in RUNS:
            raw_ref = require_unique(indexed[run]["raw"], "RAW", participant, run, a.allow_incomplete)
            fused_ref = require_unique(indexed[run]["fused"], "FUSED_TRENDS", participant, run, a.allow_incomplete)
            sp_ref = require_unique(indexed[run]["sp20"], "SP-20", participant, run, a.allow_incomplete)

            row = {
                "participant": participant,
                "run": run,
                "source_kind": source.kind,
                "source_path": str(source.path),
                "offset_seconds": float(offsets.get((participant, run), 0.0)),
                "alignment_method": methods.get((participant, run),
                                                "zero_offset_no_verified_offset"),
                "accepted_postgate_pairs": int(accepted_counts.get((participant, run), 0)),
            }

            if raw_ref is None or fused_ref is None or sp_ref is None:
                row["status"] = "INCOMPLETE_SOURCE"
                run_rows.append(row)
                continue