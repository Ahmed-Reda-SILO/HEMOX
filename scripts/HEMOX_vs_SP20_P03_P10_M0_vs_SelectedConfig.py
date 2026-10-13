from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timedelta
from io import BytesIO
from pathlib import Path, PurePosixPath
import re
from zipfile import BadZipFile, ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

# ==================== USER SETTINGS ====================
DATA_FOLDER = r""                    # Root containing P01..P10 folders/ZIP files
PARTICIPANT = "ALL"                   # ALL or P01 ... P10
PARTICIPANT_GROUP = "heldout"        # default for this script
RUN = "ALL"                           # ALL or S01/S02/NS01/NS02
OUTPUT_FOLDER = r""                   # Blank: beside the Python script
SP20_FILE = r""                       # Optional external SP-20 file, single-run mode only
START_S, END_S = 0.0, 300.0
SPO2_Y_LIMITS = None                  # None or (e.g. 90, 100)
HR_Y_LIMITS = None                    # None or (e.g. 55, 110)
HEMOX_OFFSET_S = 0.0                  # Added to exported HEMOX elapsed_s
SP20_START_ELAPSED_S = 0.0            # First SP-20 row's elapsed-time origin
SP20_TIME_SHIFT_S = 0.0               # Additional SP-20 adjustment
SPO2_MODEL = "m1"                     # default R2 held-out overlay uses M1 for SpO2; HR remains M0
FROZEN_M1_OFFSET_PP = 0.906754        # P01/P02-calibrated frozen M1 correction
DPI = 600
IMAGE_FORMAT = "jpg"                  # jpg or png
FIGURE_SET = "none"                   # overlay-only script; per-run figures disabled by default
GENERATE_ALLRUNS_OVERVIEW = True      # global overview across all selected participants/runs
ALLRUNS_OVERVIEW_MODE = "overlay"    # overlay -> one 2-panel figure
SMOOTH_WINDOW_S = 11                  # centered moving-average window for grand-mean traces
# =======================================================

SUBJECTS = tuple(f"P{i:02d}" for i in range(1, 11))
DEVELOPMENT_SUBJECTS = tuple(f"P{i:02d}" for i in range(1, 3))
HELDOUT_SUBJECTS = tuple(f"P{i:02d}" for i in range(3, 11))
RUNS = ("S01", "S02", "NS01", "NS02")
BLUE, RED, GREY = "#0072B2", "#C62828", "#777F87"
GRID, BORDER = "#CFD5DC", "#B0BAC3"


@dataclass(frozen=True)
class Source:
    kind: str                        # 'directory' or 'zip'
    path: Path
    member: str | None = None

    @property
    def filename(self):
        return PurePosixPath((self.member or self.path.name).replace("\\", "/")).name

    @property
    def description(self):
        return str(self.path) + (" :: " + self.member if self.member else "")


def parser(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, default=None, help="Root directory, participant ZIP, or export CSV")
    p.add_argument("--sp20", type=Path, default=None,
                   help="External SP-20 CSV (requires one --participant and one --run)")
    p.add_argument("--participant", default=PARTICIPANT, help="ALL, P01, ..., P10")
    p.add_argument("--participant-group", choices=("all", "development", "heldout"), default=PARTICIPANT_GROUP,
                   help="When --participant is ALL, select all P01-P10, development P01-P02, or heldout P03-P10")
    p.add_argument("--run", default=RUN, help="ALL, S01, S02, NS01, NS02")
    p.add_argument("--start", type=float, default=START_S)
    p.add_argument("--end", type=float, default=END_S)
    p.add_argument("--spo2-model", choices=("m0", "m1"), default=SPO2_MODEL)
    p.add_argument("--hemox-offset", type=float, default=HEMOX_OFFSET_S)
    p.add_argument("--sp20-start-elapsed", type=float, default=SP20_START_ELAPSED_S)
    p.add_argument("--sp20-shift", type=float, default=SP20_TIME_SHIFT_S)
    p.add_argument("--manifest", type=Path, default=None,
                   help="Optional participant,archive CSV selects ZIP revision for each participant")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--format", choices=("jpg", "png"), default=IMAGE_FORMAT)
    p.add_argument("--figure-set", choices=("none", "both", "combined", "singles"), default=FIGURE_SET,
                   help="none recommended here: skip per-run figures and create only the overlay figure")
    p.add_argument("--allruns-overview", choices=("overlay", "none"),
                   default=ALLRUNS_OVERVIEW_MODE if GENERATE_ALLRUNS_OVERVIEW else "none",
                   help="Create global overview across all selected participants/runs")
    p.add_argument("--smooth-window", type=int, default=SMOOTH_WINDOW_S,
                   help="Centered moving-average window in seconds for thick grand-mean traces (default: 11)")
    p.add_argument("--dpi", type=int, default=DPI)
    return p.parse_args(argv)


def normal_subject(value):
    m = re.fullmatch(r"P\s*0*(\d{1,2})", str(value).strip(), flags=re.I)
    return f"P{int(m.group(1)):02d}" if m and 1 <= int(m.group(1)) <= 10 else None


def selected_subjects(value, group="all"):
    if str(value).upper().strip() == "ALL":
        g = str(group).lower().strip()
        if g == "development":
            return DEVELOPMENT_SUBJECTS
        if g == "heldout":
            return HELDOUT_SUBJECTS
        return SUBJECTS
    v = normal_subject(value)
    if v is None: raise ValueError("--participant must be ALL or P01 through P10")
    return (v,)


def selected_runs(value):
    v = str(value).upper().strip()
    if v == "ALL": return RUNS
    if v not in RUNS: raise ValueError("--run must be ALL, S01, S02, NS01, NS02")
    return (v,)


def choose_folder():
    try:
        from tkinter import Tk, filedialog
        root = Tk(); root.withdraw(); root.attributes("-topmost", True)
        folder = filedialog.askdirectory(title="Select root folder with HEMOX P01–P10 exports",
                                         initialdir=str(Path(__file__).resolve().parent))
        root.destroy()
        return Path(folder) if folder else None
    except (ImportError, OSError, RuntimeError):
        return None


def root_folder(args):
    root = args.input or (Path(DATA_FOLDER) if DATA_FOLDER.strip() else None)
    if root is None:
        root = Path(__file__).resolve().parent
        if not any(root.glob("*.zip")) and not any(root.rglob("*_FUSED_TRENDS.csv")):
            root = choose_folder()
    if root is None: raise SystemExit("No source folder selected; configure DATA_FOLDER or --input")
    root = root.expanduser().resolve()
    if not root.exists(): raise FileNotFoundError(root)
    if root.is_file() and root.suffix.lower() not in (".zip", ".csv"):
        raise ValueError("--input must be a directory, ZIP, or CSV")
    return root


def sources(root):
    if root.is_file():
        files = [root]
    else:
        files = sorted(x for x in root.rglob("*") if x.is_file())
    for path in files:
        if path.suffix.lower() == ".csv":
            yield Source("directory", path)
        elif path.suffix.lower() == ".zip":
            try:
                with ZipFile(path) as z:
                    for member in z.namelist():
                        if member.lower().endswith(".csv") and not member.endswith("/"):
                            yield Source("zip", path, member)
            except BadZipFile:
                print("WARNING: invalid ZIP skipped:", path)


def source_key(src):
    name = src.filename.upper()
    if name.endswith("_FUSED_TRENDS.CSV") and "CALIBRATED" not in name:
        role = "hemox"
    elif name.startswith("SP20") and name.endswith(".CSV"):
        role = "sp20"
    else:
        return None
    m = re.search(r"P0*(\d{1,2})[_\-]?(NS0[12]|S0[12])(?=$|[^A-Z0-9])", name)
    if not m or not 1 <= int(m.group(1)) <= 10:
        return None
    return (f"P{int(m.group(1)):02d}", m.group(2), role)


def load_manifest(path):
    if path is None: return {}
    data = pd.read_csv(path)
    if not {"participant", "archive"}.issubset(data.columns):
        raise ValueError("Manifest requires columns participant,archive")
    chosen = {}
    for row in data.itertuples(index=False):
        who = normal_subject(row.participant)
        archive = Path(str(row.archive)).expanduser()
        if not archive.is_absolute(): archive = path.resolve().parent / archive
        archive = archive.resolve()
        if who is None or who in chosen or not archive.is_file() or archive.suffix.lower() != ".zip":
            raise ValueError(f"Invalid/duplicate participant or missing ZIP in manifest: {row}")
        chosen[who] = archive
    return chosen


def choose_source(options, who, manifest):
    # Explicitly selected archive takes precedence over extracted versions.
    if who in manifest:
        options = [s for s in options if s.kind == "zip" and s.path == manifest[who]]
    else:
        extracted = [s for s in options if s.kind == "directory"]
        if extracted: options = extracted
    # Deduplicate identical source records, but never silently select revisions.
    options = list(dict.fromkeys(options))
    if not options: return None
    if len(options) != 1:
        raise ValueError("Ambiguous exports: " + " | ".join(s.description for s in options) +
                         ". Supply --manifest or point --input to one version.")
    return options[0]


def read_bytes(source):
    if source.kind == "directory": return source.path.read_bytes()
    with ZipFile(source.path) as z: return z.read(source.member)


def parse_hemox(source, start, end, offset, model):
    frame = pd.read_csv(BytesIO(read_bytes(source)), low_memory=False)
    need = ("elapsed_s", "fused_spo2_percent", "fused_hr_bpm", "quality_valid")
    missing = [x for x in need if x not in frame.columns]
    if missing: raise ValueError(f"HEMOX FUSED_TRENDS lacks columns: {missing}")
    t = pd.to_numeric(frame.elapsed_s, errors="coerce").to_numpy(float)
    if not len(t) or not np.isfinite(t).all() or np.any(np.diff(t) <= 0):
        raise ValueError("HEMOX elapsed_s must be strictly increasing, unique and finite")
    valid = frame.quality_valid.astype(str).str.strip().str.lower().isin(("true", "1", "yes", "valid"))
    # Explicit writable copies are required with some pandas/NumPy builds
    # (notably Python 3.13 environments), where to_numpy() may return
    # a read-only view. We modify these arrays below when masking invalid data.
    spo = pd.to_numeric(frame.fused_spo2_percent, errors="coerce").to_numpy(dtype=float, copy=True)
    hr = pd.to_numeric(frame.fused_hr_bpm, errors="coerce").to_numpy(dtype=float, copy=True)
    valid_arr = valid.to_numpy(dtype=bool, copy=True)
    if model == "m1":
        spo = np.clip(spo - FROZEN_M1_OFFSET_PP, 0, 100)
    spo[(~valid_arr) | (~np.isfinite(spo)) | (spo < 0) | (spo > 100)] = np.nan
    hr[(~valid_arr) | (~np.isfinite(hr)) | (hr <= 0)] = np.nan
    out = pd.DataFrame({"time_s": t + offset, "spo2": spo, "hr": hr,
                        "quality_valid": valid_arr})
    out = out[out.time_s.between(start, end)].reset_index(drop=True)
    if out.empty: raise ValueError("No HEMOX trend rows in plotting interval")
    return out


def decode_sp20(blob):
    if blob.startswith((b"\xff\xfe", b"\xfe\xff")) or b"\x00" in blob[:150]:
        return blob.decode("utf-16")
    try: return blob.decode("utf-8-sig")
    except UnicodeDecodeError: return blob.decode("cp1252")


def clock_parse(raw):
    raw = " ".join(str(raw).strip().split())
    patterns = ("%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %H:%M:%S",
                "%d/%m/%Y %I:%M:%S %p", "%d/%m/%Y %H:%M:%S",
                "%Y-%m-%d %I:%M:%S %p", "%Y-%m-%d %H:%M:%S")
    for fmt in patterns:
        try: return datetime.strptime(raw, fmt)
        except ValueError: pass
    return None


def parse_sp20(source, start, end, start_elapsed, shift):
    text = decode_sp20(read_bytes(source))
    meta = next((line for line in text.splitlines() if re.search(r"\bStart\s*:", line, re.I)), "")
    hdr = re.search(r"Start:\s*(.+?),\s*End:", meta, re.I)
    start_clock = clock_parse(hdr.group(1)) if hdr else None
    if hdr and start_clock is None: raise ValueError("Unrecognized SP-20 header recording start time")
    stamp = r"(?:\d{1,4}[/.-]){2}\d{1,4}\s+\d{1,2}:\d{2}:\d{2}(?:\s*[AP]M)?"
    row_regex = re.compile(r"^\s*(?P<clock>" + stamp + r")\s*[,;\t]\s*(?P<values>.*)$", re.I)
    records = []
    for line in text.splitlines():
        m = row_regex.match(line)
        if not m: continue
        values = re.split(r"\s*[,;\t]\s*", m.group("values"))
        if len(values) < 2: continue
        dt = clock_parse(m.group("clock"))
        if dt is None: continue
        if not re.search(r"\b[AP]M\b", m.group("clock"), re.I):
            expected = (start_clock + timedelta(seconds=len(records)) if start_clock else
                        records[-1][0] + timedelta(seconds=1) if records else dt)
            dt = min((dt + timedelta(hours=k) for k in (-12, 0, 12)),
                     key=lambda item: abs((item - expected).total_seconds()))
        records.append((dt, values[0], values[1]))
    if not records: raise ValueError("No readable timestamped SP-20 measurement rows")
    frame = pd.DataFrame(records, columns=("stamp", "spo2", "hr"))
    if frame.stamp.duplicated().any() or not frame.stamp.is_monotonic_increasing:
        raise ValueError("SP-20 time is duplicated or nonmonotone")
    frame["time_s"] = (frame.stamp - frame.stamp.iloc[0]).dt.total_seconds() + start_elapsed + shift
    for col in ("spo2", "hr"): frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame.loc[~frame.spo2.between(0, 100), "spo2"] = np.nan
    frame.loc[frame.hr <= 0, "hr"] = np.nan
    frame = frame[frame.time_s.between(start, end)][["time_s", "spo2", "hr"]].reset_index(drop=True)
    if frame.empty: raise ValueError("No SP-20 rows in selected plotting interval")
    return frame


def style(dpi):
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((x for x in ("Arial Narrow", "Arial", "Liberation Sans", "DejaVu Sans")
                 if x in installed), "DejaVu Sans")
    plt.rcParams.update({"font.family": font, "font.weight": "bold",
                         "axes.labelweight": "bold", "axes.labelsize": 12,
                         "xtick.labelsize": 10, "ytick.labelsize": 10,
                         "axes.linewidth": 1.3, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white",
                         "axes.facecolor": "white", "savefig.facecolor": "white",
                         "savefig.dpi": dpi})


def subplot(ax, he, ref, metric, panel, start, end, limits=None, legend=True):
    if ref is not None:
        ax.plot(ref.time_s, ref[metric], color=RED, lw=1.35, label="SP-20", zorder=3)
    if he is not None:
        ax.plot(he.time_s, he[metric], color=BLUE, lw=1.65, label="HEMOX", zorder=2)
    ax.set_xlim(start, end)
    if limits is not None: ax.set_ylim(*limits)
    ax.grid(True, color=GRID, lw=.55, alpha=.55)
    ax.set_axisbelow(True)
    ax.text(.018, .96, panel, transform=ax.transAxes, va="top", fontsize=10.5,
            fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", alpha=.78))
    if legend:
        leg = ax.legend(loc="upper right", ncol=2, fontsize=9, facecolor="white",
                        edgecolor=BORDER, framealpha=.92)
        for entry in leg.get_texts(): entry.set_fontweight("bold")


def save(fig, path, fmt, dpi):
    settings = {"dpi": dpi, "facecolor": "white"}
    if fmt == "jpg": settings["pil_kwargs"] = {"quality": 96, "subsampling": 0}
    fig.savefig(path, **settings)
    plt.close(fig)


def per_run_figures(he, ref, who, run, dest, args):
    if args.figure_set == "none":
        return []
    dest.mkdir(parents=True, exist_ok=True)
    modelsuffix = "M1" if args.spo2_model == "m1" else "M0"
    spec = (("spo2", r"SpO$_2$ (%)", SPO2_Y_LIMITS, "SpO2"),
            ("hr", "HR/PR (bpm)", HR_Y_LIMITS, "HR"))
    images = []
    if args.figure_set in ("both", "singles"):
        for metric, ylabel, limits, slug in spec:
            fig, ax = plt.subplots(figsize=(9.1, 4.25))
            subplot(ax, he, ref, metric, "(a)", args.start, args.end, limits=limits)
            ax.set_ylabel(ylabel); ax.set_xlabel("Elapsed time (s)")
            fig.subplots_adjust(left=.105, right=.985, bottom=.17, top=.975)
            target = dest / f"{who}_{run}_{slug}_vs_SP20_{modelsuffix}.{args.format}"
            save(fig, target, args.format, args.dpi); images.append(target)
    if args.figure_set in ("both", "combined"):
        fig, axes = plt.subplots(2, 1, figsize=(9.1, 7), sharex=True)
        for letter, ax, (metric, ylabel, limits, _) in zip("ab", axes, spec):
            subplot(ax, he, ref, metric, f"({letter})", args.start, args.end, limits=limits)
            ax.set_ylabel(ylabel)
        axes[-1].set_xlabel("Elapsed time (s)")
        fig.subplots_adjust(left=.11, right=.985, top=.97, bottom=.09, hspace=.25)
        target = dest / f"{who}_{run}_SpO2_HR_vs_SP20_{modelsuffix}.{args.format}"
        save(fig, target, args.format, args.dpi); images.append(target)
    return images


def automatic_limits(records, metric):
    arrays = []
    for he, ref in records.values():
        for frame in (he, ref):
            vals = frame[metric].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals): arrays.append(vals)
    if not arrays: return None
    allvalues = np.concatenate(arrays)
    low, high = float(allvalues.min()), float(allvalues.max())
    span = max(high-low, 1.0 if metric == "hr" else .5)
    return (low-.09*span, high+.09*span)


def overview_figure(who, records, statuses, dest, args):
    if not records: return None
    fig, axs = plt.subplots(2, 4, figsize=(18.0, 7.6), sharex=True, sharey="row")
    limits = {"spo2": SPO2_Y_LIMITS or automatic_limits(records, "spo2"),
              "hr": HR_Y_LIMITS or automatic_limits(records, "hr")}
    for col, run in enumerate(RUNS):
        for row, (metric, label) in enumerate((("spo2", r"SpO$_2$ (%)"),
                                                ("hr", "HR/PR (bpm)"))):
            ax = axs[row, col]
            letter = chr(ord('a') + row*4 + col)
            if run in records:
                he, ref = records[run]
                subplot(ax, he, ref, metric, f"({letter}) {run}" if row == 0 else f"({letter})",
                        args.start, args.end, limits=limits[metric], legend=(row == 0 and col == 0))
            else:
                ax.set_xlim(args.start, args.end)
                if limits[metric] is not None: ax.set_ylim(*limits[metric])
                ax.grid(True, color=GRID, lw=.55, alpha=.55)
                ax.text(.018, .96, f"({letter}) {run}" if row == 0 else f"({letter})",
                        transform=ax.transAxes, va="top", fontsize=10.5, fontweight="bold")
                ax.text(.5, .5, statuses.get(run, "Not available"), transform=ax.transAxes,
                        ha="center", va="center", color=GREY, fontsize=11)
            if col == 0: ax.set_ylabel(label)
            if row == 1: ax.set_xlabel("Elapsed time (s)")
    fig.subplots_adjust(left=.065, right=.99, top=.985, bottom=.105, wspace=.18, hspace=.22)
    target = dest / f"{who}_FourRun_SpO2_HR_Overview_{args.spo2_model.upper()}.{args.format}"
    save(fig, target, args.format, args.dpi)
    return target



def global_limits(records_by_key, metric):
    arrays = []
    for he, ref in records_by_key.values():
        for frame in (he, ref):
            vals = frame[metric].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals):
                arrays.append(vals)
    if not arrays:
        return None
    allvalues = np.concatenate(arrays)
    low, high = float(allvalues.min()), float(allvalues.max())
    span = max(high - low, 1.0 if metric == "hr" else 0.5)
    return (low - 0.07 * span, high + 0.07 * span)


def allruns_overview_separate(records_by_key, status_by_key, target_subjects, target_runs, dest, args):
    if not records_by_key:
        return []
    created = []
    specs = (("spo2", r"SpO$_2$ (%)", SPO2_Y_LIMITS or global_limits(records_by_key, "spo2"), "SpO2"),
             ("hr", "HR/PR (bpm)", HR_Y_LIMITS or global_limits(records_by_key, "hr"), "HR"))
    for metric, ylabel, limits, slug in specs:
        nrows, ncols = len(target_subjects), len(target_runs)
        fig_w = max(10.5, 3.35 * ncols)
        fig_h = max(8.0, 1.70 * nrows)
        fig, axs = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), sharex=True, sharey=True, squeeze=False)
        for col, run in enumerate(target_runs):
            axs[0, col].set_title(run, fontweight="bold", fontsize=11, pad=7)
        first_with_legend = True
        for row, who in enumerate(target_subjects):
            for col, run in enumerate(target_runs):
                ax = axs[row, col]
                key = (who, run)
                if key in records_by_key:
                    he, ref = records_by_key[key]
                    subplot(ax, he, ref, metric, "", args.start, args.end, limits=limits, legend=first_with_legend)
                    first_with_legend = False
                else:
                    ax.set_xlim(args.start, args.end)
                    if limits is not None:
                        ax.set_ylim(*limits)
                    ax.grid(True, color=GRID, lw=.55, alpha=.55)
                    ax.set_axisbelow(True)
                    ax.text(.5, .5, status_by_key.get(key, "Not available"), transform=ax.transAxes,
                            ha="center", va="center", color=GREY, fontsize=9.5, fontweight="bold")
                if col == 0:
                    ax.set_ylabel(f"{who}\n{ylabel}")
                if row == nrows - 1:
                    ax.set_xlabel("Elapsed time (s)")
        fig.subplots_adjust(left=.085, right=.99, top=.965, bottom=.065, wspace=.10, hspace=.20)
        target = dest / f"AllParticipants_AllRuns_{slug}_vs_SP20_{args.spo2_model.upper()}.{args.format}"
        save(fig, target, args.format, args.dpi)
        created.append(target)
    return created


def allruns_overview_combined(records_by_key, status_by_key, target_subjects, target_runs, dest, args):
    if not records_by_key:
        return None
    nrows, ncols = 2 * len(target_subjects), len(target_runs)
    fig_w = max(10.5, 3.35 * ncols)
    fig_h = max(12.0, 1.50 * nrows)
    fig, axs = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), sharex=True, squeeze=False)
    limits = {"spo2": SPO2_Y_LIMITS or global_limits(records_by_key, "spo2"),
              "hr": HR_Y_LIMITS or global_limits(records_by_key, "hr")}
    for col, run in enumerate(target_runs):
        axs[0, col].set_title(run, fontweight="bold", fontsize=11, pad=7)
    first_with_legend = True
    for subj_idx, who in enumerate(target_subjects):
        for col, run in enumerate(target_runs):
            key = (who, run)
            ax_spo2 = axs[2*subj_idx, col]
            ax_hr = axs[2*subj_idx + 1, col]
            if key in records_by_key:
                he, ref = records_by_key[key]
                subplot(ax_spo2, he, ref, "spo2", "", args.start, args.end, limits=limits["spo2"], legend=first_with_legend)
                first_with_legend = False
                subplot(ax_hr, he, ref, "hr", "", args.start, args.end, limits=limits["hr"], legend=False)
            else:
                for ax, metric in ((ax_spo2, "spo2"), (ax_hr, "hr")):
                    ax.set_xlim(args.start, args.end)
                    if limits[metric] is not None:
                        ax.set_ylim(*limits[metric])
                    ax.grid(True, color=GRID, lw=.55, alpha=.55)
                    ax.set_axisbelow(True)
                    ax.text(.5, .5, status_by_key.get(key, "Not available"), transform=ax.transAxes,
                            ha="center", va="center", color=GREY, fontsize=9.0, fontweight="bold")
            if col == 0:
                ax_spo2.set_ylabel(f"{who}\nSpO$_2$ (%)")
                ax_hr.set_ylabel(f"{who}\nHR/PR (bpm)")
            if subj_idx == len(target_subjects) - 1:
                ax_hr.set_xlabel("Elapsed time (s)")
    fig.subplots_adjust(left=.095, right=.99, top=.97, bottom=.05, wspace=.10, hspace=.18)
    target = dest / f"AllParticipants_AllRuns_SpO2_HR_vs_SP20_{args.spo2_model.upper()}.{args.format}"
    save(fig, target, args.format, args.dpi)
    return target



def smooth_mean_trace(frame, metric, window_s):
    if frame is None or frame.empty or window_s <= 1:
        return frame
    window_s = int(window_s)
    if window_s % 2 == 0:
        window_s += 1
    out = frame.copy()
    out[metric] = (out[metric]
                   .rolling(window=window_s, center=True, min_periods=1)
                   .mean())
    return out


def allruns_overlay_figure(records_by_key, target_subjects, target_runs, dest, args):
    if not records_by_key:
        return None
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 7.3), sharex=True)
    # Limits use the full set of curves that may appear in each panel.
    spo2_arrays, hr_arrays = [], []
    for he_m0, he_sel, ref in records_by_key.values():
        for frame in (he_m0, he_sel, ref):
            if frame is not None and "spo2" in frame:
                vals = frame["spo2"].to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if len(vals): spo2_arrays.append(vals)
            if frame is not None and "hr" in frame:
                vals = frame["hr"].to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if len(vals): hr_arrays.append(vals)
    def limits_from(arrays, fallback_span):
        if not arrays: return None
        allv = np.concatenate(arrays)
        low, high = float(allv.min()), float(allv.max())
        span = max(high - low, fallback_span)
        return (low - 0.07 * span, high + 0.07 * span)
    limits = {"spo2": SPO2_Y_LIMITS or limits_from(spo2_arrays, 0.5),
              "hr": HR_Y_LIMITS or limits_from(hr_arrays, 1.0)}

    # Faint individual traces: SP-20 and HEMOX M0 only (selected SpO2 differs only by frozen offset,
    # so plotting its individual traces as well adds clutter without new shape information).
    first = True
    run_count = 0
    for who in target_subjects:
        for run in target_runs:
            key = (who, run)
            if key not in records_by_key:
                continue
            he_m0, he_sel, ref = records_by_key[key]
            run_count += 1
            label_ref = "SP-20" if first else None
            label_he_m0 = "HEMOX M0" if first else None
            axes[0].plot(ref.time_s, ref["spo2"], color=RED, lw=.75, alpha=.16, label=label_ref, zorder=2)
            axes[0].plot(he_m0.time_s, he_m0["spo2"], color=BLUE, lw=.75, alpha=.16, label=label_he_m0, zorder=1)
            axes[1].plot(ref.time_s, ref["hr"], color=RED, lw=.75, alpha=.16, label=label_ref, zorder=2)
            axes[1].plot(he_m0.time_s, he_m0["hr"], color=BLUE, lw=.75, alpha=.16, label=label_he_m0, zorder=1)
            first = False

    def mean_trace(metric, source):
        per_run = []
        for he_m0, he_sel, ref in records_by_key.values():
            frame = {"sp20": ref, "hemox_m0": he_m0, "hemox_sel": he_sel}[source][["time_s", metric]].dropna().copy()
            if frame.empty:
                continue
            frame["time_bin_s"] = np.rint(frame["time_s"].to_numpy(float)).astype(int)
            frame = frame[frame["time_bin_s"].between(int(np.floor(args.start)), int(np.ceil(args.end)))]
            if frame.empty:
                continue
            one_run = frame.groupby("time_bin_s", sort=True)[metric].mean().reset_index()
            per_run.append(one_run)
        if not per_run:
            return None
        merged = pd.concat(per_run, ignore_index=True)
        agg = merged.groupby("time_bin_s", sort=True)[metric].mean().reset_index()
        agg.rename(columns={"time_bin_s": "time_s"}, inplace=True)
        return agg

    sp20_spo2 = smooth_mean_trace(mean_trace("spo2", "sp20"), "spo2", args.smooth_window)
    hemox_m0_spo2 = smooth_mean_trace(mean_trace("spo2", "hemox_m0"), "spo2", args.smooth_window)
    hemox_sel_spo2 = smooth_mean_trace(mean_trace("spo2", "hemox_sel"), "spo2", args.smooth_window)
    sp20_hr = smooth_mean_trace(mean_trace("hr", "sp20"), "hr", args.smooth_window)
    hemox_m0_hr = smooth_mean_trace(mean_trace("hr", "hemox_m0"), "hr", args.smooth_window)

    if sp20_spo2 is not None:
        axes[0].plot(sp20_spo2.time_s, sp20_spo2.spo2, color=RED, lw=3.4, alpha=1.0, zorder=6, label="SP-20")
    if hemox_m0_spo2 is not None:
        axes[0].plot(hemox_m0_spo2.time_s, hemox_m0_spo2.spo2, color=BLUE, lw=3.1, alpha=1.0, zorder=5, label="HEMOX M0")
    if hemox_sel_spo2 is not None:
        axes[0].plot(hemox_sel_spo2.time_s, hemox_sel_spo2.spo2, color="#009E73", lw=3.1, alpha=1.0, zorder=7, label="HEMOX M1")
    if sp20_hr is not None:
        axes[1].plot(sp20_hr.time_s, sp20_hr.hr, color=RED, lw=3.4, alpha=1.0, zorder=6, label="SP-20")
    if hemox_m0_hr is not None:
        axes[1].plot(hemox_m0_hr.time_s, hemox_m0_hr.hr, color=BLUE, lw=3.1, alpha=1.0, zorder=5, label="HEMOX M0")

    for ax, metric, ylabel, panel in ((axes[0], "spo2", r"SpO$_2$ (%)", "(a)"),
                                      (axes[1], "hr", "HR/PR (bpm)", "(b)")):
        ax.set_xlim(args.start, args.end)
        if limits[metric] is not None:
            ax.set_ylim(*limits[metric])
        ax.grid(True, color=GRID, lw=.55, alpha=.55)
        ax.set_axisbelow(True)
        ax.set_ylabel(ylabel)
        ax.text(.018, .96, panel, transform=ax.transAxes, va="top", fontsize=10.5,
                fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", alpha=.78))
    axes[-1].set_xlabel("Elapsed time (s)")
    leg0 = axes[0].legend(loc="upper right", ncol=3, fontsize=8.8, facecolor="white",
                          edgecolor=BORDER, framealpha=.92)
    for entry in leg0.get_texts():
        entry.set_fontweight("bold")
    leg1 = axes[1].legend(loc="upper right", ncol=2, fontsize=8.8, facecolor="white",
                          edgecolor=BORDER, framealpha=.92)
    for entry in leg1.get_texts():
        entry.set_fontweight("bold")
    fig.subplots_adjust(left=.095, right=.985, top=.975, bottom=.09, hspace=.23)
    cohort_slug = "P03_P10" if tuple(target_subjects) == HELDOUT_SUBJECTS else ("P01_P02" if tuple(target_subjects) == DEVELOPMENT_SUBJECTS else "P01_P10")
    target = dest / f"{cohort_slug}_M0_vs_SelectedConfig_SpO2M1_HRM0.{args.format}"
    save(fig, target, args.format, args.dpi)
    return target


def validate(args):
    values = (args.start, args.end, args.hemox_offset, args.sp20_start_elapsed, args.sp20_shift)
    if not np.isfinite(values).all() or args.start < 0 or args.end <= args.start:
        raise ValueError("Require finite 0 <= start < end and finite time adjustments")
    if args.dpi < 72: raise ValueError("DPI must be >= 72")
    if int(args.smooth_window) < 1: raise ValueError("--smooth-window must be >= 1")
    for label, bounds in (("SPO2_Y_LIMITS", SPO2_Y_LIMITS), ("HR_Y_LIMITS", HR_Y_LIMITS)):
        if bounds is not None and (len(bounds) != 2 or not np.isfinite(bounds).all() or bounds[0] >= bounds[1]):
            raise ValueError(f"{label} must be None or finite (min,max)")


def main(argv=None):
    args = parser(argv)
    validate(args)
    target_subjects, target_runs = selected_subjects(args.participant, args.participant_group), selected_runs(args.run)
    root = root_folder(args)
    manifest = load_manifest(args.manifest)
    external_sp20 = args.sp20 or (Path(SP20_FILE) if SP20_FILE.strip() else None)
    if external_sp20 is not None:
        if len(target_subjects) != 1 or len(target_runs) != 1:
            raise ValueError("--sp20 / SP20_FILE requires a single participant and run")
        external_sp20 = external_sp20.expanduser().resolve()
        if not external_sp20.is_file() or external_sp20.suffix.lower() != ".csv":
            raise FileNotFoundError(f"External SP-20 CSV does not exist: {external_sp20}")
    out = (args.out_dir or (Path(OUTPUT_FOLDER) if OUTPUT_FOLDER.strip() else
           Path(__file__).resolve().parent / "HEMOX_vs_SP20_AllParticipants_outputs")).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    style(args.dpi)
    catalog = {}
    for src in sources(root):
        key = source_key(src)
        if key: catalog.setdefault(key, []).append(src)
    if external_sp20 is not None:
        catalog[(target_subjects[0], target_runs[0], "sp20")] = [Source("directory", external_sp20)]
    reports, saved = [], []
    all_records, all_status = {}, {}
    for who in target_subjects:
        participant_data, status_by_run = {}, {}
        for run in target_runs:
            row = {"participant": who, "run": run, "model_spo2": "M0_vs_selected",
                   "model_hr": "M0", "alignment": "PROVISIONAL first-row elapsed origin; check offsets",
                   "hemox_offset_s": args.hemox_offset,
                   "sp20_first_row_s": args.sp20_start_elapsed + args.sp20_shift}
            try:
                hsrc = choose_source(catalog.get((who, run, "hemox"), []), who, manifest)
                rsrc = choose_source(catalog.get((who, run, "sp20"), []), who, manifest)
                row["hemox_source"] = hsrc.description if hsrc else ""
                row["sp20_source"] = rsrc.description if rsrc else ""
                if not hsrc or not rsrc:
                    row["status"] = ("not_reported" if (who, run) == ("P08", "NS02") else
                                     "no_matching_source")
                    status_by_run[run] = "Not reported" if row["status"] == "not_reported" else "Not available"
                    all_status[(who, run)] = status_by_run[run]
                    print(who, run, row["status"])
                else:
                    # Build BOTH configurations from the same native export:
                    #   M0: native HEMOX SpO2 / native HR
                    #   selected: SpO2 M1 (frozen P01-P02 offset), HR M0
                    he_m0 = parse_hemox(hsrc, args.start, args.end, args.hemox_offset, "m0")
                    he_sel = parse_hemox(hsrc, args.start, args.end, args.hemox_offset, "m1")
                    ref = parse_sp20(rsrc, args.start, args.end, args.sp20_start_elapsed, args.sp20_shift)
                    dest = out / who / run
                    images = []
                    if args.figure_set != "none":
                        dest.mkdir(parents=True, exist_ok=True)
                        he_m0.to_csv(dest / "HEMOX_M0_plotted_trends.csv", index=False)
                        he_sel.to_csv(dest / "HEMOX_Selected_SpO2M1_HRM0_plotted_trends.csv", index=False)
                        ref.to_csv(dest / "SP20_plotted_trends.csv", index=False)
                    row.update(status="ok", figure_count=len(images),
                               hemox_rows=len(he_sel), sp20_rows=len(ref),
                               hemox_spo2_valid=int(he_sel.spo2.notna().sum()),
                               hemox_hr_valid=int(he_sel.hr.notna().sum()),
                               output_dir=(str(dest) if args.figure_set != "none" else str(out)))
                    participant_data[run] = (he_sel, ref)
                    all_records[(who, run)] = (he_m0, he_sel, ref)
                    all_status[(who, run)] = "ok"
                    saved.extend(images)
                    print(who, run, "OK: loaded M0 + selected configuration + SP-20")
            except (ValueError, OSError, KeyError, BadZipFile) as exc:
                row["status"] = "ERROR: " + str(exc)
                status_by_run[run] = "Source error"
                all_status[(who, run)] = "Source error"
                print(who, run, row["status"])
            reports.append(row)
        # Overlay-only workflow: intentionally do NOT generate per-participant
        # four-run overview/atlas figures here. All valid runs are accumulated
        # into all_records and plotted together once after the batch loop.
    if args.allruns_overview != "none" and all_records:
        extra = allruns_overlay_figure(all_records, target_subjects, target_runs, out, args)
        if extra:
            saved.append(extra)
            print("All-runs overlay:", extra)
    manifest_path = out / "run_manifest.csv"
    pd.DataFrame(reports).to_csv(manifest_path, index=False)
    (out / "README_QA.txt").write_text(
        "Figures represent exported HEMOX and SP-20 time trends, not timestamp-paired "
        "accuracy statistics. Zero-offset first-row comparison is PROVISIONAL unless "
        "independently verified. HEMOX quality-invalid outputs stay NaN, with no "
        "interpolation/carry-forward. SpO2 M0 native or frozen M1; HR M0. "
        "Unavailable runs are not imputed. This script creates one global overlay figure only. A participant's repeated runs/windows "
        "are not independent participants.\n", encoding="utf-8")
    print("Saved", len(saved), "figure(s) to", out)
    print("Run QA:", manifest_path)
    return out


if __name__ == "__main__":
    main()
