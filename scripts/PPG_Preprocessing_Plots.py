
from __future__ import annotations

import argparse
import hashlib
from io import BytesIO
from pathlib import Path
import re
from zipfile import ZipFile, BadZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

# ------------------------- USER SETTINGS ------------------------------------
# A valid Windows example of the path in the supplied script. Change as needed.
# For all P01--P10, point to the containing analysis folder rather than only P01.
DATA_FOLDER = r"E:\post_doc\research proposal\32-Iot based healthcare system\1-proposal start 4-2026\12-measurment\analysis\P09"
OUTPUT_FOLDER = r""  # blank -> HEMOX_PPG_R2_Styled_outputs beside this .py
START_S = 100.0
END_S = 150.0
DPI = 600
FIGURE_VARIANTS = "clean"  # clean (default), full, or both
IMAGE_FORMAT = "jpg"  # png or jpg (both 600 dpi by default)
FS = 50
PARTICIPANTS = tuple(f"P{i:02d}" for i in range(1, 11))
RUNS = ("S01", "S02", "NS01", "NS02")
MOTION_RUNS = ("NS01", "NS02")
STATIONARY_RUNS = ("S01", "S02")
GENERATE_MOTION_FIGURE = True
GENERATE_STATIONARY_FIGURE = True
MOTION_LAYOUT = "2panel"  # 2panel (CH1, CH2 only) or 3panel (adds 2-s IR AC RMS)
MOTION_MIN_DURATION_S = 295.0  # avoid annotating protocol on truncated runs
STATIONARY_MIN_DURATION_S = 295.0  # avoid annotating full-run stationary figure on truncated runs
MOTION_WINDOWS = (
    (120.0, 130.0, "RH", "Scheduled RH motion", "#D55E00"),
    (180.0, 190.0, "LH", "Scheduled LH motion", "#009E73"),
    (240.0, 250.0, "Both", "Scheduled both-hands motion", "#C62828"),
)

# Use None for autoscale, or (minimum, maximum) for a fixed axis. Panel names:
# CH1-RAW, CH1-IR, CH1-RED, CH2-RAW, CH2-IR, CH2-RED.
# Raw panels: left=IR, right=RED. IR/RED panels: left=AC, right=absolute DC.
LEFT_Y_LIMITS = {f"CH{ch}-{part}": None for ch in (1, 2) for part in ("RAW", "IR", "RED")}
RIGHT_Y_LIMITS = {
    "CH1-RAW": None, "CH1-IR": (0, 120000), "CH1-RED": (0, 120000),
    "CH2-RAW": None, "CH2-IR": (0, 120000), "CH2-RED": (0, 120000),
}
# ---------------------------------------------------------------------------

# Colors and typography locked to the accompanying R2 correlation script.
BLUE = "#0072B2"    # IR and normal waveform
RED = "#C62828"     # RED waveform (also regression in correlation figures)
GREY = "#777F87"    # DC/reference overlay
ORANGE = "#D55E00"
GREEN = "#009E73"
GRID = "#CFD5DC"
INK = "#253348"
BORDER = "#B0BAC3"


def setup_plot_style() -> str:
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((name for name in ("Arial Narrow", "Arial", "Liberation Sans", "DejaVu Sans")
                 if name in installed), "DejaVu Sans")
    plt.rcParams.update({
        "font.family": font, "font.weight": "bold", "axes.labelweight": "bold",
        "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11,
        "axes.linewidth": 1.3, "axes.spines.top": False,
        "axes.spines.right": False, "savefig.dpi": DPI,
        "axes.titleweight": "bold", "legend.fontsize": 9,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "mathtext.default": "regular",
    })
    return font


def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, help="RAW CSV, ZIP or root folder; overrides DATA_FOLDER")
    p.add_argument("--data-dir", type=Path, help="Alias for --input")
    p.add_argument("--participant", help="Optional participant filter, e.g. P09")
    p.add_argument("--run", help="Optional session filter: S01, S02, NS01, NS02")
    p.add_argument("--start", type=float, default=START_S)
    p.add_argument("--end", type=float, default=END_S)
    p.add_argument("--left-ylim", action="append", nargs=3, default=[], metavar=("PANEL", "MIN", "MAX"),
                   help="Override left y-axis; repeat for additional panels")
    p.add_argument("--right-ylim", action="append", nargs=3, default=[], metavar=("PANEL", "MIN", "MAX"),
                   help="Override right y-axis; e.g. CH1-IR 58500 58900")
    p.add_argument("--out-dir", type=Path, help="Output root directory")
    p.add_argument("--dpi", type=int, default=DPI)
    p.add_argument("--variants", choices=("both", "full", "clean"), default=FIGURE_VARIANTS,
                   help="Save panel-labeled, minimal-label, or both versions; neither has a title/footer")
    p.add_argument("--format", choices=("png", "jpg"), default=IMAGE_FORMAT, dest="image_format")
    p.add_argument("--no-motion", action="store_true", help="Skip controlled-motion figure for NS01/NS02")
    p.add_argument("--motion-layout", choices=("2panel", "3panel"), default=MOTION_LAYOUT,
                   help="Two filtered IR traces, optionally plus 2-s RMS (default: 2panel)")
    p.add_argument("--no-stationary", action="store_true",
                   help="Skip full-run stationary IR-AC figure for S01/S02")
    return p.parse_args(argv)


def axis_limit_overrides(specifications):
    valid = set(LEFT_Y_LIMITS)
    changed = {}
    for panel, minimum, maximum in specifications:
        panel = panel.strip().upper()
        if panel not in valid:
            raise ValueError(f"Unknown axis panel {panel!r}; valid: {', '.join(sorted(valid))}")
        try:
            lo, hi = float(minimum), float(maximum)
        except ValueError as exc:
            raise ValueError(f"{panel}: MIN and MAX must be numeric") from exc
        if not np.isfinite([lo, hi]).all() or lo >= hi:
            raise ValueError(f"{panel}: axis limits require finite MIN < MAX")
        changed[panel] = (lo, hi)
    return changed


def validate_axis_limits(settings, label):
    if set(settings) != set(LEFT_Y_LIMITS):
        raise ValueError(f"{label} requires exactly {sorted(LEFT_Y_LIMITS)}")
    for key, bounds in settings.items():
        if bounds is not None and (len(bounds) != 2 or not np.isfinite(bounds).all() or bounds[0] >= bounds[1]):
            raise ValueError(f"{label}[{key!r}] must be None or finite (MIN, MAX)")


def identity(text):
    m = re.search(r"P0*(\d{1,2})[_\-]?(NS0[12]|S0[12])(?=$|[^A-Z0-9])", str(text).upper())
    if m and 1 <= int(m.group(1)) <= 10:
        return (f"P{int(m.group(1)):02d}", m.group(2))
    return None


def choose_data_folder(initial):
    try:
        from tkinter import Tk, filedialog
        root = Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        name = filedialog.askdirectory(title="Choose folder with HEMOX P01-P10 RAW exports", initialdir=str(initial))
        root.destroy()
        return Path(name) if name else None
    except (ImportError, OSError, RuntimeError):
        return None


def raw_candidates(path):
    """Enumerate real CSVs and ZIP members without interpreting ZIP paths as local files."""
    path = path.expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Input does not exist: {path}. Check DATA_FOLDER or pass --input.")
    if path.is_file() and path.suffix.lower() == ".csv":
        return [(path, None)] if path.name.upper().endswith("_RAW.CSV") else []
    if path.is_file() and path.suffix.lower() == ".zip":
        with ZipFile(path) as z:
            return [(path, item) for item in z.namelist()
                    if not item.endswith("/") and Path(item).name.upper().endswith("_RAW.CSV")]
    if not path.is_dir():
        raise ValueError(f"Input must be a directory, ZIP or *_RAW.csv: {path}")
    paths = sorted(p for p in path.rglob("*") if p.is_file())
    matches = [(p, None) for p in paths if p.name.upper().endswith("_RAW.CSV")]
    for archive in (p for p in paths if p.suffix.lower() == ".zip"):
        try:
            with ZipFile(archive) as z:
                matches.extend((archive, item) for item in z.namelist()
                               if not item.endswith("/") and Path(item).name.upper().endswith("_RAW.CSV"))
        except BadZipFile:
            print(f"Warning: skipping invalid ZIP: {archive}")
    return matches


def source_identity(path, member):
    """Find subject and run from member/file name or enclosing ZIP/folder context."""
    if member:
        for piece in (Path(member).name, member):
            key = identity(piece)
            if key:
                return key
    key = identity(path.name)
    if key:
        return key
    pieces = ((Path(member).stem, *Path(member).parts) if member else ()) + (path.stem,) + tuple(p.name for p in path.parents)
    subject = None
    for piece in pieces:
        m = re.match(r"P0*(\d{1,2})(?!\d)", piece, re.I)
        if m and 1 <= int(m.group(1)) <= 10:
            subject = f"P{int(m.group(1)):02d}"
            break
    session = None
    for piece in pieces:
        m = re.search(r"(?<![A-Z0-9])(NS0[12]|S0[12])(?=$|[^A-Z0-9])", piece, re.I)
        if m:
            session = m.group(1).upper()
            break
        # Some source names use e.g. P09S02 without a separator.
        m = re.search(r"P0*\d{1,2}[_\-]?(NS0[12]|S0[12])(?=$|[^A-Z0-9])", piece, re.I)
        if m:
            session = m.group(1).upper()
            break
    return (subject, session) if subject and session else None


def select_sources(args):
    data_path = args.input or args.data_dir or (Path(DATA_FOLDER) if DATA_FOLDER.strip() else None)
    if data_path is None:
        data_path = choose_data_folder(Path(__file__).resolve().parent)
    if data_path is None:
        raise SystemExit("No input selected. Set DATA_FOLDER or pass --input.")
    candidates = raw_candidates(data_path)
    if args.participant:
        m = re.fullmatch(r"P0*(\d{1,2})", args.participant.strip(), re.I)
        if not m or not 1 <= int(m.group(1)) <= 10:
            raise ValueError("--participant must be P01 through P10")
        args.participant = f"P{int(m.group(1)):02d}"
    if args.run:
        args.run = args.run.strip().upper()
        if args.run not in RUNS:
            raise ValueError(f"--run must be one of {RUNS}")
    jobs = []
    for path, member in candidates:
        key = source_identity(path, member)
        if (key and key[0] in PARTICIPANTS and key[1] in RUNS
                and (not args.participant or args.participant == key[0])
                and (not args.run or args.run == key[1])):
            jobs.append((key[0], key[1], path, member))
    jobs.sort(key=lambda j: (j[0], RUNS.index(j[1]), str(j[2]), j[3] or ""))
    if not jobs:
        raise FileNotFoundError(f"No matching P01-P10 *_RAW.csv found inside {data_path}. "
                                "Check --input and participant/session names.")
    return jobs


def read_raw(path, member):
    if member is None:
        frame = pd.read_csv(path, low_memory=False)
    else:
        with ZipFile(path) as z:
            frame = pd.read_csv(BytesIO(z.read(member)), low_memory=False)
    fields = ("common_time_us", "ch1_red_raw", "ch1_ir_raw", "ch2_red_raw", "ch2_ir_raw")
    absent = [field for field in fields if field not in frame]
    if absent:
        raise ValueError(f"RAW CSV lacks columns: {absent}")
    time = pd.to_numeric(frame["common_time_us"], errors="coerce").to_numpy(float) / 1e6
    if len(time) < 3 or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError("common_time_us needs at least 3 finite strictly increasing timestamps")
    if not .018 <= np.median(np.diff(time)) <= .022:
        raise ValueError("Expected HEMOX 50 Hz RAW export (median sample interval ~0.02 s)")
    raw = {}
    for ch in (1, 2):
        for wave in ("red", "ir"):
            name = f"ch{ch}_{wave}"
            arr = pd.to_numeric(frame[name + "_raw"], errors="coerce").to_numpy(float)
            if not np.isfinite(arr).all():
                raise ValueError(f"{name}_raw has missing/non-numeric samples")
            raw[name] = arr
    return time, raw


def rolling_median(values, radius):
    n = len(values)
    result = np.empty(n, dtype=float)
    if n >= 2 * radius + 1:
        result[radius:n-radius] = np.median(
            np.lib.stride_tricks.sliding_window_view(values, 2*radius+1), axis=1)
        edges = list(range(radius)) + list(range(n-radius, n))
    else:
        edges = range(n)
    for i in edges:
        result[i] = np.median(values[max(0, i-radius):min(n, i+radius+1)])
    return result


def hampel(values):
    radius = 6
    median = rolling_median(values, radius)
    mad = rolling_median(np.abs(values - median), radius)
    scale = 1.4826 * mad
    replaced = (scale > 0) & (np.abs(values - median) > 6 * scale)
    replaced[0] = replaced[-1] = False
    result = values.copy()
    result[replaced] = median[replaced]
    return result


def biquad_coefficients(low_pass, cutoff):
    omega = 2 * np.pi * cutoff / FS
    cosine, sine = np.cos(omega), np.sin(omega)
    alpha = sine / (2 * np.sqrt(.5))
    if low_pass:
        b0, b1, b2 = (1-cosine)/2, 1-cosine, (1-cosine)/2
    else:
        b0, b1, b2 = (1+cosine)/2, -(1+cosine), (1+cosine)/2
    a0 = 1 + alpha
    return b0/a0, b1/a0, b2/a0, -2*cosine/a0, (1-alpha)/a0


def biquad(values, coefficients):
    out = np.empty(len(values), dtype=float)
    x1 = x2 = float(values[0])
    y1 = y2 = 0.0
    for i, x0 in enumerate(values):
        y = (coefficients[0]*x0 + coefficients[1]*x1 + coefficients[2]*x2
             - coefficients[3]*y1 - coefficients[4]*y2)
        out[i] = y
        x2, x1, y2, y1 = x1, x0, y1, y
    return out


def zero_phase(values, coefficients):
    return biquad(biquad(values, coefficients)[::-1], coefficients)[::-1]


def process(raw):
    smooth, ac = {}, {}
    for key, signal in raw.items():
        smooth[key] = hampel(signal)
        highpass = zero_phase(smooth[key], biquad_coefficients(False, .5))
        ac[key] = zero_phase(highpass, biquad_coefficients(True, 4.0))
    return smooth, ac


def dc_windows(raw):
    n = len(next(iter(raw.values())))
    width, step = 8*FS, 2*FS
    starts = range(0, n-width+1, step)
    data = {"time_s": [(s + width/2)/FS for s in starts]}
    for key, values in raw.items():
        data[key] = [abs(float(np.mean(values[s:s+width]))) for s in starts]
    return pd.DataFrame(data)


def format_axes(ax, right=None):
    ax.grid(True, lw=.55, color=GRID, alpha=.5)
    ax.set_axisbelow(True)
    ax.spines["left"].set_linewidth(1.3)
    ax.spines["bottom"].set_linewidth(1.3)
    ax.spines["top"].set_visible(False)
    if right is not None:
        right.spines["top"].set_visible(False)
        right.spines["right"].set_visible(True)
        right.spines["right"].set_color(BORDER)
        right.spines["right"].set_linewidth(1.3)
        right.grid(False)


def legend_for(ax, handles):
    legend = ax.legend(handles=handles, loc="upper right", ncol=2, fontsize=9,
                       facecolor="white", framealpha=.94, edgecolor=BORDER,
                       handlelength=2.1, columnspacing=1.1)
    for label in legend.get_texts():
        label.set_fontweight("bold")


def plot_channel(ch, time, raw, ac, dc, start, end, out, dpi, left_limits,
                 right_limits, image_format, participant, run, clean=False):
    show = (time >= start) & (time <= end)
    if not np.any(show):
        raise ValueError(f"No samples in plotting interval {start:g}-{end:g} s")
    fig, axes = plt.subplots(3, 1, figsize=(10.2, 8.6), sharex=True)
    ir, red = f"ch{ch}_ir", f"ch{ch}_red"
    # (a) Native raw waveforms; colors here identify wavelengths, not channels.
    raw_red_axis = axes[0].twinx()
    line_ir, = axes[0].plot(time[show], raw[ir][show], color=BLUE, lw=1.4, label="Raw IR")
    line_red, = raw_red_axis.plot(time[show], raw[red][show], color=RED, lw=1.4, label="Raw RED")
    axes[0].set_ylabel("IR raw (counts)", color=BLUE)
    raw_red_axis.set_ylabel("RED raw (counts)", color=RED)
    axes[0].tick_params(axis="y", colors=BLUE)
    raw_red_axis.tick_params(axis="y", colors=RED)
    if left_limits[f"CH{ch}-RAW"] is not None:
        axes[0].set_ylim(*left_limits[f"CH{ch}-RAW"])
    if right_limits[f"CH{ch}-RAW"] is not None:
        raw_red_axis.set_ylim(*right_limits[f"CH{ch}-RAW"])
    axes[0].text(.012, .96, "(a)" if clean else "(a) Raw IR and RED",
                 transform=axes[0].transAxes, ha="left", va="top", weight="bold", fontsize=11,
                 bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
    legend_for(axes[0], [line_ir, line_red])
    format_axes(axes[0], raw_red_axis)

    subset = dc.loc[dc.time_s.between(start, end)]
    for ax, key, name, color, letter in (
            (axes[1], ir, "IR", BLUE, "b"), (axes[2], red, "RED", RED, "c")):
        dc_axis = ax.twinx()
        ac_line, = ax.plot(time[show], ac[key][show], color=color, lw=1.45,
                           label=f"{name} AC")
        dc_line, = dc_axis.plot(subset.time_s, subset[key], ls="--", lw=1.5,
                                marker="o", markersize=4.0, markerfacecolor="white",
                                markeredgewidth=1.2, color=GREY,
                                label=f"{name} DC")
        ax.set_ylabel(f"{name} AC (counts)", color=color)
        dc_axis.set_ylabel(f"{name} DC (counts)", color=GREY)
        ax.tick_params(axis="y", colors=color)
        dc_axis.tick_params(axis="y", colors=GREY)
        if left_limits[f"CH{ch}-{name}"] is not None:
            ax.set_ylim(*left_limits[f"CH{ch}-{name}"])
        if right_limits[f"CH{ch}-{name}"] is not None:
            dc_axis.set_ylim(*right_limits[f"CH{ch}-{name}"])
        ax.text(.012, .96, f"({letter})" if clean else f"({letter}) {name} AC and DC",
                transform=ax.transAxes, ha="left", va="top", weight="bold", fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
        legend_for(ax, [ac_line, dc_line])
        format_axes(ax, dc_axis)
        if subset.empty:
            ax.text(.5, .08, "No DC window center within selected interval",
                    transform=ax.transAxes, ha="center", fontsize=9,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
    for ax in axes:
        ax.set_xlim(start, end)
        ax.tick_params(axis="x", labelsize=11)
    axes[-1].set_xlabel("Elapsed time (s)")
    # All variants: no global title or explanatory footer. Retain panel IDs, axes and legends.
    fig.subplots_adjust(left=.145, right=.855, top=.975, bottom=.075, hspace=.30)
    suffix = "_Clean" if clean else ""
    output = out / f"CH{ch}_Raw_IR_RED_AC_DC_{start:g}-{end:g}s{suffix}.{image_format}"
    save_options = dict(dpi=dpi, facecolor="white")
    if image_format == "jpg":
        save_options["pil_kwargs"] = {"quality": 96, "subsampling": 0}
    fig.savefig(output, **save_options)
    plt.close(fig)
    return output



def ir_ac_rms(time_elapsed, values, bin_seconds=2.0):
    """Observed 2-s RMS of filtered IR AC; NOT motion detection or signal-quality score.

    Bin membership uses actual timestamps, so a small timing gap does not shift
    subsequent window boundaries. Retains every finite sample inside each bin.
    """
    if len(time_elapsed) != len(values):
        raise ValueError("RMS time/value lengths differ")
    bins = np.floor(np.maximum(0, time_elapsed) / bin_seconds).astype(int)
    result = []
    for index in np.unique(bins):
        active = bins == index
        wave = np.asarray(values[active], dtype=float)
        if wave.size:
            result.append((index*bin_seconds + bin_seconds/2,
                           float(np.sqrt(np.mean(np.square(wave))))))
    return pd.DataFrame(result, columns=("elapsed_s", "ir_ac_rms_counts"))


def plot_controlled_motion(time, ac, out, dpi, image_format, layout="3panel"):
    """Full-run protocol visualization. Bands are SCHEDULED, not measured motion.

    Original 50-Hz filtered IR AC samples are plotted (no plot-time smoothing).
    Third panel summarizes amplitude as 2-s RMS; it is not an SQI or eligibility flag.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    elapsed = np.asarray(time, dtype=float) - float(time[0])
    if elapsed[-1] < MOTION_MIN_DURATION_S:
        raise ValueError(f"NS run covers only {elapsed[-1]:.2f} s; need at least "
                         f"{MOTION_MIN_DURATION_S:g} s for full-protocol annotation")
    if layout not in ("2panel", "3panel"):
        raise ValueError(f"Unsupported motion layout: {layout}")
    count = 3 if layout == "3panel" else 2
    fig, axes = plt.subplots(count, 1, figsize=(11.3, 7.35 if count == 3 else 5.65),
                             sharex=True, gridspec_kw={"height_ratios": ([1, 1, .82] if count == 3 else [1, 1])})
    handles = [Line2D([0], [0], color=BLUE, lw=1.4, label="CH1 IR AC"),
               Line2D([0], [0], color=RED, lw=1.4, label="CH2 IR AC")]
    handles += [Patch(facecolor=color, alpha=.24, edgecolor="none", label=legend)
                for _, _, _, legend, color in MOTION_WINDOWS]
    for ax, key, color, letter in ((axes[0], "ch1_ir", BLUE, "a"),
                                    (axes[1], "ch2_ir", RED, "b")):
        ax.plot(elapsed, ac[key], color=color, lw=.7, rasterized=True)
        ax.set_ylabel(f"CH{1 if letter == 'a' else 2} IR AC (counts)")
        ax.text(.009, .92, f"({letter})", transform=ax.transAxes, ha="left", va="top",
                weight="bold", fontsize=11, bbox=dict(facecolor="white", edgecolor="none", alpha=.87))
    rms1 = rms2 = None
    if count == 3:
        rms1 = ir_ac_rms(elapsed, ac["ch1_ir"])
        rms2 = ir_ac_rms(elapsed, ac["ch2_ir"])
        ax = axes[2]
        ax.plot(rms1.elapsed_s, rms1.ir_ac_rms_counts, color=BLUE, lw=1.4, label="CH1")
        ax.plot(rms2.elapsed_s, rms2.ir_ac_rms_counts, color=RED, lw=1.4, label="CH2")
        ax.set_ylabel("IR AC RMS (counts)")
        ax.text(.009, .92, "(c)", transform=ax.transAxes, ha="left", va="top",
                weight="bold", fontsize=11, bbox=dict(facecolor="white", edgecolor="none", alpha=.87))
    for ax in axes:
        for start, stop, _, _, color in MOTION_WINDOWS:
            ax.axvspan(start, stop, color=color, alpha=.15, lw=0, zorder=0)
            ax.axvline(start, color=color, alpha=.62, lw=.75, ls="--", zorder=1)
            ax.axvline(stop, color=color, alpha=.62, lw=.75, ls="--", zorder=1)
        ax.set_xlim(0, 300)
        ax.set_xticks(np.arange(0, 301, 30))
        format_axes(ax)
    axes[-1].set_xlabel("Elapsed time from first RAW sample (s)")
    legend = fig.legend(handles=handles, ncol=5, loc="upper center", bbox_to_anchor=(.5, .995),
                        fontsize=9, facecolor="white", framealpha=.96, edgecolor=BORDER,
                        columnspacing=1.3, handlelength=2)
    for label in legend.get_texts():
        label.set_fontweight("bold")
    fig.subplots_adjust(left=.122, right=.985, top=.90, bottom=.095,
                        hspace=.26 if count == 3 else .18)
    output = out / f"ControlledMotion_IR_AC_{layout}.{image_format}"
    save_kwargs = {"dpi": dpi, "facecolor": "white"}
    if image_format == "jpg":
        save_kwargs["pil_kwargs"] = {"quality": 96, "subsampling": 0}
    fig.savefig(output, **save_kwargs)
    plt.close(fig)
    if count == 3:
        rms = rms1.rename(columns={"ir_ac_rms_counts": "ch1_ir_ac_rms_counts"}).merge(
            rms2.rename(columns={"ir_ac_rms_counts": "ch2_ir_ac_rms_counts"}),
            on="elapsed_s", how="outer", validate="one_to_one")
        rms.to_csv(out / "ControlledMotion_IR_AC_RMS_2s.csv", index=False)
    return output



def plot_stationary_ir_ac(time, ac, out, dpi, image_format):
    """Full-run stationary visualization for S01/S02 using filtered IR AC traces only."""
    elapsed = np.asarray(time, dtype=float) - float(time[0])
    if elapsed[-1] < STATIONARY_MIN_DURATION_S:
        raise ValueError(f"Stationary run covers only {elapsed[-1]:.2f} s; need at least "
                         f"{STATIONARY_MIN_DURATION_S:g} s for full-run stationary figure")
    fig, axes = plt.subplots(2, 1, figsize=(11.3, 5.65), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1]})
    handles = []
    for ax, key, color, label, letter in ((axes[0], "ch1_ir", BLUE, "CH1 IR AC", "a"),
                                          (axes[1], "ch2_ir", RED, "CH2 IR AC", "b")):
        line, = ax.plot(elapsed, ac[key], color=color, lw=.7, rasterized=True, label=label)
        handles.append(line)
        ax.set_ylabel(f"{label.split()[0]} IR AC (counts)")
        ax.text(.009, .92, f"({letter})", transform=ax.transAxes, ha="left", va="top",
                weight="bold", fontsize=11, bbox=dict(facecolor="white", edgecolor="none", alpha=.87))
        ax.set_xlim(0, 300)
        ax.set_xticks(np.arange(0, 301, 30))
        format_axes(ax)
    axes[-1].set_xlabel("Elapsed time from first RAW sample (s)")
    legend = fig.legend(handles=handles, ncol=2, loc="upper center", bbox_to_anchor=(.5, .995),
                        fontsize=9, facecolor="white", framealpha=.96, edgecolor=BORDER,
                        columnspacing=1.3, handlelength=2)
    for label in legend.get_texts():
        label.set_fontweight("bold")
    fig.subplots_adjust(left=.122, right=.985, top=.90, bottom=.095, hspace=.18)
    output = out / f"Stationary_IR_AC_2panel.{image_format}"
    save_kwargs = {"dpi": dpi, "facecolor": "white"}
    if image_format == "jpg":
        save_kwargs["pil_kwargs"] = {"quality": 96, "subsampling": 0}
    fig.savefig(output, **save_kwargs)
    plt.close(fig)
    return output


def process_run(participant, run, path, member, args, output_root, left_limits, right_limits):
    time, raw = read_raw(path, member)
    start, end = max(args.start, float(time[0])), min(args.end, float(time[-1]))
    if end <= start:
        raise ValueError(f"Requested interval does not overlap recording {time[0]:.3f}–{time[-1]:.3f} s")
    smooth, ac = process(raw)    # Process the full run before time selection.
    dc = dc_windows(raw)         # Unmasked raw 8-s DC visualization estimate.
    source_label = Path(member).stem if member else path.stem
    if member and path.stem.upper() not in source_label.upper():
        source_label = path.stem + "_" + source_label
    source_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", source_label).strip("_.")
    # Distinguish identically named extracted/ZIP sources without overwriting outputs.
    source_fingerprint = hashlib.sha1((str(path.resolve()) + "::" + (member or "")).encode()).hexdigest()[:8]
    out = output_root / participant / run / f"{source_label}_{source_fingerprint}"
    out.mkdir(parents=True, exist_ok=True)
    images = []
    modes = ((False, True) if args.variants == "both" else
             (False,) if args.variants == "full" else (True,))
    for ch in (1, 2):
        for clean in modes:
            images.append(plot_channel(ch, time, raw, ac, dc, start, end, out, args.dpi,
                                       left_limits, right_limits, args.image_format,
                                       participant, run, clean=clean))
    table = pd.DataFrame({"time_s": time})
    for key in raw:
        table[key+"_raw"] = raw[key]
        table[key+"_hampel"] = smooth[key]
        table[key+"_ac_0p5_4hz"] = ac[key]
    table.to_csv(out / "PPG_samples.csv", index=False)
    dc.to_csv(out / "DC_window_levels.csv", index=False)
    motion_status = "not_applicable_stationary"
    stationary_status = "not_applicable_motion"
    if run in MOTION_RUNS and GENERATE_MOTION_FIGURE and not args.no_motion:
        try:
            motion = plot_controlled_motion(time, ac, out, args.dpi, args.image_format,
                                            layout=args.motion_layout)
            images.append(motion)
            motion_status = "scheduled_protocol_figure_created"
        except ValueError as exc:
            # The existing short-window plots are still useful for incomplete NS exports.
            motion_status = f"skipped: {exc}"
            print("  MOTION FIGURE SKIPPED:", exc, flush=True)
    elif run in MOTION_RUNS:
        motion_status = "disabled"

    if run in STATIONARY_RUNS and GENERATE_STATIONARY_FIGURE and not args.no_stationary:
        try:
            stationary = plot_stationary_ir_ac(time, ac, out, args.dpi, args.image_format)
            images.append(stationary)
            stationary_status = "full_run_stationary_figure_created"
        except ValueError as exc:
            stationary_status = f"skipped: {exc}"
            print("  STATIONARY FIGURE SKIPPED:", exc, flush=True)
    elif run in STATIONARY_RUNS:
        stationary_status = "disabled"
    for item in images:
        print("  Saved:", item, flush=True)
    return {"participant": participant, "run": run, "source_file": str(path),
            "zip_member": member or "", "samples": len(time), "plot_start_s": start,
            "plot_end_s": end, "output_directory": str(out),
            "figure_variants": args.variants, "figure_count": len(images),
            "motion_figure_status": motion_status,
            "stationary_figure_status": stationary_status,
            "motion_intervals": "scheduled_from_first_raw_sample_not_motion_detected",
            "status": "ok"}


def main(argv=None):
    args = arguments(argv)
    if (not np.isfinite([args.start, args.end]).all()
            or args.start < 0 or args.end <= args.start):
        raise ValueError("Require finite 0 <= --start < --end")
    if args.dpi < 72:
        raise ValueError("--dpi must be at least 72")
    left_limits = dict(LEFT_Y_LIMITS)
    right_limits = dict(RIGHT_Y_LIMITS)
    left_limits.update(axis_limit_overrides(args.left_ylim))
    right_limits.update(axis_limit_overrides(args.right_ylim))
    validate_axis_limits(left_limits, "LEFT_Y_LIMITS")
    validate_axis_limits(right_limits, "RIGHT_Y_LIMITS")
    jobs = select_sources(args)
    root = (args.out_dir or (Path(OUTPUT_FOLDER) if OUTPUT_FOLDER.strip() else
            Path(__file__).resolve().parent / "HEMOX_PPG_R2_Styled_outputs")).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    print(f"Using manuscript font: {setup_plot_style()} | DPI: {args.dpi}", flush=True)
    print(f"Found {len(jobs)} matching RAW export(s)", flush=True)
    record, failures = [], []
    for idx, (participant, run, path, member) in enumerate(jobs, 1):
        print(f"[{idx}/{len(jobs)}] {participant}/{run}: {path.name} / {Path(member).name if member else 'RAW.csv'}", flush=True)
        try:
            record.append(process_run(participant, run, path, member, args, root,
                                      left_limits, right_limits))
        except (ValueError, KeyError, OSError, BadZipFile) as exc:
            failures.append(f"{participant}/{run} {path}: {exc}")
            record.append({"participant": participant, "run": run, "source_file": str(path),
                           "zip_member": member or "", "status": f"ERROR: {exc}"})
            print("  FAILED:", exc, flush=True)
    pd.DataFrame(record).to_csv(root / "source_manifest.csv", index=False)
    print("Output folder:", root, flush=True)
    if failures:
        raise RuntimeError(f"{len(failures)} of {len(jobs)} export(s) failed; see source_manifest.csv. "
                           + "; ".join(failures))
    return root


if __name__ == "__main__":
    main()
