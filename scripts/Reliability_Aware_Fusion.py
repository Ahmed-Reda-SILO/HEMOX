from __future__ import annotations

import argparse
from pathlib import Path
from io import BytesIO
from zipfile import ZipFile
import re
import warnings
import gc

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.backends.backend_agg import FigureCanvasAgg

DURATION_S = 300.0
FROZEN_M1_OFFSET_PP = 0.906754  # P01/P02; never refit from the displayed run
DPI = 600
FALLBACK_DPIS = (450, 300)
SAVED_DPI_LOG = []
CURRENT_EXPORT_DPI = DPI  # session-level adaptive export DPI; downgraded after the first OOM
CH1, CH2, FUSED = '#0072B2', '#D55E00', '#009E73'  # Explorer palette
DIFF_COLOR, LIMIT_COLOR = '#7A3E9D', '#4D4D4D'
MOTION_BANDS = (
    (120.0, 130.0, 'RH', CH1),
    (180.0, 190.0, 'LH', CH2),
    (240.0, 250.0, 'Both', '#7A3E9D'),
)
FIELD_MAP = {
    'RAW': '_RAW.csv',
    'WINDOWS': '_WINDOWS.csv',
    'PARAMETER_TRENDS': '_PARAMETER_TRENDS.csv',
}


def parse_args():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description='HEMOX Explorer-compatible native signal-trends and optional IMU storyboard')
    p.add_argument('--data-dir', type=Path, default=here, help='Folder containing P01 ... P10 or participant ZIP files')
    p.add_argument('--archive', type=Path, help='Optional explicit participant ZIP, useful when several versions exist')
    p.add_argument('--participant', default='P09', help='P01 ... P10 (default: P09)')
    p.add_argument('--run', default='NS02', help='NS01 or NS02 (default: NS02)')
    p.add_argument('--fusion-source', choices=('explorer', 'parameter'), default='explorer',
                   help='explorer: reconstruct display from WINDOWS; parameter: use exported PARAMETER_TRENDS')
    p.add_argument('--spo2-display', choices=('native','m1'), default='native',
                   help='native (default) preserves the original reliability-aware fusion figure; m1 is retained only as an explicit optional illustration')
    p.add_argument('--spo2-offset', type=float, default=FROZEN_M1_OFFSET_PP,
                   help='Frozen P01/P02 M1 offset in pp; only used with --spo2-display m1')
    p.add_argument('--no-spo2-calibration', action='store_true',
                   help='Backward-compatible alias for --spo2-display native')
    p.add_argument('--only-explorer', action='store_true',
                   help='Output only the four Explorer-style charts plus one combined 2x2 plot')
    p.add_argument('--no-clean-variants', action='store_true',
                   help='Skip supplementary no-shadow, all-solid variants (normally generated)')
    p.add_argument('--duration', type=float, default=DURATION_S)
    p.add_argument('--imu-max-points', type=int, default=6500,
                   help='Max approximate plotted IMU samples; min/max-preserving reduction, full values saved to CSV')
    p.add_argument('--out-dir', type=Path,
                   help='Output folder; default Study_Reliability_Aware_Fusion_outputs/Pxx_NSxx next to code')
    return p.parse_args()


def canonical_subject(v):
    m = re.fullmatch(r'[Pp]\s*0*(\d{1,2})', str(v).strip())
    if not m or not (1 <= int(m.group(1)) <= 10):
        raise ValueError(f'Invalid participant {v!r}; expected P01 to P10')
    return f'P{int(m.group(1)):02d}'


def canonical_run(v):
    run = str(v).strip().upper()
    if run not in ('NS01', 'NS02'):
        raise ValueError('Representative controlled-motion run must be NS01 or NS02')
    return run


def matches_run(path, subject, run):
    # Find the complete subject/run token; ensures P04S01 is not misread as P04NS01.
    return bool(re.search(r'P0*' + str(int(subject[1:])) + re.escape(run) + r'(?=[_./\\]|$)',
                          str(path).upper()))


def resolve_source(root, subject, archive):
    if archive is not None:
        p = archive.expanduser().resolve()
        if not p.is_file() or p.suffix.lower() != '.zip':
            raise FileNotFoundError(f'--archive must point to an existing ZIP file: {p}')
        return 'zip', p
    root = root.expanduser().resolve()
    if root.is_file() and root.suffix.lower() == '.zip':
        return 'zip', root
    if not root.is_dir():
        raise FileNotFoundError(f'Data directory does not exist: {root}')
    if root.name.upper() == subject:
        return 'directory', root
    folders = [x for x in root.iterdir() if x.is_dir() and x.name.upper() == subject]
    if len(folders) > 1:
        raise ValueError(f'Duplicate participant folders: {folders}')
    if folders:
        return 'directory', folders[0]
    # Search exact 'P04.zip' first. Support versioned P10(1).zip, but refuse
    # ambiguously choosing among versions when both P10(1) and P10(3) exist.
    zips = [x for x in root.iterdir() if x.is_file() and x.suffix.lower() == '.zip'
            and re.fullmatch(re.escape(subject) + r'(?:\(\d+\))?', x.stem.upper())]
    exact = [x for x in zips if x.stem.upper() == subject]
    if len(exact) == 1:
        return 'zip', exact[0]
    if len(zips) == 1:
        return 'zip', zips[0]
    if len(zips) > 1:
        raise ValueError(f'Several archives for {subject}: {[x.name for x in zips]}; use --archive.')
    raise FileNotFoundError(f'No extracted {subject} folder or {subject}.zip found in {root}')


def resolve_members(kind, source, subject, run):
    if kind == 'directory':
        entries = [x for x in source.rglob('*.csv') if matches_run(x.name, subject, run)]
    else:
        with ZipFile(source) as z:
            entries = [x for x in z.namelist() if not x.endswith('/') and x.lower().endswith('.csv')
                       and matches_run(Path(x).name, subject, run)]
    found = {}
    for key, ending in FIELD_MAP.items():
        options = [x for x in entries if str(x).upper().endswith(ending.upper())]
        if len(options) > 1:
            raise ValueError(f'{subject}/{run}: duplicate {ending} files in source: {options}')
        found[key] = options[0] if options else None
    if not (found['WINDOWS'] or found['PARAMETER_TRENDS']):
        raise FileNotFoundError(f'{subject}/{run}: need WINDOWS.csv or PARAMETER_TRENDS.csv')
    return found


def load_csv(kind, source, member):
    if kind == 'directory':
        return pd.read_csv(member, low_memory=False)
    with ZipFile(source) as z:
        return pd.read_csv(BytesIO(z.read(member)), low_memory=False)


def numeric(df, col, *, required=True):
    if col not in df.columns:
        if required:
            raise ValueError(f'Missing expected export column: {col}')
        return pd.Series(np.full(len(df), np.nan), index=df.index)
    return pd.to_numeric(df[col], errors='coerce')


def valid_time(df, col, duration):
    t = numeric(df, col)
    if t.isna().any() or not t.is_monotonic_increasing or t.duplicated().any():
        raise ValueError(f'Invalid/non-monotonic {col} in HEMOX export')
    mask = t.ge(0) & t.le(duration)
    return df.loc[mask].copy(), t.loc[mask].to_numpy(float)


def bool_column(df, col):
    if col not in df.columns:
        raise ValueError(f'Missing eligibility column: {col}')
    return df[col].astype(str).str.strip().str.lower().isin(('true','yes','1','valid','accepted')).to_numpy(bool)


def weighted(a, b, w1, w2):
    a, b = np.asarray(a, float), np.asarray(b, float)
    good1 = np.isfinite(a) & np.isfinite(w1) & (w1 > 0)
    good2 = np.isfinite(b) & np.isfinite(w2) & (w2 > 0)
    output = np.full(len(a), np.nan)
    both = good1 & good2
    output[both] = (a[both] * w1[both] + b[both] * w2[both]) / (w1[both] + w2[both])
    only1 = good1 & ~good2
    only2 = good2 & ~good1
    output[only1], output[only2] = a[only1], b[only2]
    return output


def explorer_windows_trends(windows, duration):
    """Mirrors Explorer v1.7.12 v172AndroidWindows display reconstruction.

    This produces a visualization from analysis windows, not a firmware re-run.
    CH1/CH2 displayed values remain the exported filtered window values, as in
    renderSignalEvidence; eligibility/status rules gate the fused estimate.
    """
    w, t = valid_time(windows, 'end_s', duration)
    accepted1 = w['ch1_status'].astype(str).str.strip().str.lower().eq('accepted').to_numpy(bool)
    accepted2 = w['ch2_status'].astype(str).str.strip().str.lower().eq('accepted').to_numpy(bool)
    eligible1 = bool_column(w, 'ch1_spo2_eligible')
    eligible2 = bool_column(w, 'ch2_spo2_eligible')
    sqi1, sqi2 = numeric(w,'ch1_sqi').to_numpy(float), numeric(w,'ch2_sqi').to_numpy(float)
    mov1, mov2 = numeric(w,'ch1_motion_percent').to_numpy(float), numeric(w,'ch2_motion_percent').to_numpy(float)
    w1 = np.where(accepted1 & np.isfinite(sqi1) & np.isfinite(mov1),
                  np.maximum(0, sqi1) * np.maximum(0, 1-mov1/100), 0.0)
    w2 = np.where(accepted2 & np.isfinite(sqi2) & np.isfinite(mov2),
                  np.maximum(0, sqi2) * np.maximum(0, 1-mov2/100), 0.0)
    h1, h2 = numeric(w,'ch1_hr_bpm').to_numpy(float), numeric(w,'ch2_hr_bpm').to_numpy(float)
    o1, o2 = numeric(w,'ch1_spo2_percent').to_numpy(float), numeric(w,'ch2_spo2_percent').to_numpy(float)
    ghr1, ghr2 = np.where(accepted1, h1, np.nan), np.where(accepted2, h2, np.nan)
    go1, go2 = np.where(eligible1, o1, np.nan), np.where(eligible2, o2, np.nan)
    hr = weighted(ghr1, ghr2, w1, w2)
    spo2 = weighted(go1, go2, w1, w2)
    # Explorer v1.7.12 display rule: both numerically eligible and difference >3 -> no SpO2 fusion.
    # Note: this is the browser's display reconstruction, NOT a claim about firmware fallback.
    mismatch = np.isfinite(go1) & np.isfinite(go2) & (np.abs(go1-go2) > 3)
    spo2[mismatch] = np.nan
    return pd.DataFrame(dict(time_s=t,ch1_hr_bpm=h1,ch2_hr_bpm=h2,
                             fused_hr_native_bpm=hr, ch1_spo2_percent=o1,
                             ch2_spo2_percent=o2,fused_spo2_native_percent=spo2,
                             ch1_status=w['ch1_status'].to_numpy(),ch2_status=w['ch2_status'].to_numpy(),
                             ch1_spo2_eligible=eligible1,ch2_spo2_eligible=eligible2,
                             ch1_weight=w1,ch2_weight=w2,interchannel_spo2_disagreement=mismatch))


def exported_parameter_trends(parameter, duration):
    p,t = valid_time(parameter, 'end_s', duration)
    cols = ('ch1_hr_bpm','ch2_hr_bpm','fused_hr_bpm','ch1_spo2_percent',
            'ch2_spo2_percent','fused_spo2_percent')
    for c in cols:
        numeric(p,c)
    return pd.DataFrame(dict(time_s=t,
        ch1_hr_bpm=numeric(p,'ch1_hr_bpm').to_numpy(float),
        ch2_hr_bpm=numeric(p,'ch2_hr_bpm').to_numpy(float),
        fused_hr_native_bpm=numeric(p,'fused_hr_bpm').to_numpy(float),
        ch1_spo2_percent=numeric(p,'ch1_spo2_percent').to_numpy(float),
        ch2_spo2_percent=numeric(p,'ch2_spo2_percent').to_numpy(float),
        fused_spo2_native_percent=numeric(p,'fused_spo2_percent').to_numpy(float)))


def make_imu(raw, duration):
    # Explorer reads common_time_us, converts microseconds to seconds; gravity stays
    # in accelerometer magnitude. No per-run reference-device offset is applied.
    r, t = valid_time(raw,'common_time_us',duration * 1e6)
    t = t / 1e6
    result = pd.DataFrame(dict(time_s=t))
    for ch in (1,2):
        a_cols = [f'ch{ch}_acc_{v}_g' for v in 'xyz']
        g_cols = [f'ch{ch}_gyro_{v}_dps' for v in 'xyz']
        a = np.column_stack([numeric(r,col).to_numpy(float) for col in a_cols])
        gyro = np.column_stack([numeric(r,col).to_numpy(float) for col in g_cols])
        medians = []
        for j in range(3):
            finite = gyro[np.isfinite(gyro[:,j]),j]
            if len(finite) == 0:
                raise ValueError(f'CH{ch} gyroscope axis {j} contains no finite values')
            medians.append(float(np.median(finite)))
        result[f'ch{ch}_acc_magnitude_g'] = np.sqrt(np.sum(a*a,axis=1))
        result[f'ch{ch}_gyro_magnitude_dps'] = np.sqrt(np.sum((gyro-np.asarray(medians))**2,axis=1))
        for name, value in zip('xyz', medians):
            result[f'ch{ch}_gyro_median_{name}_dps'] = value
    return result


def select_minmax_indices(table, ycols, maxpts):
    """Min/max-preserving display reduction; full-rate IMU is always exported separately."""
    n = len(table)
    if n <= maxpts:
        return np.arange(n)
    bins = np.array_split(np.arange(n),max(1,maxpts//10))
    candidates = {0,n-1}
    for b in bins:
        if not len(b):
            continue
        candidates.update((int(b[0]),int(b[-1])))
        for col in ycols:
            y=table[col].to_numpy(float)[b]
            if np.isfinite(y).any():
                candidates.add(int(b[np.nanargmin(y)]))
                candidates.add(int(b[np.nanargmax(y)]))
    return np.asarray(sorted(candidates), dtype=int)


def setup_plot_style():
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ('Arial Narrow','Arial','Liberation Sans','DejaVu Sans') if f in installed),'DejaVu Sans')
    plt.rcParams.update({
        'font.family':font,'font.weight':'bold','axes.labelweight':'bold',
        'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,
        'legend.fontsize':10.5,'axes.linewidth':1.45,'lines.linewidth':1.6,
        'axes.spines.top':False,'axes.spines.right':False,
        'savefig.dpi':DPI,'figure.dpi':120,
    })
    return font


def shade(ax, first_panel=False):
    for start,end,label,color in MOTION_BANDS:
        ax.axvspan(start,end,facecolor=color,alpha=.095,edgecolor='none',zorder=0)
        # The three movement periods are identified in the manuscript caption.
        # Avoid internal labels that could obscure transient peaks or legends.


def breaks_at_gaps(x,y,gap):
    x=np.asarray(x,float); y=np.asarray(y,float)
    if len(x)<2:
        return x,y
    # Plotting NaNs prevents visually connecting a real missing time interval.
    large=np.flatnonzero(np.diff(x)>gap)
    if not len(large):
        return x,y
    return np.insert(x,large+1,np.nan),np.insert(y,large+1,np.nan)


def draw_explorer_fused(ax, x, native_y, *, color=FUSED, label='Fused', lw=2.2,
                        return_diagnostics=False, solid_hold=False):
    """Explorer fused rendering. Default dashed display hold matches Explorer;
    optional solid_hold changes presentation ONLY, never eligibility/data.
    Returns actual segment counts for QA when requested.
    """
    xx=np.asarray(x,float); yy=np.asarray(native_y,float)
    previous=None; holding=False; last_missing=None
    counts={'solid':0,'held_dashed':0,'resume_vertical':0}
    first_label=True
    for t,y in zip(xx,yy):
        if not np.isfinite(t):
            continue
        if np.isfinite(y):
            current=(float(t),float(y))
            if previous is not None:
                if holding:
                    ax.plot((previous[0],current[0]),(previous[1],previous[1]),
                            color=color,lw=max(1,lw*.8),alpha=.52,
                            ls='-' if solid_hold else (0,(6,4)),
                            zorder=3,label='_nolegend_')
                    ax.plot((current[0],current[0]),(previous[1],current[1]),
                            color=color,lw=lw,alpha=.96,zorder=4,label='_nolegend_')
                    counts['held_dashed']+=1
                    counts['resume_vertical']+=1
                else:
                    ax.plot((previous[0],current[0]),(previous[1],current[1]),
                            color=color,lw=lw,alpha=.96,zorder=4,
                            label=label if first_label else '_nolegend_')
                    counts['solid']+=1
                    first_label=False
            if previous is None:
                ax.plot([current[0]],[current[1]],marker='.',ms=3,color=color,
                        label=label if first_label else '_nolegend_',zorder=4)
                first_label=False
            previous=current;holding=False;last_missing=None
        elif previous is not None:
            holding=True;last_missing=float(t)
    if holding and previous is not None and last_missing is not None:
        ax.plot((previous[0],last_missing),(previous[1],previous[1]),
                color=color,lw=max(1,lw*.8),alpha=.52,
                ls='-' if solid_hold else (0,(6,4)),
                zorder=3,label='_nolegend_')
        counts['held_dashed']+=1
    return counts if return_diagnostics else None


def plot_panel(ax, panel, trends, imu_plot, duration, first_panel, letter=None,
               show_differences=True, add_legend=True, clean=False):
    if not clean:
        shade(ax,first_panel)
    ax2 = None
    if panel in ('Acceleration','Gyroscope'):
        isacc=panel=='Acceleration'
        field='acc_magnitude_g' if isacc else 'gyro_magnitude_dps'
        ylabel=('Acceleration (g)' if isacc else 'Gyroscope (°/s)') if letter else (
            'Acceleration magnitude (g)' if isacc else 'Gyroscope magnitude (°/s)')
        for ch,col in ((1,CH1),(2,CH2)):
            x,y=breaks_at_gaps(imu_plot.time_s,imu_plot[f'ch{ch}_{field}'],1.0)
            ax.plot(x,y,color=col,lw=1.05,alpha=.86,label=f'CH{ch}',zorder=2)
        ax.set_ylabel(ylabel)
    else:
        ishr=panel=='HR'
        cols=(('ch1_hr_bpm','ch2_hr_bpm','fused_hr_selected') if ishr else
              ('ch1_spo2_percent','ch2_spo2_percent','fused_spo2_selected_percent'))
        fused_label = 'Fused'
        for key,name,col,width in zip(cols, ('CH1','CH2',fused_label),
                                      (CH1,CH2,FUSED),(1.55,1.55,2.2)):
            if name.startswith('Fused'):
                draw_explorer_fused(ax,trends.time_s,trends[key],color=col,label=name,
                                    lw=width,solid_hold=clean)
            else:
                x,y=breaks_at_gaps(trends.time_s,trends[key],4.25)
                ax.plot(x,y,'-',color=col,lw=width,alpha=.83,
                        marker='o',markersize=1.9,label=name,zorder=3)
        if ishr:
            ax.set_ylabel('HR (bpm)')
            diff_col='abs_hr_diff_bpm'
            diff_label='|ΔHR|'
            diff_ylabel=r'$|\Delta \mathrm{HR}|$ (bpm)'
            threshold=5.0
        else:
            ax.set_ylabel(r'$\mathbf{SpO}_{\mathbf{2}}$ (%)')
            diff_col='abs_spo2_diff_pp'
            diff_label=r'$|\Delta \mathrm{SpO}_{2}|$'
            diff_ylabel=r'$|\Delta \mathrm{SpO}_{2}|$ (%)'
            threshold=3.0
        if show_differences:
            ax2 = ax.twinx()
            x2,y2=breaks_at_gaps(trends.time_s,trends[diff_col],4.25)
            ax2.plot(x2,y2,'-',color=DIFF_COLOR,lw=1.35,alpha=.9,label=diff_label,zorder=2)
            ax2.axhline(threshold,color=LIMIT_COLOR,lw=1.1,ls='-' if clean else '--',alpha=.9,
                        label=f'Threshold ({threshold:g})',zorder=1)
            ax2.set_ylabel(diff_ylabel,color=DIFF_COLOR)
            ax2.tick_params(axis='y',colors=DIFF_COLOR)
            ax2.spines['right'].set_visible(True)
            ax2.spines['right'].set_color(DIFF_COLOR)
            finite_diff = pd.to_numeric(trends[diff_col], errors='coerce').to_numpy(float)
            finite_diff = finite_diff[np.isfinite(finite_diff)]
            ymax = max(threshold * 1.25, np.nanmax(finite_diff) * 1.08 if finite_diff.size else threshold * 1.25)
            ax2.set_ylim(0, ymax)
    ax.set_xlim(0,duration)
    ax.set_xticks(np.arange(0,duration+1,30))
    ax.grid(True,color='#CBD1D7',lw=.62,alpha=.72,zorder=0)
    ax.set_axisbelow(True)
    if letter:
        ax.text(.012,.97,f'({letter})',transform=ax.transAxes,va='top',ha='left',
                fontsize=13,weight='bold')
    handles,labels=ax.get_legend_handles_labels()
    if ax2 is not None:
        h2,l2=ax2.get_legend_handles_labels()
        handles += h2
        labels += l2
    if add_legend:
        # Standalone figures: place the entire legend above the axes, NOT on data.
        leg=ax.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,1.035),
                      ncol=min(len(labels),3),borderaxespad=0,
                      framealpha=.97,facecolor='white',edgecolor='#A9B3BD',
                      fontsize=9.4,handlelength=1.75,columnspacing=.85,
                      borderpad=.35,labelspacing=.3)
        for txt in leg.get_texts():
            txt.set_fontweight('bold')
    for t in (*ax.get_xticklabels(),*ax.get_yticklabels()):
        t.set_fontweight('bold')
    if ax2 is not None:
        for t in ax2.get_yticklabels():
            t.set_fontweight('bold')


def _fresh_agg_canvas(fig):
    """Drop any cached Agg renderer before a lower-DPI retry."""
    try:
        old = getattr(fig, 'canvas', None)
        if old is not None and hasattr(old, 'renderer'):
            old.renderer = None
    except Exception:
        pass
    FigureCanvasAgg(fig)
    gc.collect()


def _release_figure(fig):
    """Aggressively release artists/canvas after export on RAM-constrained PCs."""
    try:
        plt.close(fig)
    finally:
        try:
            fig.clear()
        except Exception:
            pass
        try:
            fig.canvas = None
        except Exception:
            pass
        gc.collect()


def save_figure(fig, path):
    """Save JPG using a session-level adaptive DPI.

    Start at 600 DPI. Once an allocation fails, permanently downgrade this Python
    run to the next safe DPI (450, then 300) instead of retrying 600 DPI for every
    subsequent figure. This prevents repeated failed RendererAgg allocations from
    fragmenting/exhausting process memory.
    """
    global CURRENT_EXPORT_DPI

    path = Path(path)
    ladder = (DPI, *FALLBACK_DPIS)
    try:
        start_idx = ladder.index(CURRENT_EXPORT_DPI)
    except ValueError:
        start_idx = 0

    last_error = None
    for dpi in ladder[start_idx:]:
        try:
            if path.exists():
                try:
                    path.unlink()
                except OSError:
                    pass

            # A failed high-DPI draw can leave a large Agg renderer cached on the
            # Figure. Always attach a fresh canvas before an export attempt.
            _fresh_agg_canvas(fig)
            fig.savefig(
                path,
                dpi=dpi,
                facecolor='white',
                format='jpg',
                pil_kwargs={'quality':92, 'subsampling':1},
            )

            SAVED_DPI_LOG.append(int(dpi))
            CURRENT_EXPORT_DPI = int(dpi)
            _release_figure(fig)
            return dpi

        except MemoryError as exc:
            last_error = exc
            try:
                if path.exists():
                    path.unlink()
            except OSError:
                pass

            next_idx = ladder.index(dpi) + 1
            if next_idx < len(ladder):
                new_dpi = ladder[next_idx]
                if CURRENT_EXPORT_DPI != new_dpi:
                    warnings.warn(
                        f'Raster allocation failed at {dpi} DPI; all remaining figures '
                        f'will use {new_dpi} DPI for this run to avoid repeated OOM.',
                        RuntimeWarning,
                        stacklevel=2,
                    )
                CURRENT_EXPORT_DPI = int(new_dpi)
                _fresh_agg_canvas(fig)
                continue
            break

    _release_figure(fig)
    raise MemoryError(
        f'Could not save {path.name} even at {FALLBACK_DPIS[-1]} DPI. '
        'The process is genuinely RAM-constrained. Re-run after closing memory-heavy '
        'applications, or use --only-explorer / --no-clean-variants.'
    ) from last_error


def _explorer_series(ax,trends,kind,with_letter=None,clean=False):
    """Signal-filtered and fused trends as Explorer: four SEPARATE chart concepts."""
    if kind in ('SpO2','HR'):
        ishr=kind=='HR'
        ch1='ch1_hr_bpm' if ishr else 'ch1_spo2_percent'
        ch2='ch2_hr_bpm' if ishr else 'ch2_spo2_percent'
        native='fused_hr_selected' if ishr else 'fused_spo2_selected_percent'
        for col,color,label in ((ch1,CH1,'CH1'),(ch2,CH2,'CH2')):
            # Explorer breaks the line on a missing CHANNEL value; no visual hold.
            ax.plot(trends.time_s,trends[col],color=color,lw=1.65,alpha=.88,label=label)
        draw_explorer_fused(ax,trends.time_s,trends[native],color=FUSED,label='Fused',
                            lw=2.3,solid_hold=clean)
        ax.set_ylabel('HR (bpm)' if ishr else r'$\mathbf{SpO}_{\mathbf{2}}$ (%)')
        hs=[Line2D([0],[0],color=c,lw=w,label=l) for c,w,l in
            ((CH1,1.65,'CH1'),(CH2,1.65,'CH2'),(FUSED,2.3,'Fused'))]
    else:
        spo2=kind=='AbsSpO2'
        col='abs_spo2_diff_pp' if spo2 else 'abs_hr_diff_bpm'
        limit=3 if spo2 else 5
        color=DIFF_COLOR if spo2 else '#B54328'
        ax.plot(trends.time_s,trends[col],color=color,lw=1.8,label='|CH1 − CH2|')
        ax.axhline(limit,color=DIFF_COLOR,lw=1.3,ls='-' if clean else '--',
                   label=f'{limit} threshold')
        ax.set_ylabel(r'$|\Delta \mathrm{SpO}_{2}|$ (pp)' if spo2 else r'$|\Delta \mathrm{HR}|$ (bpm)')
        ax.set_ylim(bottom=0)
        hs=[Line2D([0],[0],color=color,lw=1.8,label='|CH1 − CH2|'),
            Line2D([0],[0],color=DIFF_COLOR,lw=1.3,
                   ls='-' if clean else '--',label=f'{limit} limit')]
    ax.grid(True,color='#CBD1D7',lw=.62,alpha=.72)
    ax.set_axisbelow(True)
    ax.set_xlim(0,float(trends.attrs.get('duration_s',DURATION_S)))
    ax.set_xticks(np.arange(0,float(trends.attrs.get('duration_s',DURATION_S))+1,30))
    if with_letter:
        ax.text(.012,.97,f'({with_letter})',transform=ax.transAxes,va='top',ha='left',
                fontsize=12,weight='bold')
    return hs


def make_explorer_signal_figures(trends,args,figs,clean=False):
    """The actual Explorer section: four separate charts and one compact 2x2."""
    names=('SpO2','HR','AbsSpO2','AbsHR')
    stems=('Explorer_SpO2_Trends','Explorer_HR_Trends',
           'Explorer_Absolute_SpO2_Difference','Explorer_Absolute_HR_Difference')
    created=[]
    for name,stem in zip(names,stems):
        fig,ax=plt.subplots(figsize=(7.2,3.8))
        hs=_explorer_series(ax,trends,name,clean=clean)
        ax.set_xlabel('Time (s)')
        leg=fig.legend(handles=hs,loc='upper center',bbox_to_anchor=(.5,.989),
                       ncol=len(hs),fontsize=9.7,handlelength=1.85,
                       columnspacing=1.1,framealpha=.97,edgecolor='#A9B3BD')
        for t in leg.get_texts():t.set_fontweight('bold')
        fig.subplots_adjust(left=.15,right=.975,bottom=.215,top=.81)
        out=figs/f'{args.participant}_{args.run}_{stem}.jpg'
        save_figure(fig,out);created.append(out)
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.4),sharex=True)
    for ax,name,letter in zip(axs.flat,names,'abcd'):
        _explorer_series(ax,trends,name,letter,clean=clean)
        ax.set_xlabel('Time (s)')
        ax.tick_params(labelsize=10.5)
        ax.yaxis.label.set_size(11)
    handles=[Line2D([0],[0],color=c,lw=w,label=l) for c,w,l in
             ((CH1,1.6,'CH1'),(CH2,1.6,'CH2'),(FUSED,2.3,'Fused'))]
    handles.append(Line2D([0],[0],color=DIFF_COLOR,lw=1.3,
                          ls='-' if clean else '--',label='Difference limit'))
    leg=fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.997),
                   ncol=3,fontsize=9,framealpha=.97,edgecolor='#A9B3BD')
    for t in leg.get_texts():t.set_fontweight('bold')
    fig.subplots_adjust(left=.13,right=.98,bottom=.085,top=.86,wspace=.39,hspace=.40)
    out=figs/f'{args.participant}_{args.run}_Explorer_Signal_Evidence_4Panel.jpg'
    save_figure(fig,out);created.append(out)
    return created


def make_figures(trends,imu,args,figs,clean=False, include_hold_in_shared_legend=True):
    """Save 8 JPGs; combined uses one shared outside legend, standalone outside legends.

    The acceleration and gyroscope panels are shared by both figure variants.
    include_hold_in_shared_legend=False creates an additional clean combined
    figure style matching the user's requested layout: no motion shading and no
    hold entry in the shared legend; every figure now omits that entry.
    Display-only holds remain visually drawn but are not valid measurements.
    """
    fields=('ch1_acc_magnitude_g','ch2_acc_magnitude_g',
            'ch1_gyro_magnitude_dps','ch2_gyro_magnitude_dps')
    indices=select_minmax_indices(imu,fields,args.imu_max_points)
    imu_plot=imu.iloc[indices].copy()
    panels=(('Acceleration','a','Acceleration'),('Gyroscope','b','Gyroscope'),
            ('HR','c','HR'),('SpO2','d','SpO2'))
    generated=[]

    for show_differences in (True,False):
        fig,axes=plt.subplots(4,1,figsize=(7.2,10.3),sharex=True)
        for i,(name,letter,_) in enumerate(panels):
            plot_panel(axes[i],name,trends,imu_plot,args.duration,i==0,letter,
                       show_differences=show_differences,add_legend=False,clean=clean)
            if i<3:
                axes[i].tick_params(axis='x',labelbottom=False)
        axes[-1].set_xlabel('Time (s)')
        # One shared legend at the very top. Never place legend inside any plot.
        shared=[
            Line2D([0],[0],color=CH1,lw=1.6,label='CH1'),
            Line2D([0],[0],color=CH2,lw=1.6,label='CH2'),
            Line2D([0],[0],color=FUSED,lw=2.2,label='Fused'),
        ]
        if show_differences:
            shared.extend([
                Line2D([0],[0],color=DIFF_COLOR,lw=1.35,
                       label='Inter-channel |difference|'),
                Line2D([0],[0],color=LIMIT_COLOR,lw=1.1,
                       ls='-' if clean else '--',label='Difference guide'),
            ])
        shared_ncol = min(3, max(1, len(shared)))
        lg=fig.legend(handles=shared,loc='upper center',bbox_to_anchor=(.5,.994),
                      ncol=shared_ncol,framealpha=.97,facecolor='white',edgecolor='#A9B3BD',
                      fontsize=9.3,handlelength=1.85,columnspacing=1.25,
                      borderpad=.32,labelspacing=.33)
        for txt in lg.get_texts():txt.set_fontweight('bold')
        fig.subplots_adjust(left=.17,right=.85 if show_differences else .985,
                            bottom=.064,top=.885 if show_differences else .915,hspace=.23)
        suffix='' if show_differences else '_NoDiff'
        if not include_hold_in_shared_legend:
            suffix = suffix + '_NoHoldLegend'
        combined=figs/f'{args.participant}_{args.run}_Reliability_Aware_Fusion_4Panel{suffix}.jpg'
        save_figure(fig,combined)
        generated.append(combined)

    # Individual figures: legend is explicitly above and outside each axis.
    for name,letter,stem in panels:
        for show_differences in ((True,False) if name in ('HR','SpO2') else (True,)):
            fig,ax=plt.subplots(figsize=(7.2,3.8))
            plot_panel(ax,name,trends,imu_plot,args.duration,name=='Acceleration',
                       show_differences=show_differences,add_legend=True,clean=clean)
            ax.set_xlabel('Time (s)')
            two_rows = show_differences and name in ('HR','SpO2')
            fig.subplots_adjust(left=.17,right=.85 if two_rows else .985,
                                bottom=.225,top=.745 if two_rows else .82)
            suffix='' if show_differences else '_NoDiff'
            out=figs/f'{args.participant}_{args.run}_{stem}{suffix}.jpg'
            save_figure(fig,out)
            generated.append(out)
    return generated,len(imu_plot)


def summarize_motion_intervals(imu,trends):
    intervals=(('stationary_1',0,120),('right_hand',120,130),('stationary_2',130,180),
               ('left_hand',180,190),('stationary_3',190,240),('bilateral',240,250),
               ('stationary_4',250,300))
    rows=[]
    for name,start,end in intervals:
        q=imu[(imu.time_s>=start)&(imu.time_s<end)]
        p=trends[(trends.time_s>=start)&(trends.time_s<end)]
        row=dict(condition=name,start_s=start,end_s=end,imu_samples=len(q),trend_windows=len(p))
        for ch in (1,2):
            for suffix in ('acc_magnitude_g','gyro_magnitude_dps'):
                v=q[f'ch{ch}_{suffix}']
                row[f'ch{ch}_{suffix}_max']=float(v.max()) if v.notna().any() else np.nan
            for suffix in ('hr_bpm','spo2_percent'):
                if suffix=='hr_bpm':
                    col=f'ch{ch}_hr_bpm'
                else:
                    col=f'ch{ch}_spo2_percent'
                row[f'ch{ch}_{suffix}_valid_windows']=int(p[col].notna().sum())
        for col in ('fused_hr_selected','fused_spo2_selected_percent'):
            row[f'{col}_valid_windows']=int(p[col].notna().sum())
        rows.append(row)
    return pd.DataFrame(rows)


def report_motion_peak_deviations(trends):
    """Report maximum inter-channel deviations during the controlled-motion periods.

    Deviation is the existing Explorer quantity |CH1-CH2|, computed from finite
    filtered channel estimates irrespective of fusion eligibility.  Motion windows
    follow the protocol: RH 120-130 s, LH 180-190 s, and Both 240-250 s.
    The interval end is treated as exclusive, matching summarize_motion_intervals().
    """
    motion_intervals=(
        ('RH',120.0,130.0),
        ('LH',180.0,190.0),
        ('Both',240.0,250.0),
    )
    rows=[]
    for label,start,end in motion_intervals:
        p=trends[(trends.time_s>=start)&(trends.time_s<end)]
        row=dict(motion=label,start_s=start,end_s=end,trend_windows=len(p))
        for metric,col,unit in (('spo2','abs_spo2_diff_pp','pp'),
                                ('hr','abs_hr_diff_bpm','bpm')):
            vals=pd.to_numeric(p[col],errors='coerce')
            finite=vals[np.isfinite(vals.to_numpy(float))]
            if len(finite):
                idx=finite.idxmax()
                row[f'max_{metric}_deviation']=float(finite.loc[idx])
                row[f'max_{metric}_time_s']=float(trends.loc[idx,'time_s'])
            else:
                row[f'max_{metric}_deviation']=np.nan
                row[f'max_{metric}_time_s']=np.nan
            row[f'{metric}_unit']=unit
        rows.append(row)

    report=pd.DataFrame(rows)

    # Overall peak across all controlled-motion periods only.
    motion_mask=np.zeros(len(trends),dtype=bool)
    t=pd.to_numeric(trends['time_s'],errors='coerce').to_numpy(float)
    for _,start,end in motion_intervals:
        motion_mask |= (t>=start) & (t<end)
    p_all=trends.loc[motion_mask]
    overall=dict(motion='Overall motion',start_s=np.nan,end_s=np.nan,trend_windows=len(p_all))
    for metric,col,unit in (('spo2','abs_spo2_diff_pp','pp'),
                            ('hr','abs_hr_diff_bpm','bpm')):
        vals=pd.to_numeric(p_all[col],errors='coerce')
        finite=vals[np.isfinite(vals.to_numpy(float))]
        if len(finite):
            idx=finite.idxmax()
            overall[f'max_{metric}_deviation']=float(finite.loc[idx])
            overall[f'max_{metric}_time_s']=float(trends.loc[idx,'time_s'])
        else:
            overall[f'max_{metric}_deviation']=np.nan
            overall[f'max_{metric}_time_s']=np.nan
        overall[f'{metric}_unit']=unit
    report=pd.concat([report,pd.DataFrame([overall])],ignore_index=True)
    return report


def main():
    args=parse_args()
    args.participant=canonical_subject(args.participant)
    args.run=canonical_run(args.run)
    if args.duration<=0 or args.imu_max_points<100:
        raise ValueError('Duration must be positive; imu-max-points must be >=100')
    if not np.isfinite(args.spo2_offset):
        raise ValueError('--spo2-offset must be finite')
    font=setup_plot_style()
    kind,source=resolve_source(args.data_dir,args.participant,args.archive)
    members=resolve_members(kind,source,args.participant,args.run)
    raw=load_csv(kind,source,members['RAW']) if members['RAW'] else None
    windows=load_csv(kind,source,members['WINDOWS']) if members['WINDOWS'] else None
    parameter=load_csv(kind,source,members['PARAMETER_TRENDS']) if members['PARAMETER_TRENDS'] else None
    if args.fusion_source=='explorer' and windows is not None:
        trends=explorer_windows_trends(windows,args.duration)
        source_statement='Explorer v1.7.12: WINDOWS native CH + final-window fused display reconstruction'
    elif parameter is not None:
        trends=exported_parameter_trends(parameter,args.duration)
        source_statement='Explorer PARAMETER_TRENDS fallback: native app-exported values'
    else:
        raise FileNotFoundError('--fusion-source parameter requires PARAMETER_TRENDS.csv')
    if trends.empty:
        raise ValueError('Selected source contains no parameter windows within requested duration')
    if raw is None:
        if not args.only_explorer:
            raise FileNotFoundError('RAW.csv is required for the original four-panel IMU storyboard. '
                                    'Use --only-explorer to plot the Signal-filtered and fused trends without RAW.')
        imu=pd.DataFrame(columns=['time_s','ch1_acc_magnitude_g','ch2_acc_magnitude_g',
                                  'ch1_gyro_magnitude_dps','ch2_gyro_magnitude_dps'])
    else:
        imu=make_imu(raw,args.duration)
        if imu.empty and not args.only_explorer:
            raise ValueError('Selected RAW file has no samples in the requested duration')
    # Manuscript model separation:
    #   HR final model: M0 (native).
    #   SpO2 fusion FIGURE: native M0, so the green fused trace remains the
    #   actual reliability-weighted fusion of the displayed CH1/CH2 estimates.
    #   SpO2 validation/performance model: M1 offset-only, exported separately.
    trends['fused_hr_selected']=trends['fused_hr_native_bpm']
    trends['fused_spo2_m0_percent']=trends['fused_spo2_native_percent']
    trends['fused_spo2_m1_percent']=(
        trends['fused_spo2_native_percent']-args.spo2_offset
    ).clip(0,100)

    args.spo2_display='native' if args.no_spo2_calibration else args.spo2_display
    trends.attrs['spo2_native_display']=(args.spo2_display=='native')
    if args.spo2_display=='native':
        trends['fused_spo2_selected_percent']=trends['fused_spo2_m0_percent']
        applied_model='M0 native reliability-aware fusion display'
        offset=0.0
    else:
        # Explicit OPTIONAL illustration only. This does not redefine the
        # reliability-aware fusion itself; M1 is a downstream calibration.
        trends['fused_spo2_selected_percent']=trends['fused_spo2_m1_percent']
        applied_model='OPTIONAL M1 calibrated illustration (downstream of fusion)'
        offset=args.spo2_offset

    validation_model='M1 offset-only manuscript model; frozen P01/P02 development offset'
    validation_offset=args.spo2_offset
    # Hold is visual-only. ORIGINAL finite selected values are never overwritten.
    # The line renderer decides solid vs lighter dashed missing stretches using
    # selected values; these ffilled columns are provided ONLY for audit/CSV.
    trends['fused_hr_display_held_only']=trends['fused_hr_selected'].ffill()
    trends['fused_spo2_display_held_only_percent']=trends['fused_spo2_selected_percent'].ffill()
    trends['fused_spo2_m1_display_held_only_percent']=trends['fused_spo2_m1_percent'].ffill()
    trends['hr_display_held_flag']=trends['fused_hr_selected'].isna() & trends['fused_hr_display_held_only'].notna()
    trends['spo2_display_held_flag']=trends['fused_spo2_selected_percent'].isna() & trends['fused_spo2_display_held_only_percent'].notna()
    # Explorer absolute differences use FINITE filtered CH values, irrespective of
    # whether one/both were eligible for fused selection.
    trends['abs_hr_diff_bpm']=np.abs(pd.to_numeric(trends['ch1_hr_bpm'], errors='coerce') - pd.to_numeric(trends['ch2_hr_bpm'], errors='coerce'))
    trends['abs_spo2_diff_pp']=np.abs(pd.to_numeric(trends['ch1_spo2_percent'], errors='coerce') - pd.to_numeric(trends['ch2_spo2_percent'], errors='coerce'))
    trends.attrs['duration_s']=args.duration
    out=(args.out_dir if args.out_dir else
         Path(__file__).resolve().parent/'Study_Reliability_Aware_Fusion_outputs'/f'{args.participant}_{args.run}').resolve()
    figs,tables=out/'Figures',out/'Tables'
    figs.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    trends.to_csv(tables/'filtered_fused_figure_data.csv',index=False)
    imu.to_csv(tables/'raw_imu_magnitudes.csv',index=False)
    summary=summarize_motion_intervals(imu,trends)
    summary.to_csv(tables/'motion_interval_diagnostics.csv',index=False)
    motion_peaks=report_motion_peak_deviations(trends)
    motion_peaks.to_csv(tables/'motion_peak_deviations.csv',index=False)
    # Preserve the original 13 figures in Figures/ and add an independent clean
    # counterpart for each in Figures_NoShadow_Solid/. Both read identical data.
    explorer_figures=make_explorer_signal_figures(trends,args,figs,clean=False)
    if args.only_explorer:
        figures,display_count=[],0
    else:
        figures,display_count=make_figures(trends,imu,args,figs,clean=False)
    original_figures=explorer_figures+figures
    clean_figures=[]
    clean_nohold_figures=[]
    if not args.no_clean_variants:
        clean_dir=out/'Figures_NoShadow_Solid'
        clean_dir.mkdir(parents=True,exist_ok=True)
        clean_figures=make_explorer_signal_figures(trends,args,clean_dir,clean=True)
        if not args.only_explorer:
            storyboard_clean,clean_display_count=make_figures(
                trends,imu,args,clean_dir,clean=True)
            clean_figures.extend(storyboard_clean)
            assert clean_display_count == display_count

            clean_nohold_dir=out/'Figures_NoShadow_NoHoldLegend'
            clean_nohold_dir.mkdir(parents=True,exist_ok=True)
            clean_nohold_figures,clean_nohold_display_count=make_figures(
                trends,imu,args,clean_nohold_dir,clean=True,include_hold_in_shared_legend=False)
            assert clean_nohold_display_count == display_count
    figures=original_figures+clean_figures+clean_nohold_figures
    meta=pd.DataFrame([dict(participant=args.participant,run=args.run,source_type=kind,source_file=str(source),
                            fusion_source=args.fusion_source,source_description=source_statement,
                            raw_file=str(members['RAW']),windows_file=str(members['WINDOWS']),
                            parameter_file=str(members['PARAMETER_TRENDS']),
                            spo2_fusion_plot_model=applied_model,
                            spo2_validation_model=validation_model,
                            spo2_m1_subtracted_offset_pp=validation_offset,
                            hr_model='M0 (unchanged)',source_imu_rows=len(imu),display_imu_rows=display_count,
                            source_trend_rows=len(trends),valid_fused_spo2_windows=int(trends.fused_spo2_selected_percent.notna().sum()),
                            valid_fused_hr_windows=int(trends.fused_hr_selected.notna().sum()),
                            font=font,figure_dpi_target=DPI,
                            figure_dpi_min_saved=(min(SAVED_DPI_LOG) if SAVED_DPI_LOG else np.nan),
                             figure_dpi_session_final=CURRENT_EXPORT_DPI)])
    meta.to_csv(tables/'plot_source_and_configuration.csv',index=False)
    findings=[f'# Explorer-style HEMOX trends: {args.participant}/{args.run}','',
              f'Source: {source}; visualization: {source_statement}.',
              f'Physiological windows: {len(trends)}; IMU samples: {len(imu)}.',
              'Explorer source rule: native filtered CH1/CH2 from WINDOWS if present, with an independently reconstructed window-display fused series; PARAMETER_TRENDS is the fallback.',
              f'SpO2 fusion figure: {applied_model}; displayed offset = {offset:+.6f} pp.',
              f'SpO2 validation/performance model: {validation_model}; M1 offset = {validation_offset:+.6f} pp.',
              f'Raster export target: {DPI} DPI; minimum actual saved DPI this run = {min(SAVED_DPI_LOG) if SAVED_DPI_LOG else "N/A"}.',
              'Model separation: HR remains M0. Reliability-aware fusion figures default to native fused SpO2 M0 so the fused curve matches the original fusion output. Final manuscript SpO2 validation/performance uses the separately exported M1 offset-only series.',
              'The original Explorer fused trace is solid between successive finite selected estimates, but lighter dashed and horizontal across missing estimates.',
              'ADDITIONAL Figures_NoShadow_Solid provides identical data in new no-motion-shadow, all-solid plots. The lighter horizontal hold remains DISPLAY-ONLY, despite being solid.',
              'Every legend contains channel/fusion labels only as applicable; the display-hold legend entry is omitted from ALL figures, without changing plotted data.',
              'Held values are visualization-only: `fused_*_selected` retain missingness and must be used for eligibility/analysis; `_display_held_only` must NOT be used as valid observations.',
              'CH1/CH2 filtered estimates remain plotted even in non-eligible windows, as in the Explorer. This does not imply they were eligible for fusion.',
              'The separate absolute inter-channel difference charts use raw finite filtered CH values, independent of eligibility, with visual guides 3 pp and 5 bpm.',
              'IMU: acceleration norm including gravity; gyroscope axes each median-subtracted before norm.',
              'Storyboard movement annotations (not included on the four exact-section charts): RH 120-130s; LH 180-190s; bilateral 240-250s.',
              'Motion-period peak deviations are reported in Tables/motion_peak_deviations.csv as the maximum finite inter-channel |CH1-CH2| SpO2 and HR differences and their exact trend-window times.',
              *[('Motion peak ' + str(r.motion) + ': max |dSpO2| = ' +
                  (f'{r.max_spo2_deviation:.3f} pp at {r.max_spo2_time_s:.1f} s' if np.isfinite(r.max_spo2_deviation) else 'N/A') +
                  '; max |dHR| = ' +
                  (f'{r.max_hr_deviation:.3f} bpm at {r.max_hr_time_s:.1f} s' if np.isfinite(r.max_hr_deviation) else 'N/A') + '.')
                 for r in motion_peaks.itertuples(index=False)],
              'WINDOWS reconstruction is visualization only, not the native FUSED_TRENDS validation analysis; do not calculate the manuscript RMSE from displayed or held values.',
              '',f'## Original figures ({len(original_figures)})',
              *[f'- Figures/{f.name}' for f in original_figures],
              '',f'## Additional no-shadow / solid figures ({len(clean_figures)})',
              *[f'- Figures_NoShadow_Solid/{f.name}' for f in clean_figures],
              '',f'## Additional no-shadow / no-hold-legend figures ({len(clean_nohold_figures)})',
              *[f'- Figures_NoShadow_NoHoldLegend/{f.name}' for f in clean_nohold_figures], '']
    (out/'Findings_reliability_aware_fusion.md').write_text('\n'.join(findings),encoding='utf-8')
    print(f'Participant={args.participant}, run={args.run}, fusion-source={args.fusion_source}')
    print(f'IMU rows={len(imu)}, plotted IMU rows={display_count}; windows={len(trends)}')
    print(f'Selected fused SpO2 windows={trends.fused_spo2_selected_percent.notna().sum()}; HR={trends.fused_hr_selected.notna().sum()}')
    print(f'Fusion-figure SpO2 display={args.spo2_display}; displayed offset={offset:.6f} pp; HR M0 unchanged.')
    print(f'Validation/performance SpO2 model=M1 offset-only; frozen offset={validation_offset:.6f} pp.')
    print('\nMotion-period maximum inter-channel deviations (|CH1-CH2|):')
    for r in motion_peaks.itertuples(index=False):
        spo2_txt=(f'{r.max_spo2_deviation:.3f} pp at {r.max_spo2_time_s:.1f} s'
                  if np.isfinite(r.max_spo2_deviation) else 'N/A')
        hr_txt=(f'{r.max_hr_deviation:.3f} bpm at {r.max_hr_time_s:.1f} s'
                if np.isfinite(r.max_hr_deviation) else 'N/A')
        print(f'  {r.motion:>14}: SpO2 {spo2_txt}; HR {hr_txt}')
    print(f'Motion peak table: {tables / "motion_peak_deviations.csv"}')
    if SAVED_DPI_LOG:
        print(f'Raster export target={DPI} DPI; minimum actual saved DPI={min(SAVED_DPI_LOG)} DPI; session final DPI={CURRENT_EXPORT_DPI}.')
    print(f'Original figures={len(original_figures)}; additional no-shadow/solid figures={len(clean_figures)}; no-hold-legend figures={len(clean_nohold_figures)}')
    for f in figures: print('Saved:',f)
    print('Output:',out)


if __name__=='__main__':
    main()
