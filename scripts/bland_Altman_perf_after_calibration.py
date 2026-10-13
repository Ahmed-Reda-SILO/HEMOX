
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

SUBJECTS = tuple(f'P{i:02d}' for i in range(1, 11))
CALIBRATION_IDS = ('P01', 'P02')
VALIDATION_IDS = tuple(s for s in SUBJECTS if s not in CALIBRATION_IDS)
SESSIONS = ('S01', 'S02', 'NS01', 'NS02')
COLORS = {'data': '#1E00FF', 'analysis': '#C62828', 'zero': '#929AA2'}
DURATION_S = 300.0
TOLERANCE_S = 1.0
DPI = 600


def parse_args():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description='P01-P02 M1 SpO2 / P03-P10 held-out BA, native HR M0')
    p.add_argument('--data-dir', type=Path, default=here,
                   help='Root containing extracted P01 ... P10 folders (default: script directory)')
    p.add_argument('--out-dir', type=Path, default=here/'BA_after_calibration_outputs',
                   help='Output directory (default: BA_after_calibration_outputs beside script)')
    p.add_argument('--duration', type=float, default=DURATION_S)
    p.add_argument('--tolerance', type=float, default=TOLERANCE_S)
    p.add_argument('--offsets-csv', type=Path, default=None,
                   help='Optional columns: subject,run,offset_seconds,method,evidence; shift HEMOX elapsed time')
    return p.parse_args()


def setup_plot_style():
    installed = {f.name for f in font_manager.fontManager.ttflist}
    family = next((s for s in ('Arial Narrow','Arial','Liberation Sans','DejaVu Sans') if s in installed),
                  'DejaVu Sans')
    plt.rcParams.update({
        'font.family': family, 'font.weight': 'bold', 'axes.labelweight': 'bold',
        'axes.labelsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14,
        'legend.fontsize': 11, 'figure.figsize': (7.2, 3.8),
        'axes.linewidth': 1.6, 'lines.linewidth': 1.6,
        'axes.spines.top': False, 'axes.spines.right': False,
        'savefig.dpi': DPI, 'figure.dpi': 120,
    })
    return family


def canonical_subject(value):
    m = re.fullmatch(r'P\s*0*(\d{1,2})', str(value).strip(), re.I)
    return f'P{int(m.group(1)):02d}' if m and 1 <= int(m.group(1)) <= 10 else None


def session_from_path(path, subject):
    # Extract NS before S; a filename may contain HEMOXP01NS02 rather than an NS02 folder.
    for part in (path.name, *(parent.name for parent in path.parents)):
        u = part.upper()
        for match in re.finditer(r'P0*(\d{1,2})(NS0[12]|S0[12])', u):
            if f'P{int(match.group(1)):02d}' == subject:
                return match.group(2)
        if u in SESSIONS:
            return u
    return None


def find_participants(root):
    found = {}
    for folder in root.iterdir():
        if not folder.is_dir():
            continue
        subject = canonical_subject(folder.name)
        if subject:
            if subject in found:
                raise ValueError(f'Duplicate participant directory for {subject}: {found[subject]} / {folder}')
            found[subject] = folder
    return found


def index_sources(participant_dir, subject):
    sources = {run: {'hemox': [], 'sp20': []} for run in SESSIONS}
    for file in participant_dir.rglob('*.csv'):
        name = file.name.upper()
        if name.endswith('FUSED_TRENDS.CSV') and 'CALIBRATED' not in name:
            kind = 'hemox'
        elif name.startswith('SP20'):
            kind = 'sp20'
        else:
            continue
        run = session_from_path(file, subject)
        if run:
            sources[run][kind].append(file)
    return sources


def read_sp20(path):
    # Metadata first line, data headings second line.
    r = pd.read_csv(path, skiprows=1, skipinitialspace=True, dtype=str)
    r.columns = [str(x).strip().lstrip('#').strip() for x in r.columns]
    needed = ('Time', 'SpO2(%)', 'PR(bpm)')
    if not set(needed).issubset(r.columns):
        raise ValueError(f'SP20 columns absent in {path}: {set(needed)-set(r.columns)}')
    raw_times = r['Time'].astype(str).str.strip()
    stamps = pd.to_datetime(raw_times, errors='coerce', format='%m/%d/%Y %H:%M:%S')
    if stamps.isna().any():
        try:
            stamps = stamps.fillna(pd.to_datetime(raw_times, errors='coerce', format='mixed'))
        except (TypeError, ValueError):
            pass
    if stamps.isna().any() or len(stamps) == 0:
        raise ValueError(f'Invalid SP20 timestamps in {path}')
    if not stamps.is_monotonic_increasing or stamps.duplicated().any():
        raise ValueError(f'Duplicate/non-monotonic SP20 timestamps in {path}')
    return pd.DataFrame({
        'ref_t': (stamps-stamps.iloc[0]).dt.total_seconds().astype(float),
        'sp20_time': stamps.astype(str),
        'ref_spo2': pd.to_numeric(r['SpO2(%)'], errors='coerce'),
        'ref_hr': pd.to_numeric(r['PR(bpm)'], errors='coerce'),
    })


def read_fused(path):
    r = pd.read_csv(path)
    required = ('elapsed_s', 'fused_spo2_percent', 'fused_hr_bpm', 'quality_valid')
    if not set(required).issubset(r.columns):
        raise ValueError(f'FUSED_TRENDS columns absent in {path}: {set(required)-set(r.columns)}')
    times = pd.to_numeric(r['elapsed_s'], errors='coerce')
    if times.isna().any() or not times.is_monotonic_increasing or times.duplicated().any():
        raise ValueError(f'Invalid/non-monotonic FUSED_TRENDS elapsed_s in {path}')
    good = r['quality_valid'].astype(str).str.strip().str.lower().isin(('true','yes','valid','1'))
    return pd.DataFrame({
        'hemox_t': times, 'quality_valid': good,
        'he_spo2': pd.to_numeric(r['fused_spo2_percent'], errors='coerce'),
        'he_hr': pd.to_numeric(r['fused_hr_bpm'], errors='coerce'),
    })


def load_offsets(path):
    if path is None:
        return {}
    r = pd.read_csv(path)
    required = {'subject','run','offset_seconds','method','evidence'}
    if not required.issubset(r.columns):
        raise ValueError(f'Offsets CSV requires columns: {sorted(required)}')
    result = {}
    for row in r.itertuples(index=False):
        sub, run = canonical_subject(row.subject), str(row.run).upper().strip()
        if sub not in SUBJECTS or run not in SESSIONS or (sub,run) in result:
            raise ValueError(f'Invalid or duplicate offset entry: {row.subject}/{row.run}')
        method = str(row.method).strip().lower()
        evidence = str(row.evidence).strip()
        if method == 'verified' and (not evidence or evidence.lower() in ('nan','none')):
            raise ValueError(f'Verified alignment requires supporting evidence: {sub}/{run}')
        offset = float(row.offset_seconds)
        if not np.isfinite(offset):
            raise ValueError(f'Non-finite offset: {sub}/{run}')
        result[(sub,run)] = (offset, method, evidence)
    return result


def match_explorer_style(sp, he, tolerance):
    """Explorer's SP20-first greedy, one-to-one nearest time matcher, inclusive tolerance."""
    ref = list(sp.itertuples(index=False))
    hemox = list(he.itertuples(index=False))
    if not ref or not hemox:
        return []
    low = max(ref[0].ref_t, hemox[0].he_t_aligned)
    high = min(ref[-1].ref_t, hemox[-1].he_t_aligned)
    used, output = set(), []
    for sr in ref:
        if sr.ref_t < low or sr.ref_t > high:
            continue
        best_idx, best_dt = None, float('inf')
        for i, hr in enumerate(hemox):
            if i in used or hr.he_t_aligned < low or hr.he_t_aligned > high:
                continue
            delta = abs(hr.he_t_aligned - sr.ref_t)
            if delta < best_dt:
                best_idx, best_dt = i, delta
            if hr.he_t_aligned > sr.ref_t + tolerance and best_idx is not None:
                break
        if best_idx is not None and best_dt <= tolerance:
            used.add(best_idx)
            hr = hemox[best_idx]
            output.append({
                'sp20_elapsed_s': sr.ref_t, 'sp20_timestamp': sr.sp20_time,
                'hemox_elapsed_s': hr.hemox_t, 'hemox_aligned_s': hr.he_t_aligned,
                'delta_t_s': hr.he_t_aligned-sr.ref_t,
                'sp20_spo2': sr.ref_spo2, 'hemox_spo2': hr.he_spo2,
                'sp20_hr': sr.ref_hr, 'hemox_hr': hr.he_hr,
            })
    return output


def metric_values(prediction, reference, threshold):
    a = np.asarray(prediction, dtype=float)
    b = np.asarray(reference, dtype=float)
    good = np.isfinite(a) & np.isfinite(b)
    a, b = a[good], b[good]
    if len(a) == 0:
        raise ValueError('No valid paired observations for a requested metric')
    e = a-b
    bias = float(e.mean())
    sd = float(e.std(ddof=1)) if len(e)>1 else np.nan
    pos = abs(b) > 1e-12
    mape = float(100*np.mean(abs(e[pos]/b[pos]))) if pos.any() else np.nan
    return dict(valid_pairs=int(len(e)), bias=bias, sd_error=sd,
                loa_low=bias-1.96*sd, loa_high=bias+1.96*sd, loa_width=3.92*sd,
                mae=float(np.mean(abs(e))), rmse=float(np.sqrt(np.mean(e*e))),
                mape_pct=mape, acc_100_minus_mape_pct=100-mape,
                within_threshold_pct=float(100*np.mean(abs(e)<=threshold)))


def filter_valid(df, pred_col, ref_col):
    t = df[[pred_col,ref_col]].apply(pd.to_numeric, errors='coerce')
    return t.loc[np.isfinite(t[pred_col]) & np.isfinite(t[ref_col])].copy()


def plot_split_ba(data, variable, output, show_legend=True):
    if variable == 'SpO2':
        pcol, rcol, unit = 'selected_spo2', 'sp20_spo2', 'pp'
    elif variable == 'HR':
        pcol, rcol, unit = 'selected_hr', 'sp20_hr', 'bpm'
    else:
        raise ValueError(variable)
    calibration = filter_valid(data.loc[data.split.eq('Calibration')], pcol, rcol)
    validation = filter_valid(data.loc[data.split.eq('Validation')], pcol, rcol)
    if calibration.empty or validation.empty:
        raise ValueError(f'{variable}: both calibration and validation must have paired points')
    stats = metric_values(validation[pcol], validation[rcol], 3 if variable=='SpO2' else 5)
    xm_cal = ((calibration[pcol]+calibration[rcol])/2).to_numpy(float)
    yd_cal = (calibration[pcol]-calibration[rcol]).to_numpy(float)
    xm_val = ((validation[pcol]+validation[rcol])/2).to_numpy(float)
    yd_val = (validation[pcol]-validation[rcol]).to_numpy(float)
    fig, ax = plt.subplots(figsize=(7.2,3.8))
    # Draw calibration first and unobtrusively; validation stays visually prominent.
    h_cal = ax.scatter(xm_cal, yd_cal, s=16, c=COLORS['data'], alpha=.10,
                       edgecolors='none', rasterized=True, zorder=2,
                       label=f'Calibration: P01-P02 (n = {len(calibration):,})')
    h_val = ax.scatter(xm_val, yd_val, s=16, c=COLORS['data'], alpha=.80,
                       edgecolors='none', rasterized=True, zorder=3,
                       label=f'Validation: P03-P10 (n = {len(validation):,})')
    h_bias = ax.axhline(stats['bias'], color=COLORS['analysis'], lw=1.6,
                        linestyle='-', zorder=4,
                        label=f'Validation bias = {stats["bias"]:+.2f} {unit}')
    h_loa = ax.axhline(stats['loa_low'], color=COLORS['analysis'], lw=1.6,
                       linestyle='--', zorder=4,
                       label=f'Validation 95% LoA = [{stats["loa_low"]:.2f}, {stats["loa_high"]:.2f}] {unit}')
    ax.axhline(stats['loa_high'], color=COLORS['analysis'], lw=1.6, linestyle='--', zorder=4)
    h_sd = Line2D([], [], color='none', linestyle='none',
                  label=f'Validation SD = {stats["sd_error"]:.2f} {unit}')
    ax.axhline(0, color=COLORS['zero'], linestyle=':', linewidth=1, zorder=0)
    if variable=='SpO2':
        ax.set_xlabel(r'Mean(HEMOX, SP-20) $\mathbf{SpO}_{\mathbf{2}}$ (%)')
        ax.set_ylabel(r'(HEMOX - SP-20) $\mathbf{SpO}_{\mathbf{2}}$ (pp)')
    else:
        ax.set_xlabel('Mean(HEMOX, SP-20) HR (bpm)')
        ax.set_ylabel('(HEMOX - SP-20) HR (bpm)')
    ax.grid(True, color='#C8CDD2', linewidth=.6, alpha=.6)
    ax.set_axisbelow(True)
    allx = np.r_[xm_cal,xm_val]
    ally = np.r_[yd_cal,yd_val]
    xmin, xmax = float(allx.min()), float(allx.max())
    dx = max(xmax-xmin, .25 if variable=='SpO2' else 1)
    ax.set_xlim(xmin-.07*dx,xmax+.07*dx)
    ymin = min(float(ally.min()),float(stats['loa_low']),0)
    ymax = max(float(ally.max()),float(stats['loa_high']),0)
    dy = max(ymax-ymin, .5 if variable=='SpO2' else 1)
    ax.set_ylim(ymin-.12*dy,ymax+.12*dy)
    if show_legend:
        leg = ax.legend(handles=[h_cal,h_val,h_bias,h_sd,h_loa], loc='best',
                        facecolor='white',edgecolor='#A9B3BD',framealpha=.97,
                        fontsize=9.6,handlelength=1.65,labelspacing=.33)
        for text in leg.get_texts():
            text.set_fontweight('bold')
    for tick in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick.set_fontweight('bold')
    fig.subplots_adjust(left=.17,right=.98,bottom=.24,top=.97)
    fig.savefig(output, dpi=DPI, facecolor='white',
                pil_kwargs={'quality':96,'subsampling':0})
    plt.close(fig)
    return stats


def main():
    args = parse_args()
    if args.duration<=0 or args.tolerance<=0:
        raise ValueError('Duration and pairing tolerance must be positive')
    root = args.data_dir.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f'Data directory not found: {root}')
    out = args.out_dir.resolve()
    figs, tables = out/'Figures', out/'Tables'
    figs.mkdir(parents=True,exist_ok=True)
    tables.mkdir(parents=True,exist_ok=True)
    font = setup_plot_style()
    participants = find_participants(root)
    missing = [sub for sub in SUBJECTS if sub not in participants]
    if missing:
        raise FileNotFoundError(f'Expected P01-P10; missing participant folder(s): {missing}')
    offsets = load_offsets(args.offsets_csv)
    paired, audit = [], []
    for subject in SUBJECTS:
        indexed = index_sources(participants[subject],subject)
        split = 'Calibration' if subject in CALIBRATION_IDS else 'Validation'
        for run in SESSIONS:
            h_files, s_files = indexed[run]['hemox'], indexed[run]['sp20']
            if len(h_files)!=1 or len(s_files)!=1:
                raise ValueError(f'{subject}/{run}: expected 1 native FUSED_TRENDS and 1 SP20; '
                                 f'got {len(h_files)} and {len(s_files)}')
            he, sp = read_fused(h_files[0]), read_sp20(s_files[0])
            he = he[(he.hemox_t>=0)&(he.hemox_t<=args.duration)].copy()
            sp = sp[(sp.ref_t>=0)&(sp.ref_t<=args.duration)].copy()
            valid = he[he.quality_valid].copy()
            offset, method, evidence = offsets.get((subject,run),(0.,'provisional_zero_offset',''))
            valid['he_t_aligned']=valid.hemox_t+offset
            valid=valid[(valid.he_t_aligned>=0)&(valid.he_t_aligned<=args.duration)].copy()
            matches = match_explorer_style(sp,valid,args.tolerance)
            for row in matches:
                row.update(subject=subject,run=run,split=split,
                           alignment_offset_seconds=offset,alignment_method=method,
                           alignment_evidence=evidence)
            paired.extend(matches)
            spo_n=sum(np.isfinite(m['hemox_spo2']) and np.isfinite(m['sp20_spo2']) for m in matches)
            hr_n=sum(np.isfinite(m['hemox_hr']) and np.isfinite(m['sp20_hr']) for m in matches)
            audit.append(dict(subject=subject,run=run,split=split,
                              hemox_source=str(h_files[0]),sp20_source=str(s_files[0]),
                              source_rows=len(he),quality_valid_rows=len(valid),
                              matched_time_pairs=len(matches),matched_spo2_pairs=spo_n,
                              matched_hr_pairs=hr_n,alignment_method=method,
                              alignment_offset_seconds=offset,alignment_evidence=evidence))
            print(f'{subject}/{run} ({split}): paired={len(matches)}, SpO2={spo_n}, HR={hr_n}')
    x, qc = pd.DataFrame(paired), pd.DataFrame(audit)
    if x.empty or len(qc)!=40 or qc.subject.nunique()!=10 or qc.duplicated(['subject','run']).any():
        raise ValueError('Expected 40 unique runs with paired data from P01-P10')
    # Fit exactly once on P01-P02 ONLY, using native uncalibrated SpO2.
    cal = filter_valid(x.loc[x.split.eq('Calibration')],'hemox_spo2','sp20_spo2')
    if len(cal)<2:
        raise ValueError('Insufficient calibration SpO2 pairs')
    m1_offset = float((cal.hemox_spo2-cal.sp20_spo2).mean())
    x['selected_spo2']=np.clip(pd.to_numeric(x['hemox_spo2'],errors='coerce')-m1_offset,0,100)
    x['selected_hr']=pd.to_numeric(x['hemox_hr'],errors='coerce')
    x['spo2_before_error_pp']=x.hemox_spo2-x.sp20_spo2
    x['spo2_selected_error_pp']=x.selected_spo2-x.sp20_spo2
    x['hr_selected_error_bpm']=x.selected_hr-x.sp20_hr
    x.to_csv(tables/'paired_measurements_selected_models.csv',index=False)
    qc.to_csv(tables/'session_alignment_availability.csv',index=False)
    coeff = pd.DataFrame([dict(variable='SpO2',model='M1',calibration_ids='P01;P02',
                               calibration_participants=2,calibration_runs=8,
                               calibration_valid_pairs=len(cal),
                               mean_native_minus_reference_pp=m1_offset,
                               prediction_formula='clip(native_spo2 - offset, 0, 100)'),
                          dict(variable='HR',model='M0',calibration_ids='none',
                               calibration_participants=0,calibration_runs=0,
                               calibration_valid_pairs=0,
                               mean_native_minus_reference_pp=np.nan,
                               prediction_formula='native_hr (no additional calibration)')])
    coeff.to_csv(tables/'selected_model_coefficients.csv',index=False)
    selected_summary=[]
    for split in ('Calibration','Validation'):
        subset=x.loc[x.split.eq(split)]
        for variable,pcol,rcol,threshold,model in (
            ('SpO2','selected_spo2','sp20_spo2',3.0,'M1'),
            ('HR','selected_hr','sp20_hr',5.0,'M0')):
            s = metric_values(subset[pcol],subset[rcol],threshold)
            selected_summary.append(dict(split=split,variable=variable,selected_model=model,
                                         participants=int(subset.subject.nunique()),
                                         runs=int(subset[['subject','run']].drop_duplicates().shape[0]),
                                         **s))
    selected_summary = pd.DataFrame(selected_summary)
    selected_summary.to_csv(tables/'selected_models_by_split.csv',index=False)
    before_after=[]
    validation=x.loc[x.split.eq('Validation')]
    for variable,pred_col,ref_col,model,threshold in (
        ('SpO2','hemox_spo2','sp20_spo2','M0',3.0),
        ('SpO2','selected_spo2','sp20_spo2','M1',3.0),
        ('HR','selected_hr','sp20_hr','M0',5.0)):
        before_after.append(dict(variable=variable,model=model,split='Validation',
                                 **metric_values(validation[pred_col],validation[ref_col],threshold)))
    before_after=pd.DataFrame(before_after)
    before_after.to_csv(tables/'heldout_validation_before_after.csv',index=False)
    for variable in ('SpO2','HR'):
        filename=figs/f'BA_after_{variable}.jpg'
        plot_split_ba(x,variable,filename,show_legend=True)
        plot_split_ba(x,variable,figs/f'BA_after_{variable}_NoLegend.jpg',show_legend=False)
        print('Saved:',filename.name,'and NoLegend version')
    baseline=before_after[(before_after.variable=='SpO2')&(before_after.model=='M0')].iloc[0]
    after=before_after[(before_after.variable=='SpO2')&(before_after.model=='M1')].iloc[0]
    hr=before_after[(before_after.variable=='HR')&(before_after.model=='M0')].iloc[0]
    reduction=100*(baseline.rmse-after.rmse)/baseline.rmse
    findings=[
        '# HEMOX: after model selection, calibration versus held-out Bland-Altman', '',
        'P01-P02 calibration (20%); P03-P10 held-out validation (80%).',
        f'Calibration M1 offset = {m1_offset:+.6f} percentage points, estimated solely from {len(cal):,} valid P01-P02 pairs.',
        'Selected SpO2 = clip(native fused SpO2 - M1 offset, 0, 100) for both displayed cohorts.',
        'Selected HR = native fused HR (M0) for both displayed cohorts.',
        'Calibration points are faded, validation points opaque; ALL reported bias, SD and LoA lines come solely from validation.',
        f'Paired reference time tolerance = {args.tolerance:g} s; default zero device offset is PROVISIONAL unless externally verified.',
        f'Plot font: {font}; 600 DPI JPEG; no titles or captions.', '',
        '## Validation-only numerical findings', '',
        f'SpO2 M1: n={int(after.valid_pairs):,}, bias={after.bias:+.3f} pp, SD={after.sd_error:.3f} pp, '
        f'95% pooled LoA=[{after.loa_low:.3f}, {after.loa_high:.3f}] pp, MAE={after.mae:.3f} pp, RMSE={after.rmse:.3f} pp.',
        f'SpO2 M0 on the SAME validation observations: RMSE={baseline.rmse:.3f} pp; M1 matched RMSE reduction={reduction:.2f}%.',
        f'HR M0: n={int(hr.valid_pairs):,}, bias={hr.bias:+.3f} bpm, SD={hr.sd_error:.3f} bpm, '
        f'95% pooled LoA=[{hr.loa_low:.3f}, {hr.loa_high:.3f}] bpm, MAE={hr.mae:.3f} bpm, RMSE={hr.rmse:.3f} bpm.',
        '', '## Proposed manuscript paragraph', '',
        f'Bland-Altman analysis of the eight held-out participants (Fig. X) yielded a SpO2 bias of {after.bias:+.2f}% '
        f'and descriptive 95% LoA of [{after.loa_low:.2f}, {after.loa_high:.2f}]% following M1 calibration, '
        f'with a {reduction:.1f}% reduction in RMSE relative to the matched uncalibrated validation observations. '
        f'The native HR estimates (M0) showed a bias of {hr.bias:+.2f} bpm and descriptive 95% LoA of '
        f'[{hr.loa_low:.2f}, {hr.loa_high:.2f}] bpm. Calibration observations are displayed faintly for context '
        'and do not contribute to validation agreement statistics.', '',
        '## Methodological limitations', '',
        '- The SP20-relative zero-offset alignment must be independently verified before confirmation claims.',
        '- Pooled paired-point LoA are descriptive because within-participant observations are correlated.',
        '- SP20 comparator results are not equivalent to arterial-reference pulse-oximeter accuracy validation.',
        '- Calibration points are in-sample displays; no combined-cohort agreement line or pooled validation claim is made.',
    ]
    (out/'Findings_after_calibration.md').write_text('\n'.join(findings)+'\n',encoding='utf-8')
    print('M1 SpO2 offset (native minus SP20):',f'{m1_offset:.6f}')
    print(before_after[['variable','model','valid_pairs','bias','sd_error','loa_low','loa_high','mae','rmse']].to_string(index=False))
    print('Completed:',out)


if __name__=='__main__':
    main()
