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
SESSIONS = ('S01', 'S02', 'NS01', 'NS02')
COLORS = {'data': "#1E00FF", 'analysis': '#C62828', 'zero': '#929AA2'}
DURATION_S = 300.0
TOLERANCE_S = 1.0
DPI = 600


def parse_args():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description='N=10 pooled native-output Bland-Altman before calibration')
    p.add_argument('--data-dir', type=Path, default=here,
                   help='Directory with extracted P01 ... P10 folders (default: script directory)')
    p.add_argument('--out-dir', type=Path, default=here/'BA_before_outputs',
                   help='Results directory (default: BA_before_outputs beside script)')
    p.add_argument('--duration', type=float, default=DURATION_S)
    p.add_argument('--tolerance', type=float, default=TOLERANCE_S)
    p.add_argument('--offsets-csv', type=Path, default=None,
                   help='Optional columns: subject,run,offset_seconds,method,evidence; shift HEMOX elapsed time')
    p.add_argument('--allow-incomplete', action='store_true',
                   help='Analyze partial dataset but show missing/duplicate sessions in QC')
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


def performance(pairs, variable):
    ref_name, he_name, bound = (
        ('sp20_spo2','hemox_spo2',3.0) if variable == 'SpO2' else
        ('sp20_hr','hemox_hr',5.0))
    if pairs.empty:
        valid = pd.DataFrame(columns=[ref_name, he_name])
    else:
        valid = pairs[[ref_name,he_name]].apply(pd.to_numeric, errors='coerce').dropna()
        valid = valid[np.isfinite(valid[ref_name]) & np.isfinite(valid[he_name])]
    if valid.empty:
        return {'valid_pairs':0, 'bias':np.nan, 'sd_error':np.nan,
                'loa_low':np.nan, 'loa_high':np.nan, 'loa_width':np.nan,
                'mae':np.nan,'rmse':np.nan,'mape_pct':np.nan,
                'acc_100_minus_mape_pct':np.nan,'within_threshold_pct':np.nan}
    error = (valid[he_name]-valid[ref_name]).to_numpy(float)
    bias = float(np.mean(error))
    sd = float(np.std(error, ddof=1)) if len(error)>1 else np.nan
    ref = valid[ref_name].to_numpy(float)
    positive_ref = np.abs(ref)>1e-12
    mape = float(100*np.mean(np.abs(error[positive_ref]/ref[positive_ref]))) if positive_ref.any() else np.nan
    return {'valid_pairs':len(error), 'bias':bias, 'sd_error':sd,
            'loa_low':bias-1.96*sd, 'loa_high':bias+1.96*sd, 'loa_width':3.92*sd,
            'mae':float(np.mean(np.abs(error))), 'rmse':float(np.sqrt(np.mean(error**2))),
            'mape_pct':mape, 'acc_100_minus_mape_pct':100-mape,
            'within_threshold_pct':float(100*np.mean(np.abs(error)<=bound))}


def plot_pooled_ba(pairs, variable, output, show_legend=True):
    ref_name, he_name, unit = (
        ('sp20_spo2','hemox_spo2','pp') if variable == 'SpO2' else
        ('sp20_hr','hemox_hr','bpm'))
    data = pairs[[ref_name,he_name]].apply(pd.to_numeric,errors='coerce').dropna()
    data = data[np.isfinite(data[ref_name]) & np.isfinite(data[he_name])]
    if data.empty:
        raise ValueError(f'No finite paired observations for the {variable} plot')
    stats = performance(pairs, variable)
    means = ((data[he_name]+data[ref_name])/2).to_numpy(float)
    differences = (data[he_name]-data[ref_name]).to_numpy(float)
    fig, ax = plt.subplots(figsize=(7.2,3.8))
    # Blue: observations; red: numerical agreement lines. Line styles distinguish bias/LoA.
    points = ax.scatter(means,differences,s=17,c=COLORS['data'],alpha=.18,
               edgecolors='none',rasterized=True,
               label=f'Paired measurements: N = 10, n = {stats["valid_pairs"]:,}')
    bias_line = ax.axhline(stats['bias'],color=COLORS['analysis'],linestyle='-',linewidth=1.6,
               label=f'Bias = {stats["bias"]:+.2f} {unit}')
    loa_line = ax.axhline(stats['loa_low'],color=COLORS['analysis'],linestyle='--',linewidth=1.6,
               label=f'95% LoA = [{stats["loa_low"]:.2f}, {stats["loa_high"]:.2f}] {unit}')
    ax.axhline(stats['loa_high'],color=COLORS['analysis'],linestyle='--',linewidth=1.6)
    # SD is a numerical statistic, not an additional horizontal plot line.
    sd_label = Line2D([], [], color='none', linestyle='none',
                      label=f'SD = {stats["sd_error"]:.2f} {unit}')
    ax.axhline(0,color=COLORS['zero'],linestyle=':',linewidth=1.0,zorder=0)
    if variable == 'SpO2':
        ax.set_xlabel(r'Mean(HEMOX, SP-20) $\mathbf{SpO}_{\mathbf{2}}$ (%)')
        ax.set_ylabel(r'(HEMOX − SP-20) $\mathbf{SpO}_{\mathbf{2}}$ (%)')
    else:
        ax.set_xlabel('Mean(HEMOX, SP-20) HR (bpm)')
        ax.set_ylabel('(HEMOX-SP-20) HR (bpm)')
    ax.grid(True,color='#C8CDD2',linewidth=.6,alpha=.6)
    ax.set_axisbelow(True)
    x_min,x_max = float(np.min(means)),float(np.max(means))
    dx = max(x_max-x_min,.25 if variable=='SpO2' else 1.0)
    ax.set_xlim(x_min-.07*dx,x_max+.07*dx)
    y_min = min(float(np.min(differences)),stats['loa_low'],0)
    y_max = max(float(np.max(differences)),stats['loa_high'],0)
    dy = max(y_max-y_min,.5 if variable=='SpO2' else 1.0)
    ax.set_ylim(y_min-.12*dy,y_max+.12*dy)
    if show_legend:
        legend = ax.legend(handles=[points, bias_line, sd_label, loa_line], loc='best',
                           facecolor='white',edgecolor='#A9B3BD',framealpha=.97,
                           fontsize=10.5,handlelength=1.65,labelspacing=.35)
        for t in legend.get_texts():
            t.set_fontweight('bold')
    for t in (*ax.get_xticklabels(),*ax.get_yticklabels()):t.set_fontweight('bold')
    # No plot title, data-box caption, figure letter, or bottom explanatory footnote.
    fig.subplots_adjust(left=.17,right=.98,bottom=.24,top=.97)
    fig.savefig(output,dpi=DPI,facecolor='white',pil_kwargs={'quality':96,'subsampling':0})
    plt.close(fig)


def main():
    args = parse_args()
    if args.duration<=0 or args.tolerance<=0:
        raise ValueError('Duration and matching tolerance must be positive')
    root = args.data_dir.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f'Participant root not found: {root}')
    font = setup_plot_style()
    out = args.out_dir.resolve(); figs=out/'Figures'; tables=out/'Tables'
    figs.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    participants=find_participants(root)
    missing=[x for x in SUBJECTS if x not in participants]
    if missing and not args.allow_incomplete:
        raise FileNotFoundError(f'Missing participant folders: {missing}; expected P01 ... P10')
    offsets=load_offsets(args.offsets_csv)
    paired=[]; audit=[]; issues=[]
    for subject in SUBJECTS:
        if subject not in participants:
            continue
        found=index_sources(participants[subject],subject)
        for run in SESSIONS:
            hemox_paths,sp20_paths = found[run]['hemox'], found[run]['sp20']
            if len(hemox_paths)!=1 or len(sp20_paths)!=1:
                problem=(f'{subject}/{run}: expected exactly one FUSED_TRENDS and SP20 CSV; '
                         f'found {len(hemox_paths)}/{len(sp20_paths)}. HEMOX={hemox_paths} SP20={sp20_paths}')
                if not args.allow_incomplete:raise ValueError(problem)
                issues.append(problem)
                continue
            hemox=read_fused(hemox_paths[0]);sp20=read_sp20(sp20_paths[0])
            hemox=hemox[(hemox.hemox_t>=0)&(hemox.hemox_t<=args.duration)].copy()
            sp20=sp20[(sp20.ref_t>=0)&(sp20.ref_t<=args.duration)].copy()
            valid=hemox[hemox.quality_valid].copy()
            offset,method,evidence=offsets.get((subject,run),(0.,'provisional_zero_offset',''))
            valid['he_t_aligned']=valid.hemox_t+offset
            valid=valid[(valid.he_t_aligned>=0)&(valid.he_t_aligned<=args.duration)].copy()
            matches=match_explorer_style(sp20,valid,args.tolerance)
            for row in matches:
                row.update(subject=subject,run=run,offset_seconds=offset,
                           alignment_method=method,alignment_evidence=evidence,
                           spo2_error_pp=row['hemox_spo2']-row['sp20_spo2'],
                           hr_error_bpm=row['hemox_hr']-row['sp20_hr'])
            paired.extend(matches)
            spo2_n=int(sum(np.isfinite(r['hemox_spo2']) and np.isfinite(r['sp20_spo2']) for r in matches))
            hr_n=int(sum(np.isfinite(r['hemox_hr']) and np.isfinite(r['sp20_hr']) for r in matches))
            audit.append(dict(subject=subject,run=run,source_rows=len(hemox),
                 quality_valid_rows=len(valid),available_spo2=int(np.isfinite(valid.he_spo2).sum()),
                 available_hr=int(np.isfinite(valid.he_hr).sum()),sp20_reference_rows=len(sp20),
                 matched_time_pairs=len(matches),valid_spo2_pairs=spo2_n,valid_hr_pairs=hr_n,
                 mean_abs_dt_s=float(np.mean([abs(r['delta_t_s']) for r in matches])) if matches else np.nan,
                 max_abs_dt_s=float(np.max([abs(r['delta_t_s']) for r in matches])) if matches else np.nan,
                 offset_seconds=offset,alignment_method=method,alignment_evidence=evidence,
                 hemox_file=str(hemox_paths[0]),sp20_file=str(sp20_paths[0])))
            print(f'{subject}/{run}: FUSED={len(hemox)}, quality-valid={len(valid)}, '
                  f'matched={len(matches)}, SpO2={spo2_n}, HR={hr_n}, offset={offset:+g} s')
    if not paired:raise RuntimeError('No matched data. Check source files, time alignment, and tolerance.')
    x=pd.DataFrame(paired);qc=pd.DataFrame(audit)
    if (qc.duplicated(['subject','run']).any() or qc.subject.nunique()!=len(SUBJECTS) or len(qc)!=40) and not args.allow_incomplete:
        raise ValueError('Unexpected duplicate or missing sessions in N=10 / 40-run dataset.')
    x.to_csv(tables/'paired_measurements.csv',index=False)
    qc.to_csv(tables/'session_alignment_availability.csv',index=False)
    summaries=[]
    for variable in ('SpO2','HR'):
        m=performance(x,variable)
        available=int(qc.available_spo2.sum() if variable=='SpO2' else qc.available_hr.sum())
        rows=int(qc.source_rows.sum())
        summaries.append(dict(dataset='All P01-P10',variable=variable,
           unit='pp' if variable=='SpO2' else 'bpm',participants=int(qc.subject.nunique()),runs=len(qc),
           fused_source_rows=rows,quality_valid_rows=int(qc.quality_valid_rows.sum()),valid_outputs=available,
           output_availability_pct=100*available/rows if rows else np.nan,
           matched_time_pairs=len(x),**m))
        figure=figs/f'BA_N10_{variable}.jpg'
        plot_pooled_ba(x,variable,figure,show_legend=True)
        print(f'Saved: {figure.name}')
        figure_nl=figs/f'BA_N10_{variable}_NoLegend.jpg'
        plot_pooled_ba(x,variable,figure_nl,show_legend=False)
        print(f'Saved: {figure_nl.name}')
    summary=pd.DataFrame(summaries)
    summary.to_csv(tables/'overall_performance.csv',index=False)
    breakdown=[]
    for cols in (('subject',),('subject','run')):
        for key,subset in x.groupby(list(cols),sort=True):
            if not isinstance(key,tuple):key=(key,)
            entry=dict(zip(cols,key))
            entry['scope']='participant' if len(cols)==1 else 'run'
            for variable in ('SpO2','HR'):
                breakdown.append(dict(**entry,variable=variable,**performance(subset,variable)))
    pd.DataFrame(breakdown).to_csv(tables/'participant_and_run_performance.csv',index=False)
    counts=pd.DataFrame([dict(participants=int(qc.subject.nunique()),sessions=len(qc),
          recording_minutes=5*len(qc),source_rows=int(qc.source_rows.sum()),
          quality_valid_rows=int(qc.quality_valid_rows.sum()),matched_time_pairs=len(x),
          valid_spo2_pairs=int(summary.loc[summary.variable=='SpO2','valid_pairs'].iloc[0]),
          valid_hr_pairs=int(summary.loc[summary.variable=='HR','valid_pairs'].iloc[0]))])
    counts.to_csv(tables/'dataset_counts.csv',index=False)
    def line(row):
        return (f"n={int(row.valid_pairs):,}; Bias={row.bias:+.3f}; MAE={row.mae:.3f}; "
                f"RMSE={row.rmse:.3f}; ACC={row.acc_100_minus_mape_pct:.2f}%; "
                f"LoA=[{row.loa_low:.3f}, {row.loa_high:.3f}]; "
                f"within-threshold={row.within_threshold_pct:.2f}%; "
                f"availability={row.output_availability_pct:.2f}%")
    spo2=summary.loc[summary.variable=='SpO2'].iloc[0]
    hr=summary.loc[summary.variable=='HR'].iloc[0]
    text=["# HEMOX / SP-20 - pooled N=10 agreement before calibration", "",
          f"Subjects: {int(qc.subject.nunique())}; recordings: {len(qc)}; recording duration: {5*len(qc)} minutes.",
          f"Source: native FUSED_TRENDS (quality_valid=true); pairing tolerance: {args.tolerance:g} s; plot font: {font}.",
          "No calibration coefficients were fitted or applied. No calibration/validation split is used.", "",
          "## Overall performance", "",
          "SpO2 (pp): " + line(spo2),
          "HR (bpm): " + line(hr), "",
          "ACC = 100 - MAPE (%), as in Explorer. Separate columns report the percentage within +/-3 pp (SpO2) and +/-5 bpm (HR).",
          "Error sign: HEMOX minus SP-20. LoA = bias +/- 1.96 * sample SD.",
          "Plots contain no title or caption; blue markers denote matched samples and red lines denote bias and LoA.",
          "Two versions are saved for each variable: one with legend and one without legend, numerical annotations, titles or captions.", "",
          "## Interpretation limits", "",
          "Default matching uses a zero-offset elapsed-time assumption. This is provisional until independent synchronization evidence is supplied.",
          "Pooled paired-observation LoA are descriptive because each participant contributes repeated measurements.",
          "Comparator-based SpO2 accuracy against SP-20 is not an arterial-reference pulse-oximeter validation.",
          "Check output availability, per-run data coverage and range alongside paired-only accuracy."
    ]
    if issues:text.extend(['', '## Missing/duplicate sources', *issues])
    (out/'Findings_before_calibration.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(f'Complete. Participants={qc.subject.nunique()} runs={len(qc)}. Output: {out}')
    print(summary[['variable','valid_pairs','bias','mae','rmse','acc_100_minus_mape_pct','loa_low','loa_high']].to_string(index=False))


if __name__=='__main__':
    main()
