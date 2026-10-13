
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

SUBJECTS = tuple(f'P{i:02d}' for i in range(1, 11))
SESSIONS = ('S01', 'S02', 'NS01', 'NS02')
CALIBRATION_SIZES = (0, 1, 2, 3, 4, 5)
DURATION_S = 300.0
TOLERANCE_S = 1.0
DPI = 600
FIGSIZE = (7.2, 3.8)
COLORS = {
    'BeforeBaseline': '#C62828',   # red
    'M1': '#0072B2',               # blue
    'M2': '#2E7D32',               # green
    'M3': '#7B1FA2',               # purple
}
MODEL_ORDER = ('M1', 'M2', 'M3')


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description='N=10 participant-disjoint multi-model calibration study')
    parser.add_argument('--data-dir', type=Path, default=here,
                        help='Root containing extracted P01 ... P10 folders (default: script folder)')
    parser.add_argument('--out-dir', type=Path, default=here / 'Study_perf_Multi_calibration_outputs',
                        help='Directory for plots and tables')
    parser.add_argument('--duration', type=float, default=DURATION_S)
    parser.add_argument('--tolerance', type=float, default=TOLERANCE_S)
    parser.add_argument('--offsets-csv', type=Path, default=None,
                        help='Optional columns: subject,run,offset_seconds,method,evidence')
    parser.add_argument('--baseline-csv', type=Path, default=here / 'overall_performance.csv',
                        help=('Optional fixed full-N=10 uncalibrated baseline CSV. If the file is absent, '
                              'the baseline is derived directly from the current native paired measurements.'))
    return parser.parse_args()


def setup_plot_style():
    installed = {f.name for f in font_manager.fontManager.ttflist}
    family = next((s for s in ('Arial Narrow', 'Arial', 'Liberation Sans', 'DejaVu Sans')
                   if s in installed), 'DejaVu Sans')
    plt.rcParams.update({
        'font.family': family, 'font.weight': 'bold', 'axes.labelweight': 'bold',
        'axes.labelsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14,
        'legend.fontsize': 10.5, 'figure.figsize': FIGSIZE,
        'axes.linewidth': 1.6, 'lines.linewidth': 1.6,
        'axes.spines.top': False, 'axes.spines.right': False,
        'savefig.dpi': DPI, 'figure.dpi': 120,
    })
    return family


def canonical_subject(value):
    m = re.fullmatch(r'P\s*0*(\d{1,2})', str(value).strip(), re.I)
    return f'P{int(m.group(1)):02d}' if m and 1 <= int(m.group(1)) <= 10 else None


def session_from_path(path, subject):
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
        if folder.is_dir():
            subject = canonical_subject(folder.name)
            if subject:
                if subject in found:
                    raise ValueError(f'Duplicate participant folder {subject}: {found[subject]}, {folder}')
                found[subject] = folder
    missing = sorted(set(SUBJECTS) - set(found))
    if missing:
        raise FileNotFoundError(f'Missing participant folders: {missing}; expected P01 ... P10')
    return found


def index_sources(folder, subject):
    output = {run: {'hemox': [], 'sp20': []} for run in SESSIONS}
    for path in folder.rglob('*.csv'):
        name = path.name.upper()
        if name.endswith('FUSED_TRENDS.CSV') and 'CALIBRATED' not in name:
            kind = 'hemox'
        elif name.startswith('SP20'):
            kind = 'sp20'
        else:
            continue
        run = session_from_path(path, subject)
        if run:
            output[run][kind].append(path)
    for run in SESSIONS:
        for kind in ('hemox', 'sp20'):
            if len(output[run][kind]) != 1:
                raise ValueError(f'{subject}/{run}: expected exactly one {kind} CSV; found {output[run][kind]}')
    return output


def read_sp20(path):
    r = pd.read_csv(path, skiprows=1, skipinitialspace=True, dtype=str)
    r.columns = [str(c).strip().lstrip('#').strip() for c in r.columns]
    needed = ('Time', 'SpO2(%)', 'PR(bpm)')
    if not set(needed).issubset(r.columns):
        raise ValueError(f'Missing SP20 fields in {path}: {set(needed)-set(r.columns)}')
    raw = r['Time'].astype(str).str.strip()
    times = pd.to_datetime(raw, errors='coerce', format='%m/%d/%Y %H:%M:%S')
    if times.isna().any():
        try:
            times = times.fillna(pd.to_datetime(raw, errors='coerce', format='mixed'))
        except (TypeError, ValueError):
            pass
    if not len(times) or times.isna().any() or not times.is_monotonic_increasing or times.duplicated().any():
        raise ValueError(f'Invalid/non-monotonic SP20 timestamps in {path}')
    return pd.DataFrame({
        'ref_t': (times-times.iloc[0]).dt.total_seconds().astype(float),
        'sp20_time': times.astype(str),
        'ref_spo2': pd.to_numeric(r['SpO2(%)'], errors='coerce'),
        'ref_hr': pd.to_numeric(r['PR(bpm)'], errors='coerce'),
    })


def read_fused(path):
    r = pd.read_csv(path)
    needed = ('elapsed_s', 'fused_spo2_percent', 'fused_hr_bpm', 'quality_valid')
    if not set(needed).issubset(r.columns):
        raise ValueError(f'Missing native fused fields in {path}: {set(needed)-set(r.columns)}')
    t = pd.to_numeric(r['elapsed_s'], errors='coerce')
    if t.isna().any() or not t.is_monotonic_increasing or t.duplicated().any():
        raise ValueError(f'Invalid/non-monotonic elapsed_s in {path}')
    good = r['quality_valid'].astype(str).str.strip().str.lower().isin(('true', 'yes', 'valid', '1'))
    return pd.DataFrame({
        'hemox_t': t, 'quality_valid': good,
        'he_spo2': pd.to_numeric(r['fused_spo2_percent'], errors='coerce'),
        'he_hr': pd.to_numeric(r['fused_hr_bpm'], errors='coerce'),
    })


def load_offsets(path):
    if path is None:
        return {}
    df = pd.read_csv(path)
    required = {'subject', 'run', 'offset_seconds', 'method', 'evidence'}
    if not required.issubset(df.columns):
        raise ValueError(f'Offsets CSV missing: {sorted(required-set(df.columns))}')
    if df.duplicated(['subject', 'run']).any():
        raise ValueError('Duplicate subject/run keys in offsets CSV')
    result = {}
    for row in df.itertuples(index=False):
        subject, run = canonical_subject(row.subject), str(row.run).strip().upper()
        if subject not in SUBJECTS or run not in SESSIONS:
            raise ValueError(f'Invalid offset key {row.subject}/{row.run}')
        method, evidence = str(row.method).strip().lower(), str(row.evidence).strip()
        if method == 'verified' and (not evidence or evidence.lower() in ('none', 'nan')):
            raise ValueError(f'Verified alignment requires evidence: {subject}/{run}')
        offset = float(row.offset_seconds)
        if not np.isfinite(offset):
            raise ValueError(f'Invalid offset: {subject}/{run}')
        result[(subject, run)] = (offset, method, evidence)
    return result


def match_explorer_style(sp, he, tolerance):
    refs = list(sp.itertuples(index=False))
    hems = list(he.itertuples(index=False))
    if not refs or not hems:
        return []
    low = max(refs[0].ref_t, hems[0].he_t_aligned)
    high = min(refs[-1].ref_t, hems[-1].he_t_aligned)
    used, result = set(), []
    for sr in refs:
        if sr.ref_t < low or sr.ref_t > high:
            continue
        best_idx, best_dt = None, float('inf')
        for i, hr in enumerate(hems):
            if i in used or hr.he_t_aligned < low or hr.he_t_aligned > high:
                continue
            dt = abs(hr.he_t_aligned - sr.ref_t)
            if dt < best_dt:
                best_idx, best_dt = i, dt
            if hr.he_t_aligned > sr.ref_t + tolerance and best_idx is not None:
                break
        if best_idx is not None and best_dt <= tolerance:
            used.add(best_idx)
            hr = hems[best_idx]
            result.append({
                'sp20_elapsed_s': sr.ref_t, 'sp20_timestamp': sr.sp20_time,
                'hemox_elapsed_s': hr.hemox_t, 'hemox_aligned_s': hr.he_t_aligned,
                'delta_t_s': hr.he_t_aligned-sr.ref_t,
                'sp20_spo2': sr.ref_spo2, 'hemox_spo2': hr.he_spo2,
                'sp20_hr': sr.ref_hr, 'hemox_hr': hr.he_hr,
            })
    return result


def build_pairs(data_dir, duration, tolerance, offsets):
    participants = find_participants(data_dir)
    paired, audit = [], []
    for subject in SUBJECTS:
        indexed = index_sources(participants[subject], subject)
        for run in SESSIONS:
            h = indexed[run]['hemox'][0]
            s = indexed[run]['sp20'][0]
            he, sp = read_fused(h), read_sp20(s)
            he = he[(he.hemox_t>=0)&(he.hemox_t<=duration)].copy()
            sp = sp[(sp.ref_t>=0)&(sp.ref_t<=duration)].copy()
            valid = he[he.quality_valid].copy()
            offset, method, evidence = offsets.get((subject,run),(0.,'provisional_zero_offset',''))
            valid['he_t_aligned'] = valid.hemox_t + offset
            valid = valid[(valid.he_t_aligned>=0)&(valid.he_t_aligned<=duration)].copy()
            records = match_explorer_style(sp, valid, tolerance)
            for row in records:
                row.update(subject=subject,run=run,alignment_status=method,
                           alignment_offset_seconds=offset,alignment_evidence=evidence)
            paired.extend(records)
            audit.append(dict(subject=subject,run=run,source_rows=len(he),
                              quality_valid_rows=len(valid),sp20_reference_rows=len(sp),
                              available_spo2=int(np.isfinite(valid.he_spo2).sum()),
                              available_hr=int(np.isfinite(valid.he_hr).sum()),
                              matched_time_pairs=len(records),
                              matched_spo2_pairs=sum(np.isfinite(q['sp20_spo2']) and np.isfinite(q['hemox_spo2']) for q in records),
                              matched_hr_pairs=sum(np.isfinite(q['sp20_hr']) and np.isfinite(q['hemox_hr']) for q in records),
                              alignment_status=method,offset_seconds=offset,
                              hemox_source=str(h),sp20_source=str(s)))
            print(f'{subject}/{run}: FUSED={len(he)}, quality-valid={len(valid)}, matched={len(records)}')
    if not paired:
        raise RuntimeError('No matched observations; inspect pairing tolerance and alignment offsets')
    x, qc = pd.DataFrame(paired), pd.DataFrame(audit)
    if len(qc)!=40 or qc.subject.nunique()!=10 or qc.duplicated(['subject','run']).any():
        raise ValueError('Expected 40 unique sessions from 10 participants')
    return x, qc


def metric_values(prediction, reference, threshold):
    a=np.asarray(prediction,dtype=float);b=np.asarray(reference,dtype=float)
    valid=np.isfinite(a)&np.isfinite(b)
    a,b=a[valid],b[valid]
    if not len(a):
        return {k:np.nan for k in ('bias','absolute_bias','sd_error','loa_low','loa_high','loa_width',
                                  'mae','rmse','mape_pct','acc_100_minus_mape_pct','within_threshold_pct')} | {'valid_pairs':0}
    e=a-b; bias=float(e.mean()); sd=float(e.std(ddof=1)) if len(e)>1 else np.nan
    good=np.abs(b)>1e-12
    mape=float(100*np.mean(np.abs(e[good]/b[good]))) if good.any() else np.nan
    return dict(valid_pairs=int(len(e)),bias=bias,absolute_bias=abs(bias),sd_error=sd,
                loa_low=bias-1.96*sd, loa_high=bias+1.96*sd, loa_width=3.92*sd,
                mae=float(np.mean(np.abs(e))), rmse=float(np.sqrt(np.mean(e*e))),
                mape_pct=mape, acc_100_minus_mape_pct=100-mape,
                within_threshold_pct=100*float(np.mean(np.abs(e)<=threshold)))


def arrays_by_subject(paired, variable):
    pcol, rcol = (('hemox_spo2','sp20_spo2') if variable=='SpO2' else ('hemox_hr','sp20_hr'))
    results = {}
    for sub in SUBJECTS:
        t = paired.loc[paired.subject.eq(sub), [pcol, rcol]].apply(pd.to_numeric, errors='coerce').dropna()
        t = t.loc[np.isfinite(t[pcol])&np.isfinite(t[rcol])]
        if t.empty:
            raise ValueError(f'{sub} has zero valid {variable} pairs')
        results[sub] = (t[pcol].to_numpy(float), t[rcol].to_numpy(float))
    return results


def fit_m1(pred, ref):
    return {'offset': float(np.mean(pred-ref))}


def apply_m1(pred, pars):
    return pred - pars['offset']


def fit_m2(pred, ref):
    a, b = np.polyfit(pred, ref, 1)
    return {'a': float(a), 'b': float(b)}


def apply_m2(pred, pars):
    return pars['a']*pred + pars['b']


def fit_m3(pred, ref):
    a, b, c = np.polyfit(pred, ref, 2)
    return {'a2': float(a), 'a1': float(b), 'a0': float(c)}


def apply_m3(pred, pars):
    return pars['a2']*pred*pred + pars['a1']*pred + pars['a0']


MODEL_FUNCS = {
    'M1': (fit_m1, apply_m1, 'Offset-only'),
    'M2': (fit_m2, apply_m2, 'Full linear'),
    'M3': (fit_m3, apply_m3, 'Quadratic'),
}


def exhaustive_study(pairs):
    fold_rows = []
    coef_rows = []
    for variable in ('SpO2', 'HR'):
        group = arrays_by_subject(pairs, variable)
        bound = 3.0 if variable=='SpO2' else 5.0
        for k in CALIBRATION_SIZES:
            fold_count = 0
            for fold_index, selected in enumerate(combinations(SUBJECTS, k), start=1):
                fold_count = fold_index
                calib = tuple(selected)
                validation = tuple(sub for sub in SUBJECTS if sub not in calib)
                assert not set(calib).intersection(validation)
                vp = np.concatenate([group[sub][0] for sub in validation])
                vr = np.concatenate([group[sub][1] for sub in validation])
                before = metric_values(vp, vr, bound)
                base = dict(variable=variable, k_calibration_participants=k,
                            calibration_runs=4*k, validation_participants=10-k,
                            validation_runs=4*(10-k), fold_index=fold_index,
                            calibration_ids=';'.join(calib), validation_ids=';'.join(validation),
                            validation_valid_pairs=len(vp))
                fold_rows.append(dict(**base, model='Before', fit_label='Native uncalibrated',
                                      calibration_valid_pairs=0 if k==0 else np.nan, **before))
                for model_name in MODEL_ORDER:
                    if k == 0:
                        after = before.copy()
                        pars = {}
                        fit_label = f'{MODEL_FUNCS[model_name][2]}; k=0 no fitting, native predictions unchanged'
                        c_pairs = 0
                    else:
                        cp = np.concatenate([group[sub][0] for sub in calib])
                        cr = np.concatenate([group[sub][1] for sub in calib])
                        fit_fn, apply_fn, desc = MODEL_FUNCS[model_name]
                        pars = fit_fn(cp, cr)
                        pred = apply_fn(vp.copy(), pars)
                        if variable == 'SpO2':
                            pred = np.clip(pred, 0, 100)
                        after = metric_values(pred, vr, bound)
                        fit_label = desc
                        c_pairs = len(cp)
                    fold_rows.append(dict(**base, model=model_name, fit_label=fit_label,
                                          calibration_valid_pairs=c_pairs, **after))
                    coef_row = dict(variable=variable, model=model_name, k_calibration_participants=k,
                                    calibration_runs=4*k, validation_participants=10-k,
                                    validation_runs=4*(10-k), fold_index=fold_index,
                                    calibration_ids=';'.join(calib), validation_ids=';'.join(validation))
                    coef_row.update(pars)
                    coef_rows.append(coef_row)
            print(f'{variable}, k={k}: {fold_count} participant-disjoint calibration combinations')
    return pd.DataFrame(fold_rows), pd.DataFrame(coef_rows)


METRICS = [
    ('bias','Bias','pp','bpm'),
    ('absolute_bias','Absolute_Bias','pp','bpm'),
    ('sd_error','SD_error','pp','bpm'),
    ('loa_width','LoA_width','pp','bpm'),
    ('mae','MAE','pp','bpm'),
    ('rmse','RMSE','pp','bpm'),
    ('acc_100_minus_mape_pct','ACC','%','%'),
]


def summarize_folds(folds):
    rows=[]
    for (variable,k,model),group in folds.groupby(['variable','k_calibration_participants','model'], sort=True):
        row=dict(variable=variable, k_calibration_participants=int(k), model=model,
                 n_combinations=int(group.fold_index.nunique()), calibration_runs=4*int(k),
                 validation_participants=10-int(k), validation_runs=4*(10-int(k)),
                 calibration_n_pairs_min=int(np.nanmin(group.calibration_valid_pairs)) if np.isfinite(group.calibration_valid_pairs).any() else 0,
                 calibration_n_pairs_max=int(np.nanmax(group.calibration_valid_pairs)) if np.isfinite(group.calibration_valid_pairs).any() else 0,
                 validation_n_pairs_min=int(group.validation_valid_pairs.min()),
                 validation_n_pairs_max=int(group.validation_valid_pairs.max()))
        for col,_,_,_ in METRICS:
            z = group[col].to_numpy(float)
            row[col+'_mean'] = float(np.nanmean(z))
            row[col+'_sd_across_combinations'] = float(np.nanstd(z, ddof=1)) if len(z)>1 else 0.0
            row[col+'_median'] = float(np.nanmedian(z))
            row[col+'_q25'] = float(np.nanquantile(z, .25))
            row[col+'_q75'] = float(np.nanquantile(z, .75))
            row[col+'_min'] = float(np.nanmin(z))
            row[col+'_max'] = float(np.nanmax(z))
        row['within_threshold_pct_mean']=float(group.within_threshold_pct.mean())
        row['within_threshold_pct_sd_across_combinations']=float(group.within_threshold_pct.std(ddof=1)) if len(group)>1 else 0.0
        row['validation_valid_pairs_mean']=float(group.validation_valid_pairs.mean())
        rows.append(row)
    return pd.DataFrame(rows)


def load_fixed_baseline(path, pairs):
    """Return the fixed full-N=10 native (M0) baseline.

    Preferred behavior:
      1) If ``path`` exists, validate it exactly against the current paired data
         and use it. This preserves reproducibility with a previously frozen
         ``overall_performance.csv``.
      2) If it does not exist, derive the same baseline directly from the
         current native paired measurements instead of aborting. The derived
         baseline is subsequently exported to Tables/fixed_overall_baseline.csv.

    No calibration is applied in either case.
    """
    specs = (
        ('SpO2', 'hemox_spo2', 'sp20_spo2', 3.0),
        ('HR',   'hemox_hr',   'sp20_hr',   5.0),
    )

    # Always compute the baseline from the current native pair population so
    # that an external file, when supplied, can be audited against it.
    derived_by_variable = {
        v: metric_values(pairs[pred], pairs[ref], threshold)
        for v, pred, ref, threshold in specs
    }

    if path.is_file():
        baseline = pd.read_csv(path)
        needed = {
            'variable','participants','runs','valid_pairs','bias','sd_error',
            'loa_width','mae','rmse','acc_100_minus_mape_pct'
        }
        absent = needed - set(baseline.columns)
        if absent:
            raise ValueError(
                f'Invalid overall_performance.csv; missing columns: {sorted(absent)}'
            )
        if len(baseline) != 2 or set(baseline.variable.astype(str)) != {'SpO2','HR'}:
            raise ValueError(
                'Overall baseline must contain exactly one SpO2 row and one HR row.'
            )

        audit = []
        for v, pred, ref, threshold in specs:
            given = baseline.loc[baseline.variable.eq(v)].iloc[0]
            derived = derived_by_variable[v]
            if int(given.valid_pairs) != derived['valid_pairs']:
                raise ValueError(
                    f'{v}: baseline valid_pairs mismatch '
                    f'({int(given.valid_pairs)} vs {derived["valid_pairs"]})'
                )
            for metric in (
                'bias','sd_error','loa_width','mae','rmse',
                'acc_100_minus_mape_pct'
            ):
                if not np.isclose(
                    float(given[metric]), float(derived[metric]),
                    rtol=1e-8, atol=1e-8
                ):
                    raise ValueError(
                        f'{v}: baseline metric mismatch for {metric}. '
                        'Use matching overall_performance.csv.'
                    )
            audit.append(dict(
                variable=v, participants=10, runs=40,
                baseline_source=str(path.resolve()),
                baseline_type='fixed overall N=10, before calibration',
                baseline_absolute_bias=abs(float(given.bias)),
                **{metric: float(given[metric]) for metric in (
                    'bias','sd_error','loa_width','mae','rmse',
                    'acc_100_minus_mape_pct'
                )},
                valid_pairs=int(given.valid_pairs),
            ))
        print(f'[BASELINE] Validated external fixed baseline: {path}')
        return pd.DataFrame(audit)

    # Missing external CSV: derive the fixed M0 baseline from exactly the same
    # native paired measurements used by this run. This is mathematically the
    # k=0 full-N=10 baseline and avoids a needless hard dependency on a CSV.
    print(
        f'[BASELINE] {path.name} not found; deriving fixed full-N=10 native '
        'baseline directly from current paired measurements.'
    )
    audit = []
    for v, pred, ref, threshold in specs:
        d = derived_by_variable[v]
        audit.append(dict(
            variable=v, participants=10, runs=40,
            baseline_source='derived_from_current_native_pairs',
            baseline_type='fixed overall N=10, before calibration (derived)',
            baseline_absolute_bias=float(d['absolute_bias']),
            **{metric: float(d[metric]) for metric in (
                'bias','sd_error','loa_width','mae','rmse',
                'acc_100_minus_mape_pct'
            )},
            valid_pairs=int(d['valid_pairs']),
        ))
    return pd.DataFrame(audit)


def make_chart(summary, fixed_baseline, variable, col, filename, ylabel, unit):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    base = fixed_baseline.loc[fixed_baseline.variable.eq(variable)].iloc[0]
    const = float(base['baseline_absolute_bias'] if col=='absolute_bias' else base[col])
    ax.axhline(const, color=COLORS['BeforeBaseline'], linestyle='--', linewidth=1.6,
               label='Before calibration (overall N=10)')
    markers = {'M1':'o', 'M2':'s', 'M3':'^'}
    for model in MODEL_ORDER:
        data = summary[(summary.variable.eq(variable)) & (summary.model.eq(model))].sort_values('k_calibration_participants')
        x = data.k_calibration_participants.to_numpy(float)
        y = data[col+'_mean'].to_numpy(float)
        ax.plot(x, y, '-', marker=markers[model], color=COLORS[model], linewidth=1.6, markersize=6, label=model)
    ax.set_xlabel('Number of calibration participants')
    ax.set_ylabel(f'{ylabel} ({unit})')
    ax.set_xticks(list(CALIBRATION_SIZES)); ax.set_xlim(-.2, 5.2)
    ax.grid(True, color='#D0D4D9', linewidth=.65, alpha=.7)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    lg = ax.legend(loc='best', facecolor='white', edgecolor='#A9B3BD', framealpha=.96)
    for t in lg.get_texts():
        t.set_fontweight('bold')
    for tick in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick.set_fontweight('bold')
    fig.subplots_adjust(left=.165, right=.98, bottom=.225, top=.97)
    fig.savefig(filename, dpi=DPI, facecolor='white', pil_kwargs={'quality':96,'subsampling':0})
    plt.close(fig)


def plot_joint_rmse(summary, filename):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    markers = {'M1':'o', 'M2':'s', 'M3':'^'}
    linestyles = {'SpO2':'-', 'HR':'--'}
    for model in MODEL_ORDER:
        for variable in ('SpO2', 'HR'):
            data = summary[(summary.variable.eq(variable)) & (summary.model.eq(model))].sort_values('k_calibration_participants')
            label = f'{model} {variable} RMSE'
            ax.plot(data.k_calibration_participants.to_numpy(float), data.rmse_mean.to_numpy(float),
                    linestyle=linestyles[variable], marker=markers[model], color=COLORS[model],
                    linewidth=1.6, markersize=5.8, label=label)
    ax.set_xlabel('Number of calibration participants')
    ax.set_ylabel('Held-out validation RMSE')
    ax.set_xticks(list(CALIBRATION_SIZES)); ax.set_xlim(-.2, 5.2)
    ax.grid(True, color='#D0D4D9', linewidth=.65, alpha=.7)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    leg = ax.legend(loc='best', facecolor='white', edgecolor='#A9B3BD', framealpha=.96, ncol=2)
    for t in leg.get_texts():
        t.set_fontweight('bold')
    for tick in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        tick.set_fontweight('bold')
    fig.subplots_adjust(left=.17, right=.98, bottom=.225, top=.97)
    fig.savefig(filename, dpi=DPI, facecolor='white', pil_kwargs={'quality':96,'subsampling':0})
    plt.close(fig)


def build_model_decision_table(summary):
    rows=[]
    for k in CALIBRATION_SIZES:
        for variable in ('SpO2','HR'):
            block = summary[(summary.variable.eq(variable)) & (summary.k_calibration_participants.eq(k)) & (summary.model.isin(MODEL_ORDER))]
            best_rmse = block.loc[block.rmse_mean.idxmin()]
            best_mae = block.loc[block.mae_mean.idxmin()]
            best_abs_bias = block.loc[block.absolute_bias_mean.idxmin()]
            rows.append(dict(k_calibration_participants=k, variable=variable,
                             best_rmse_model=best_rmse.model, best_rmse=float(best_rmse.rmse_mean),
                             best_mae_model=best_mae.model, best_mae=float(best_mae.mae_mean),
                             best_abs_bias_model=best_abs_bias.model, best_abs_bias=float(best_abs_bias.absolute_bias_mean)))
    return pd.DataFrame(rows)


def write_findings(out, summary, decision, fixed_baseline, qc, args, font):
    lines=[
        '# HEMOX: participant-disjoint multi-model calibration study (M1, M2, M3)', '',
        f'Dataset: N={qc.subject.nunique()}, sessions={len(qc)}, {5*len(qc)} recording minutes; native FUSED_TRENDS values with quality_valid=true. Font: {font}.',
        f'Reference pairing: Explorer-style one-to-one nearest elapsed time, tolerance {args.tolerance:g} s.',
        'Without independently documented timestamps, default zero-offset matching is provisional.',
        'Models: M1 offset-only, M2 full linear, M3 quadratic. SpO2 calibrated outputs clipped to [0,100].',
        'Red horizontal reference in the figures: fixed full-N=10 native uncalibrated baseline; validated from overall_performance.csv when available, otherwise derived from the current paired native measurements.', '',
        '## Fixed overall N=10 baseline', ''
    ]
    for var in ('SpO2','HR'):
        unit='pp' if var=='SpO2' else 'bpm'
        bas=fixed_baseline.loc[fixed_baseline.variable.eq(var)].iloc[0]
        lines.append(f'- {var}: bias={bas.bias:+.3f} {unit}, SD={bas.sd_error:.3f} {unit}, LoA width={bas.loa_width:.3f} {unit}, MAE={bas.mae:.3f} {unit}, RMSE={bas.rmse:.3f} {unit}, ACC={bas.acc_100_minus_mape_pct:.2f}%.')
    lines.extend(['', '## Average held-out validation RMSE by model', ''])
    for var in ('SpO2','HR'):
        unit='pp' if var=='SpO2' else 'bpm'
        lines.append(f'### {var}')
        for k in CALIBRATION_SIZES:
            block = summary[(summary.variable.eq(var)) & (summary.k_calibration_participants.eq(k)) & (summary.model.isin(MODEL_ORDER))].sort_values('model')
            chunks = [f"{r.model}={r.rmse_mean:.3f} {unit}" for r in block.itertuples()]
            best = block.loc[block.rmse_mean.idxmin()]
            lines.append(f'- k={k}: ' + '; '.join(chunks) + f'; best RMSE: {best.model}.')
        lines.append('')
    # Final recommendation heuristic
    sp_best = summary[(summary.variable=='SpO2') & (summary.model.isin(MODEL_ORDER))].groupby('model')['rmse_mean'].mean().sort_values()
    hr_best = summary[(summary.variable=='HR') & (summary.model.isin(MODEL_ORDER))].groupby('model')['rmse_mean'].mean().sort_values()
    avg_norm_rows=[]
    for model in MODEL_ORDER:
        sp = summary[(summary.variable=='SpO2') & (summary.model==model)].set_index('k_calibration_participants').rmse_mean
        hr = summary[(summary.variable=='HR') & (summary.model==model)].set_index('k_calibration_participants').rmse_mean
        avg_norm_rows.append((model, float((sp/sp.min()).mean() + (hr/hr.min()).mean())))
    joint_rank = sorted(avg_norm_rows, key=lambda x: x[1])
    recommended = joint_rank[0][0]
    lines.extend([
        '## Suggested decision logic', '',
        f'- Mean held-out SpO2 RMSE across k is lowest for: {sp_best.index[0]}.',
        f'- Mean held-out HR RMSE across k is lowest for: {hr_best.index[0]}.',
        f'- Simple joint heuristic (sum of RMSE values normalized to each output\'s per-k minimum) selects: {recommended}.',
        '- Prefer M1 if its held-out RMSE/MAE are close to M2/M3 because M1 is simpler, easier to interpret, and less prone to overfitting.',
        '- Prefer M2 only if it yields a consistent validation improvement over M1 for both SpO2 and HR or for the primary endpoint without degrading the other materially.',
        '- Prefer M3 only if its held-out validation performance is reproducibly better than both M1 and M2; otherwise treat any calibration-set fitting advantage as likely overfitting.',
        '', '## Best model by metric and calibration size', ''
    ])
    for row in decision.itertuples(index=False):
        unit='pp' if row.variable=='SpO2' else 'bpm'
        lines.append(f'- {row.variable}, k={row.k_calibration_participants}: best RMSE {row.best_rmse_model} ({row.best_rmse:.3f} {unit}); best MAE {row.best_mae_model} ({row.best_mae:.3f} {unit}); best |Bias| {row.best_abs_bias_model} ({row.best_abs_bias:.3f} {unit}).')
    lines.extend(['', '## Interpretation limits', '',
                  '- The model selection remains exploratory because the same family of held-out folds is used to compare M1, M2 and M3.',
                  '- Participant combinations change with k; compare models at the same k rather than across different k alone.',
                  '- Across-combination SD in the CSV tables is descriptive only and is not a confidence interval.',
                  '- Comparing to a commercial SP20 does not constitute arterial-reference pulse-oximeter accuracy validation.',
                  '- Confirm the selected model and calibration size on an untouched external test cohort if possible.'])
    (out/'Findings_multi_calibration_study.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def main():
    args = parse_args()
    if args.duration<=0 or args.tolerance<=0:
        raise ValueError('Duration and pairing tolerance must be positive')
    root = args.data_dir.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f'Data root not found: {root}')
    font = setup_plot_style()
    out = args.out_dir.resolve(); figs = out/'Figures'; tables = out/'Tables'
    figs.mkdir(parents=True, exist_ok=True); tables.mkdir(parents=True, exist_ok=True)
    offsets = load_offsets(args.offsets_csv)
    pairs, qc = build_pairs(root, args.duration, args.tolerance, offsets)
    fixed_baseline = load_fixed_baseline(args.baseline_csv, pairs)
    fixed_baseline.to_csv(tables/'fixed_overall_baseline.csv', index=False)
    pairs.to_csv(tables/'paired_native_measurements.csv', index=False)
    qc.to_csv(tables/'source_alignment_availability.csv', index=False)
    folds, coeffs = exhaustive_study(pairs)
    folds.to_csv(tables/'all_participant_combinations_models.csv', index=False)
    coeffs.to_csv(tables/'model_coefficients_by_fold.csv', index=False)
    summary = summarize_folds(folds)
    summary.to_csv(tables/'validation_metrics_by_calibration_size_and_model.csv', index=False)
    decision = build_model_decision_table(summary)
    decision.to_csv(tables/'best_model_by_metric_and_calibration_size.csv', index=False)
    figure_rows=[]
    for variable in ('SpO2','HR'):
        bas = fixed_baseline.loc[fixed_baseline.variable.eq(variable)].iloc[0]
        for metric,_,_,_ in METRICS:
            fixed = float(bas.baseline_absolute_bias if metric=='absolute_bias' else bas[metric])
            for model in MODEL_ORDER:
                for k in CALIBRATION_SIZES:
                    rec = summary[(summary.variable.eq(variable)) & (summary.model.eq(model)) & (summary.k_calibration_participants.eq(k))].iloc[0]
                    figure_rows.append(dict(variable=variable, metric=metric, model=model,
                                            k_calibration_participants=k, fixed_overall_N10_baseline=fixed,
                                            heldout_mean=float(rec[metric+'_mean']),
                                            heldout_sd_across_combinations=float(rec[metric+'_sd_across_combinations']),
                                            n_combinations=int(rec.n_combinations)))
    pd.DataFrame(figure_rows).to_csv(tables/'figure_data_fixed_baseline_vs_models.csv', index=False)
    for variable in ('SpO2','HR'):
        for col, name, spunit, hrunit in METRICS:
            unit = spunit if variable=='SpO2' else hrunit
            filename = figs / f'Study_{variable}_{name}.jpg'
            make_chart(summary, fixed_baseline, variable, col, filename, name.replace('_',' '), unit)
            print('Saved figure:', filename.name)
    plot_joint_rmse(summary, figs/'Study_Joint_RMSE_M1_M2_M3.jpg')
    print('Saved figure: Study_Joint_RMSE_M1_M2_M3.jpg')
    write_findings(out, summary, decision, fixed_baseline, qc, args, font)
    expected_splits = sum(len(list(combinations(SUBJECTS, k))) for k in CALIBRATION_SIZES)
    assert expected_splits == 638 and len(summary) == len(CALIBRATION_SIZES)*2*4
    print(f'Completed: participants={qc.subject.nunique()}, runs={len(qc)}, participant combinations={expected_splits}, metric rows={len(folds)}, figures={len(METRICS)*2+1}')
    print('Output:', out)


if __name__=='__main__':
    main()
