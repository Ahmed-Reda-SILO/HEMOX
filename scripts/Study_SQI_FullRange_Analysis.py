
from __future__ import annotations

import argparse
import io
import re
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
import numpy as np
import pandas as pd

SUBJECTS = tuple(f'P{i:02d}' for i in range(3, 11))
RUNS = ('S01', 'S02', 'NS01', 'NS02')
RUN_KEY = re.compile(r'(?:HEMOX|SP20)?P0*(\d{1,2})(NS0[12]|S0[12])', re.I)
M1_OFFSET_PP = 0.906754
DPI = 600
BLUE, ORANGE = '#0072B2', '#D55E00'
DARK, GRID = '#333333', '#CCD3DC'
# Bins are fixed BEFORE looking at agreement: left-closed, right-open except final.
BIN_EDGES = (0.0, .65, .75, .80, .85, .90, .95, 1.0)


def cli():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description='HEMOX absolute error vs fused SQI and participant-bootstrap bin agreement')
    p.add_argument('--data-dir', type=Path, default=here, help='Folder containing P03..P10 folders or ZIPs')
    p.add_argument('--manifest', type=Path, help='Optional participant,archive CSV selecting exactly one source for each P03..P10')
    p.add_argument('--out-dir', type=Path, default=here/'Study_Absolute_Error_vs_Fused_SQI_outputs')
    p.add_argument('--offsets-csv', type=Path, help='Optional participant,run,offset_seconds,method,evidence; HEMOX elapsed + offset. Absent run: provisional 0 offset.')
    p.add_argument('--tolerance', type=float, default=1.0, help='Max SP20/HEMOX elapsed match difference, seconds')
    p.add_argument('--spo2-model', choices=('m1','m0'), default='m1')
    p.add_argument('--spo2-offset', type=float, default=M1_OFFSET_PP)
    p.add_argument('--bootstrap', type=int, default=2000, help='Participant-level resamples (default 2000)')
    p.add_argument('--seed', type=int, default=20260928)
    p.add_argument('--zoom-min', type=float, default=0.75, help='Zoomed companion x-axis lower bound; full-range remains 0--1 (default 0.75)')
    return p.parse_args()


class FolderSource:
    def __init__(self, root):
        self.root = root.resolve()
        self.files = {f.relative_to(self.root).as_posix(): f for f in self.root.rglob('*') if f.is_file()}
    def __enter__(self): return self
    def __exit__(self, *exc): return False
    def namelist(self): return list(self.files)
    def read(self, member): return self.files[member].read_bytes()


def select_sources(root, manifest):
    root = root.expanduser().resolve()
    if not root.is_dir(): raise FileNotFoundError(f'Data directory not found: {root}')
    if manifest is not None:
        mpath = manifest if manifest.is_absolute() else root/manifest
        m = pd.read_csv(mpath, dtype=str).fillna('')
        if not {'participant','archive'}.issubset(m.columns):
            raise ValueError('Manifest columns required: participant,archive')
        subs = m.participant.astype(str).str.upper().str.strip()
        if len(m)!=8 or subs.duplicated().any() or set(subs)!=set(SUBJECTS):
            raise ValueError('Manifest must specify P03..P10 exactly once (P01/P02 excluded)')
        selections = {}
        for row in m.itertuples(index=False):
            sub = str(row.participant).strip().upper()
            src = Path(str(row.archive).strip())
            src = (src if src.is_absolute() else root/src).resolve()
            if not src.is_dir() and not (src.is_file() and src.suffix.lower()=='.zip'):
                raise FileNotFoundError(f'{sub}: missing selected folder/ZIP: {src}')
            selections[sub] = src
        return selections
    selections = {}
    for sub in SUBJECTS:
        folders = [f for f in root.iterdir() if f.is_dir() and f.name.upper()==sub]
        if len(folders)==1:
            selections[sub] = folders[0]
            continue
        if len(folders)>1: raise ValueError(f'Ambiguous folders for {sub}: {folders}')
        zips = [f for f in root.iterdir() if f.is_file() and re.fullmatch(
            re.escape(sub)+r'(?:\(\d+\))?\.zip', f.name, re.I)]
        if len(zips)!=1:
            raise FileNotFoundError(f'{sub}: missing/ambiguous source: {zips}. Use --manifest to choose corrected ZIPs.')
        selections[sub] = zips[0]
    return selections


def list_run_members(source, sub):
    out = {run: {'FUSED': [], 'WINDOWS': [], 'SP20': []} for run in RUNS}
    for member in source.namelist():
        basename = Path(member).name.upper()
        mt = RUN_KEY.search(basename)
        if not mt or f'P{int(mt.group(1)):02d}'!=sub: continue
        run = mt.group(2).upper()
        if basename.startswith('SP20') and basename.endswith('.CSV'):
            out[run]['SP20'].append(member)
        elif basename.endswith('_FUSED_TRENDS.CSV') and 'CALIBRATED' not in basename:
            out[run]['FUSED'].append(member)
        elif basename.endswith('_WINDOWS.CSV'):
            out[run]['WINDOWS'].append(member)
    for run in RUNS:
        for typ, paths in out[run].items():
            if len(paths)!=1:
                raise FileNotFoundError(f'{sub}/{run}: exactly one {typ} required, found {len(paths)}: {paths}')
    return out


def parse_bool(series, sub, run):
    vals = series.astype(str).str.strip().str.lower()
    true = {'1','true','yes','accepted','valid'}
    false = {'0','false','no','invalid','rejected'}
    unknown = sorted(set(vals)-true-false)
    if unknown: raise ValueError(f'{sub}/{run}: unexpected quality_valid tokens: {unknown}')
    return vals.isin(true)


def read_hemox(source, path, sub, run):
    h = pd.read_csv(io.BytesIO(source.read(path)),low_memory=False)
    needed = {'elapsed_s','quality_valid','fused_sqi','fused_spo2_percent','fused_hr_bpm'}
    if not needed.issubset(h): raise ValueError(f'{sub}/{run}: missing FUSED fields: {needed-set(h)}')
    t = pd.to_numeric(h.elapsed_s,errors='coerce')
    if t.isna().any() or t.duplicated().any() or not t.is_monotonic_increasing:
        raise ValueError(f'{sub}/{run}: FUSED elapsed time missing, duplicated or nonmonotonic')
    return pd.DataFrame(dict(hemox_t=t.astype(float),quality_valid=parse_bool(h.quality_valid,sub,run),
        fused_sqi=pd.to_numeric(h.fused_sqi,errors='coerce'),
        native_spo2=pd.to_numeric(h.fused_spo2_percent,errors='coerce'),
        native_hr=pd.to_numeric(h.fused_hr_bpm,errors='coerce')))


def read_windows_channel_scores(source, path, sub, run):
    """Preserve the original 2-s scheduled WINDOW records and actual CH status.

    Channel SQI is not fused SQI; missing entries stay missing and are never
    assigned an invented fused score. Accepted status is not equivalent to
    SpO2 eligibility, so we retain that separate exported flag as well.
    """
    w=pd.read_csv(io.BytesIO(source.read(path)),low_memory=False)
    shared={'window','start_s','end_s'}
    needed=shared | {f'ch{ch}_{key}' for ch in (1,2) for key in
        ('sqi','status','spo2_eligible','sqi_accept_threshold','sqi_reject_threshold')}
    if not needed.issubset(w):
        raise ValueError(f'{sub}/{run}: missing WINDOWS columns {needed-set(w)}')
    window=pd.to_numeric(w.window,errors='coerce')
    end=pd.to_numeric(w.end_s,errors='coerce')
    start=pd.to_numeric(w.start_s,errors='coerce')
    if (window.isna().any() or window.duplicated().any() or
        end.isna().any() or end.duplicated().any() or not end.is_monotonic_increasing or
        start.isna().any() or (end<=start).any()):
        raise ValueError(f'{sub}/{run}: invalid window counters or timestamps')
    rows=[]
    for ch in (1,2):
        val=pd.to_numeric(w[f'ch{ch}_sqi'],errors='coerce')
        statuses=w[f'ch{ch}_status'].astype(str).str.strip().str.lower()
        threshold_a=pd.to_numeric(w[f'ch{ch}_sqi_accept_threshold'],errors='coerce')
        threshold_r=pd.to_numeric(w[f'ch{ch}_sqi_reject_threshold'],errors='coerce')
        eligible=w[f'ch{ch}_spo2_eligible'].astype(str).str.strip().str.lower().isin(('true','1','yes'))
        rows.append(pd.DataFrame(dict(participant=sub,run=run,
            condition='Stationary' if run.startswith('S') else 'Controlled Motion',
            window=window.to_numpy(),start_s=start.to_numpy(),end_s=end.to_numpy(),
            channel=f'CH{ch}',channel_sqi=val.to_numpy(),
            status=statuses.to_numpy(),status_accepted=statuses.eq('accepted').to_numpy(),
            spo2_eligible=eligible.to_numpy(),
            sqi_accept_threshold=threshold_a.to_numpy(),
            sqi_reject_threshold=threshold_r.to_numpy())))
    return pd.concat(rows,ignore_index=True)


def source_fused_score_rows(h, sub, run):
    """All native FUSED rows, not just valid/pairable rows; preserve NaNs."""
    o=h.copy()
    o['participant']=sub; o['run']=run
    o['condition']='Stationary' if run.startswith('S') else 'Controlled Motion'
    o['source_fused_sqi_available']=np.isfinite(o.fused_sqi)&o.fused_sqi.between(0,1,inclusive='both')
    o['source_spo2_valid']=(o.quality_valid & o.source_fused_sqi_available &
        np.isfinite(o.native_spo2)&o.native_spo2.between(0,100,inclusive='both'))
    o['source_hr_valid']=(o.quality_valid & o.source_fused_sqi_available &
        np.isfinite(o.native_hr)&o.native_hr.gt(0))
    return o


def read_sp20(source,path,sub,run):
    # First SP20 CSV line contains metadata; second is '# Time, SpO2(%), PR(bpm)'.
    r = pd.read_csv(io.BytesIO(source.read(path)),skiprows=1,skipinitialspace=True,low_memory=False)
    r.columns = [str(c).strip().lstrip('#').strip() for c in r.columns]
    needed = {'Time','SpO2(%)','PR(bpm)'}
    if not needed.issubset(r): raise ValueError(f'{sub}/{run}: missing SP20 columns: {needed-set(r)}')
    absolute = pd.to_datetime(r.Time.astype(str).str.strip(),format='%m/%d/%Y %H:%M:%S',errors='coerce')
    if absolute.isna().any() or absolute.duplicated().any() or not absolute.is_monotonic_increasing:
        raise ValueError(f'{sub}/{run}: invalid/duplicated SP20 timestamp')
    r = pd.DataFrame(dict(ref_t=(absolute-absolute.iloc[0]).dt.total_seconds().astype(float),
                          ref_timestamp=absolute.astype(str),
                          ref_spo2=pd.to_numeric(r['SpO2(%)'],errors='coerce'),
                          ref_hr=pd.to_numeric(r['PR(bpm)'],errors='coerce')))
    r.loc[~r.ref_spo2.between(0,100,inclusive='both'),'ref_spo2']=np.nan
    r.loc[~(np.isfinite(r.ref_hr)&(r.ref_hr>0)),'ref_hr']=np.nan
    return r


def load_offsets(filename):
    if filename is None: return {}
    d = pd.read_csv(filename,dtype=str).fillna('')
    needed = {'participant','run','offset_seconds','method','evidence'}
    if not needed.issubset(d): raise ValueError(f'Offsets CSV requires: {sorted(needed)}')
    out = {}
    for r in d.to_dict('records'):
        sub,run = r['participant'].strip().upper(),r['run'].strip().upper()
        key=(sub,run)
        if sub not in SUBJECTS or run not in RUNS or key in out:
            raise ValueError(f'Unexpected/duplicate offsets key: {key}')
        offset=float(r['offset_seconds'])
        if not np.isfinite(offset):raise ValueError(f'{key}: non-finite offset')
        method,evidence=r['method'].strip().lower(),r['evidence'].strip()
        out[key]=(offset,method=='verified' and len(evidence)>=5,method,evidence)
    return out


def pair_reference_first(reference, hemox, tolerance):
    """SP20-first nearest unused quality-valid HEMOX, overlap limited, ±tolerance.
    Match before per-variable exclusion so SpO2/HR use the same time pairing policy.
    """
    rrows=list(reference.itertuples(index=False)); hrows=list(hemox.itertuples(index=False))
    if not rrows or not hrows:return []
    lo=max(rrows[0].ref_t,hrows[0].hemox_aligned_t)
    hi=min(rrows[-1].ref_t,hrows[-1].hemox_aligned_t)
    used=set(); matches=[]
    for r in rrows:
        if not lo<=r.ref_t<=hi:continue
        best=None; delta=float('inf')
        for i,h in enumerate(hrows):
            if i in used or not lo<=h.hemox_aligned_t<=hi:continue
            d=abs(float(h.hemox_aligned_t)-float(r.ref_t))
            if d<delta:best=i;delta=d
            if h.hemox_aligned_t>r.ref_t+tolerance and best is not None:break
        if best is not None and delta<=tolerance:
            used.add(best);h=hrows[best]
            matches.append(dict(ref_t=float(r.ref_t),ref_timestamp=r.ref_timestamp,
                hemox_t=float(h.hemox_t),hemox_aligned_t=float(h.hemox_aligned_t),
                delta_t_s=float(h.hemox_aligned_t-r.ref_t),fused_sqi=float(h.fused_sqi),
                native_spo2=float(h.native_spo2),native_hr=float(h.native_hr),
                ref_spo2=float(r.ref_spo2),ref_hr=float(r.ref_hr)))
    return matches



def bin_labels():
    labels=[]
    for i,(lo,hi) in enumerate(zip(BIN_EDGES[:-1],BIN_EDGES[1:])):
        right = ']' if i == len(BIN_EDGES)-2 else ')'
        labels.append(f'[{lo:.2f}, {hi:.2f}{right}')
    return labels


def bin_short_labels():
    return [f'{lo:.2f}–{hi:.2f}' for lo,hi in zip(BIN_EDGES[:-1],BIN_EDGES[1:])]


def add_bin_column(df):
    df=df.copy()
    df['sqi_bin']=pd.cut(df.fused_sqi,bins=BIN_EDGES,labels=bin_labels(),
                          right=False,include_lowest=True).astype('object')
    df.loc[df.fused_sqi.eq(1.0),'sqi_bin']=bin_labels()[-1]
    assert df.sqi_bin.notna().all(), 'A valid SQI was not assigned to an interval'
    return df


def bootstrap_cluster_ci(group, subjects, reps, rng):
    """Resample eight WHOLE participant clusters, not individual observations.

    Includes participants with zero observations in the bin in the bootstrap
    population; suppress CI when <3 participants contribute to that bin.
    """
    if int(group.participant.nunique())<3:
        return np.nan,np.nan,0
    s=group.groupby('participant').agg(n=('success','size'),ok=('success','sum'))
    ns=np.array([int(s.loc[p,'n']) if p in s.index else 0 for p in subjects])
    ok=np.array([int(s.loc[p,'ok']) if p in s.index else 0 for p in subjects])
    idx=rng.integers(0,len(subjects),size=(reps,len(subjects)))
    den=ns[idx].sum(axis=1);num=ok[idx].sum(axis=1)
    valid=den>0
    if valid.sum()<max(200,int(.8*reps)):
        return np.nan,np.nan,int(valid.sum())
    ci=np.percentile(100*num[valid]/den[valid],[2.5,97.5])
    return float(ci[0]),float(ci[1]),int(valid.sum())


def selected_measurements(pair, metric):
    if metric=='SpO2':
        val=pair.selected_spo2
        ref=pair.ref_spo2
        valid=np.isfinite(val)&val.between(0,100,inclusive='both')&np.isfinite(ref)
    else:
        val=pair.selected_hr
        ref=pair.ref_hr
        valid=np.isfinite(val)&(val>0)&np.isfinite(ref)
    base=pair.loc[valid].copy()
    validq=np.isfinite(base.fused_sqi)&base.fused_sqi.between(0,1,inclusive='both')
    df=base.loc[validq].copy()
    df['metric']=metric
    df['estimate']=val.loc[df.index].to_numpy(float)
    df['reference']=ref.loc[df.index].to_numpy(float)
    df['error']=df.estimate-df.reference
    df['absolute_error']=df.error.abs()
    df=add_bin_column(df)
    exclusion={'metric':metric,'paired_valid_reference':len(base),
        'excluded_missing_or_invalid_sqi':int((~validq).sum()),
        'analyzed_valid_pairs':len(df),
        'sqi_min':float(df.fused_sqi.min()) if len(df) else np.nan,
        'sqi_max':float(df.fused_sqi.max()) if len(df) else np.nan,
        'n_below_0_75':int((df.fused_sqi<.75).sum())}
    return df,exclusion


def calculate_bins(evaluated,reps,seed):
    rows=[]
    rng=np.random.default_rng(seed)
    for typ in ('SpO2','HR'):
        df=evaluated[typ]
        thresholds=(3.0,) if typ=='SpO2' else (5.0,10.0)
        for limit in thresholds:
            for label,short in zip(bin_labels(),bin_short_labels()):
                group=df.loc[df.sqi_bin==label].copy()
                group['success']=group.absolute_error.le(limit)
                n=len(group);passed=int(group.success.sum());P=int(group.participant.nunique())
                low,high,nb=bootstrap_cluster_ci(group,SUBJECTS,reps,rng)
                rows.append({'metric':typ,'threshold':limit,
                    'unit':'pp' if typ=='SpO2' else 'bpm',
                    'sqi_interval':label,'sqi_short':short,
                    'n_pairs':n,'n_agree':passed,'n_disagree':n-passed,
                    'n_participants':P,'agreement_pct':100*passed/n if n else np.nan,
                    'ci_low_pct':low,'ci_high_pct':high,
                    'bootstrap_nonempty_resamples':nb,
                    'uncertainty':'95% percentile participant-cluster bootstrap',
                    'interpretation':'EMPTY' if n==0 else ('CI unavailable: <3 participant clusters' if P<3 else 'descriptive')})
    result=pd.DataFrame(rows)
    # In each metric/boundary, the bin counts must exactly partition the same data.
    for typ in ('SpO2','HR'):
        for limit in ((3.,) if typ=='SpO2' else (5.,10.)):
            assert result.loc[(result.metric==typ)&(result.threshold==limit),'n_pairs'].sum()==len(evaluated[typ])
    return result


def style():
    installed={f.name for f in font_manager.fontManager.ttflist}
    font=next((f for f in ('Arial Narrow','Arial','Liberation Sans','DejaVu Sans') if f in installed),'DejaVu Sans')
    plt.rcParams.update({'font.family':font,'font.weight':'bold','axes.labelweight':'bold',
        'axes.linewidth':1.25,'axes.labelsize':11.3,'xtick.labelsize':9,'ytick.labelsize':9,
        'savefig.dpi':DPI,'figure.dpi':120,'axes.spines.top':False,'axes.spines.right':False})
    return font


def save(fig,path):
    fig.savefig(path,facecolor='white',dpi=DPI,pil_kwargs={'quality':96,'subsampling':0})
    plt.close(fig)


def scatter_axis(ax,df,typ,xmin,letter=''):
    for cond,c in (('Stationary',BLUE),('Controlled Motion',ORANGE)):
        sub=df[df.condition==cond]
        ax.scatter(sub.fused_sqi,sub.absolute_error,s=10,color=c,
                   alpha=.34,linewidths=0,zorder=3,rasterized=True,
                   label=f'{cond} (n={len(sub)})')
    limits=[(3.,'3 pp','--')] if typ=='SpO2' else [(5.,'5 bpm','--'),(10.,'10 bpm','--')]
    for bound,text,ls in limits:
        ax.axhline(bound,color=DARK if bound!=10 else '#777777',lw=1.25,ls=ls,zorder=2)
        ax.annotate(text+' limit',(1.0,bound),xytext=(-5,4),textcoords='offset points',
                    ha='right',va='bottom',fontsize=8.2,color=DARK)
    maximum=float(df.absolute_error.max()) if len(df) else 0.
    ax.set_ylim(0,max(maximum*1.07,limits[-1][0]*1.09))
    ax.set_xlim(xmin,1.0)
    ax.set_xticks(np.linspace(xmin,1,6))
    ax.set_xlabel('Fused SQI')
    ax.set_ylabel(r'$|\mathrm{HEMOX}-\mathrm{SP20}|$ (pp)' if typ=='SpO2' else
                  r'$|\mathrm{HEMOX}-\mathrm{SP20}|$ (bpm)')
    ax.grid(True,color=GRID,lw=.55,alpha=.7)
    ax.set_axisbelow(True)
    if letter:
        ax.text(.015,.98,f'({letter})',transform=ax.transAxes,ha='left',va='top',
                fontsize=11,weight='bold')


def make_scatter_figures(data,figs,zoom_min):
    made=[]
    legend=[Line2D([0],[0],marker='o',ls='',markerfacecolor=c,markeredgecolor='none',ms=5,label=n)
            for n,c in (('Stationary',BLUE),('Controlled motion',ORANGE))]
    for label,xmin in [('Full_0to1',0.0),('Zoom',zoom_min)]:
        fig,axs=plt.subplots(1,2,figsize=(7.2,3.8))
        for ax,typ,letter in zip(axs,('SpO2','HR'),'ab'):
            scatter_axis(ax,data[typ],typ,xmin,letter)
        fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.995),ncol=2,
                   framealpha=.97,facecolor='white',edgecolor='#A9B3BD',fontsize=9)
        fig.subplots_adjust(left=.097,right=.985,bottom=.18,top=.825,wspace=.40)
        path=figs/f'Fig1_Absolute_Error_{label}_2Panel.jpg';save(fig,path);made.append(path)
        for typ in ('SpO2','HR'):
            fig,ax=plt.subplots(figsize=(7.2,3.8))
            scatter_axis(ax,data[typ],typ,xmin)
            fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.995),ncol=2,
                       framealpha=.97,facecolor='white',edgecolor='#A9B3BD',fontsize=9)
            fig.subplots_adjust(left=.13,right=.985,bottom=.20,top=.82)
            path=figs/f'Fig1_Absolute_Error_{label}_{typ}.jpg';save(fig,path);made.append(path)
    return made


def agreement_axis(ax,stats,typ,letter=''):
    relevant=stats[stats.metric==typ]
    x=np.arange(len(bin_labels()))
    choices=[(3.,BLUE,'SpO₂: within 3 pp',0.)] if typ=='SpO2' else [
        (5.,BLUE,'HR: within 5 bpm',-.075),(10.,'#009E73','HR: within 10 bpm',.075)]
    for limit,color,name,offset in choices:
        rows=relevant[relevant.threshold.eq(limit)].reset_index(drop=True)
        vals=rows.agreement_pct.to_numpy(float)
        xx=x+offset
        ax.plot(xx,vals,'-o',markersize=4.5,color=color,lw=1.8,label=name,zorder=3)
        for i,r in enumerate(rows.itertuples(index=False)):
            if np.isfinite(r.ci_low_pct) and np.isfinite(r.ci_high_pct):
                ax.errorbar([xx[i]],[r.agreement_pct],
                    yerr=[[max(0,r.agreement_pct-r.ci_low_pct)],
                          [max(0,r.ci_high_pct-r.agreement_pct)]],
                    fmt='none',ecolor=color,lw=1.1,capsize=2,zorder=4)
    for i,r in enumerate(relevant[relevant.threshold.eq(choices[0][0])].itertuples(index=False)):
        if not r.n_pairs:
            ax.text(i,8,'No pairs',rotation=90,ha='center',va='bottom',fontsize=7.8,color=DARK)
    base=relevant[relevant.threshold.eq(choices[0][0])]
    ax.set_xticks(x,[f'{r.sqi_short}\nn={r.n_pairs}; P={r.n_participants}' for r in base.itertuples(index=False)],
                  rotation=35,ha='right',rotation_mode='anchor')
    ax.tick_params(axis='x',labelsize=7.4)
    ax.set_ylim(0,106);ax.set_yticks(np.arange(0,101,20));ax.set_xlim(-.45,len(x)-.55)
    ax.set_ylabel('Agreement rate (%)');ax.set_xlabel('Fused SQI interval')
    ax.grid(axis='y',color=GRID,alpha=.66,lw=.55);ax.set_axisbelow(True)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.025),ncol=2 if typ=='HR' else 1,
              fontsize=8.1,framealpha=.97,edgecolor='#A9B3BD')
    if letter:
        ax.text(.018,.95,f'({letter})',transform=ax.transAxes,ha='left',va='top',
                fontsize=11,weight='bold')


def make_agreement_figures(stats,figs):
    made=[]
    fig,axs=plt.subplots(1,2,figsize=(11.4,5.0))
    for ax,typ,letter in zip(axs,('SpO2','HR'),'ab'):
        agreement_axis(ax,stats,typ,letter)
    fig.text(.5,.012,'n = valid paired observations; P = contributing participants; error bars = 95% participant-cluster bootstrap CI',
             ha='center',fontsize=8.5,weight='bold')
    fig.subplots_adjust(left=.07,right=.99,bottom=.36,top=.80,wspace=.28)
    path=figs/'Fig2_Agreement_Rate_by_SQI_2Panel.jpg';save(fig,path);made.append(path)
    for typ in ('SpO2','HR'):
        fig,ax=plt.subplots(figsize=(7.2,4.5))
        agreement_axis(ax,stats,typ)
        fig.text(.5,.013,'n = valid pairs; P = participants; error bars = 95% participant-cluster bootstrap CI',
                 ha='center',fontsize=8.2)
        fig.subplots_adjust(left=.10,right=.985,bottom=.38,top=.80)
        path=figs/f'Fig2_Agreement_Rate_by_SQI_{typ}.jpg';save(fig,path);made.append(path)
    return made


def distribution_table(evaluated):
    rows=[]
    for typ in ('SpO2','HR'):
        df=evaluated[typ]
        for label,short in zip(bin_labels(),bin_short_labels()):
            g=df[df.sqi_bin==label];n=len(g)
            rows.append(dict(metric=typ,sqi_interval=label,sqi_short=short,n_pairs=n,
                n_participants=int(g.participant.nunique()),
                percent_of_metric_pairs=100*n/len(df) if len(df) else np.nan))
    return pd.DataFrame(rows)


def distribution_axis(ax,distribution,typ,letter=''):
    rows=distribution[distribution.metric==typ]
    x=np.arange(len(rows));ns=rows.n_pairs.to_numpy(int)
    bars=ax.bar(x,ns,width=.69,color=BLUE if typ=='SpO2' else ORANGE,alpha=.88,zorder=3)
    for i,(bar,n) in enumerate(zip(bars,ns)):
        ax.text(i,bar.get_height()+max(2,max(ns)*.012),f'{n}',ha='center',va='bottom',
                fontsize=8.3,color=DARK)
    ax.set_ylim(0,max(1,int(max(ns)*1.16)))
    ax.set_ylabel('Valid paired observations (n)');ax.set_xlabel('Fused SQI interval')
    ax.set_xticks(x,rows.sqi_short,rotation=35,ha='right')
    ax.tick_params(axis='x',labelsize=7.7)
    ax.grid(axis='y',color=GRID,lw=.6,alpha=.7);ax.set_axisbelow(True)
    if letter:
        ax.text(.015,.96,f'({letter})',transform=ax.transAxes,ha='left',va='top',
                fontsize=11,weight='bold')


def make_distribution_figures(distribution,figs):
    made=[]
    fig,axs=plt.subplots(1,2,figsize=(10.2,3.9))
    for ax,typ,l in zip(axs,('SpO2','HR'),'ab'):
        distribution_axis(ax,distribution,typ,l)
    fig.subplots_adjust(left=.09,right=.985,bottom=.30,top=.95,wspace=.28)
    path=figs/'Fig3_Observed_SQI_Distribution_2Panel.jpg';save(fig,path);made.append(path)
    for typ in ('SpO2','HR'):
        fig,ax=plt.subplots(figsize=(7.2,3.8))
        distribution_axis(ax,distribution,typ)
        fig.subplots_adjust(left=.13,right=.985,bottom=.31,top=.95)
        path=figs/f'Fig3_Observed_SQI_Distribution_{typ}.jpg';save(fig,path);made.append(path)
    return made


def make_eligibility_figures(window_scores, fused_source, evaluated, figs):
    """Four supplementary JPGs; never substitute channel scores for fused SQI."""
    made=[]
    channel_scored=window_scores.loc[np.isfinite(window_scores.channel_sqi) &
        window_scores.channel_sqi.between(0,1,inclusive='both')].copy()
    channel_scored=channel_scored.rename(columns={'channel_sqi':'fused_sqi'})
    channel_scored=add_bin_column(channel_scored)
    # Status accepted and modality-specific SpO2 eligibility are separate flags.

    def channel_axis(ax,channel,letter=''):
        allrows=window_scores.loc[window_scores.channel.eq(channel)]
        scored=channel_scored.loc[channel_scored.channel.eq(channel)]
        labels=bin_labels(); x=np.arange(len(labels))
        ac=np.array([int(((scored.sqi_bin==label)&scored.status_accepted).sum()) for label in labels])
        no=np.array([int(((scored.sqi_bin==label)&~scored.status_accepted).sum()) for label in labels])
        ax.bar(x,ac,width=.72,color=BLUE,alpha=.90,label='Accepted status',zorder=3)
        ax.bar(x,no,bottom=ac,width=.72,color=ORANGE,alpha=.88,label='Other status',zorder=3)
        ymax=int(max(ac+no)) if len(ac) else 0
        for i,n in enumerate(ac+no):
            if n:
                ax.text(i,n+max(1,ymax*.014),str(n),ha='center',va='bottom',fontsize=8)
        ax.set_ylim(0,max(2,ymax*1.17));ax.set_xlim(-.55,len(x)-.45)
        ax.set_xticks(x,bin_short_labels(),rotation=35,ha='right');ax.tick_params(axis='x',labelsize=8)
        ax.set_ylabel('Scheduled windows with measured channel SQI (n)')
        ax.set_xlabel(f'{channel} SQI interval')
        ax.grid(axis='y',color=GRID,alpha=.65,lw=.55);ax.set_axisbelow(True)
        if letter:ax.text(.015,.97,f'({letter})',transform=ax.transAxes,ha='left',va='top',fontsize=11,weight='bold')
        return len(allrows), len(scored), int(allrows.channel_sqi.isna().sum()), int((scored.fused_sqi<.75).sum()),int(((scored.fused_sqi<.75)&~scored.status_accepted).sum())

    fig,axs=plt.subplots(1,2,figsize=(11.4,4.6))
    details=[channel_axis(ax,ch,l) for ax,ch,l in zip(axs,('CH1','CH2'),'ab')]
    handles=[Line2D([0],[0],color=c,lw=7,label=l) for c,l in ((BLUE,'Accepted status'),(ORANGE,'Other status'))]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.99),ncol=2,framealpha=.97,edgecolor='#A9B3BD',fontsize=10)
    fig.subplots_adjust(left=.075,right=.992,bottom=.29,top=.82,wspace=.33)
    path=figs/'Fig4_All_Window_Channel_SQI_Acceptance_2Panel.jpg';save(fig,path);made.append(path)
    for ch in ('CH1','CH2'):
        fig,ax=plt.subplots(figsize=(7.2,4.2));n,sc,missing,low,low_rejected=channel_axis(ax,ch)
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.995),ncol=2,framealpha=.97,edgecolor='#A9B3BD',fontsize=9)
        fig.subplots_adjust(left=.15,right=.985,bottom=.31,top=.78)
        path=figs/f'Fig4_All_Window_Channel_SQI_Acceptance_{ch}.jpg';save(fig,path);made.append(path)

    scored=fused_source.loc[fused_source.source_fused_sqi_available].copy()
    scored=add_bin_column(scored)
    colors=('#B6BCC4',ORANGE,'#009E73')
    series_labels=('All scored source rows','Valid native output','Valid reference pairs')
    def fused_axis(ax,typ,letter=''):
        labels=bin_labels();x=np.arange(len(labels))
        valid_col='source_spo2_valid' if typ=='SpO2' else 'source_hr_valid'
        alln=np.array([(scored.sqi_bin==b).sum() for b in labels],dtype=int)
        valn=np.array([((scored.sqi_bin==b)&scored[valid_col]).sum() for b in labels],dtype=int)
        paired=evaluated[typ]
        pairn=np.array([(paired.sqi_bin==b).sum() for b in labels],dtype=int)
        # All three populations are nested, but differ in denominators. Never
        # stack the categories (which would triple-count the same observations).
        ax.bar(x,alln,color=colors[0],alpha=.75,width=.78,label=series_labels[0],zorder=2)
        ax.plot(x,valn,'-s',color=colors[1],lw=1.85,ms=5,label=series_labels[1],zorder=4)
        ax.plot(x,pairn,'-o',color=colors[2],lw=1.8,ms=4.5,label=series_labels[2],zorder=5)
        ymax=int(max(alln)) if len(alln) else 0
        for i,n in enumerate(alln):
            if n:ax.text(i,n+max(1,ymax*.015),str(n),fontsize=8,ha='center',va='bottom')
        ax.set_ylim(0,max(2,ymax*1.17));ax.set_xlim(-.55,len(x)-.45)
        ax.set_xticks(x,bin_short_labels(),rotation=35,ha='right');ax.tick_params(axis='x',labelsize=8)
        ax.set_xlabel('Native fused SQI interval');ax.set_ylabel('Source rows / reference pairs (n)')
        ax.grid(axis='y',color=GRID,alpha=.65,lw=.55);ax.set_axisbelow(True)
        if letter:ax.text(.015,.97,f'({letter}) {typ}',transform=ax.transAxes,ha='left',va='top',fontsize=11,weight='bold')
        return alln,valn,pairn
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.6))
    for ax,typ,l in zip(axs,('SpO2','HR'),'ab'):fused_axis(ax,typ,l)
    legend=[Line2D([0],[0],color=c,lw=4 if i==0 else 2,marker=None if i==0 else 'o',label=l)
            for i,(c,l) in enumerate(zip(colors,series_labels))]
    fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.99),ncol=3,framealpha=.97,edgecolor='#A9B3BD',fontsize=9.1)
    fig.subplots_adjust(left=.085,right=.991,bottom=.29,top=.82,wspace=.33)
    path=figs/'Fig5_Native_Fused_SQI_Selection_Cascade_2Panel.jpg';save(fig,path);made.append(path)
    for typ in ('SpO2','HR'):
        fig,ax=plt.subplots(figsize=(7.2,4.2));fused_axis(ax,typ)
        fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.99),ncol=3,framealpha=.97,edgecolor='#A9B3BD',fontsize=8.2)
        fig.subplots_adjust(left=.15,right=.985,bottom=.30,top=.79)
        path=figs/f'Fig5_Native_Fused_SQI_Selection_Cascade_{typ}.jpg';save(fig,path);made.append(path)
    return made


def eligibility_tables(window_scores,fused_source,evaluated):
    """All windows are the denominator; exclude missing scores only from histograms."""
    channel=window_scores.copy()
    eligible=np.isfinite(channel.channel_sqi)&channel.channel_sqi.between(0,1,inclusive='both')
    scored=channel.loc[eligible].rename(columns={'channel_sqi':'fused_sqi'})
    scored=add_bin_column(scored)
    rows=[]
    for ch in ('CH1','CH2'):
        orig=channel[channel.channel.eq(ch)]
        d=scored[scored.channel.eq(ch)]
        for label in bin_labels():
            g=d[d.sqi_bin==label]
            rows.append(dict(channel=ch,sqi_interval=label,scheduled_windows=len(orig),
                scored_windows=len(d),missing_sqi=int(orig.channel_sqi.isna().sum()),
                out_of_range_sqi=int((orig.channel_sqi.notna() & ~orig.channel_sqi.between(0,1,inclusive='both')).sum()),
                n_in_bin=len(g),n_status_accepted=int(g.status_accepted.sum()),
                n_other_status=int((~g.status_accepted).sum()),
                n_spo2_eligible=int(g.spo2_eligible.sum()),
                n_participants=int(g.participant.nunique())))
    by_channel=pd.DataFrame(rows)
    scored_fused=fused_source.loc[fused_source.source_fused_sqi_available].copy()
    scored_fused=add_bin_column(scored_fused)
    rows=[]
    for typ in ('SpO2','HR'):
        sel='source_spo2_valid' if typ=='SpO2' else 'source_hr_valid'
        df=evaluated[typ]
        for label in bin_labels():
            g=scored_fused[scored_fused.sqi_bin==label]
            pair=df[df.sqi_bin==label]
            rows.append(dict(metric=typ,sqi_interval=label,
                native_source_fused_rows=len(fused_source),
                source_fused_sqi_missing=int(fused_source.fused_sqi.isna().sum()),
                scored_native_fused_rows=len(scored_fused),
                n_source_scored_in_bin=len(g),n_valid_native_output_in_bin=int(g[sel].sum()),
                n_valid_reference_pairs_in_bin=len(pair),n_pair_participants=int(pair.participant.nunique())))
    by_fused=pd.DataFrame(rows)
    summary=[]
    for ch in ('CH1','CH2'):
        orig=channel[channel.channel.eq(ch)]
        g=scored[scored.channel.eq(ch)]
        low=g[g.fused_sqi.lt(.75)]
        summary.append(dict(source=ch,scheduled_or_source_rows=len(orig),
            finite_sqi_rows=len(g),missing_sqi_rows=int(orig.channel_sqi.isna().sum()),
            low_sqi_below_0_75=int(len(low)),
            low_sqi_accepted=int(low.status_accepted.sum()),
            low_sqi_other_status=int((~low.status_accepted).sum()),
            valid_native_spo2=np.nan,valid_native_hr=np.nan,
            paired_spo2=np.nan,paired_hr=np.nan))
    g=fused_source[fused_source.source_fused_sqi_available]
    summary.append(dict(source='FUSED_TRENDS',scheduled_or_source_rows=len(fused_source),
        finite_sqi_rows=len(g),missing_sqi_rows=int(fused_source.fused_sqi.isna().sum()),
        low_sqi_below_0_75=int(g.fused_sqi.lt(.75).sum()),low_sqi_accepted=np.nan,
        low_sqi_other_status=np.nan,valid_native_spo2=int(g.source_spo2_valid.sum()),
        valid_native_hr=int(g.source_hr_valid.sum()),paired_spo2=len(evaluated['SpO2']),
        paired_hr=len(evaluated['HR'])))
    return by_channel,by_fused,pd.DataFrame(summary)


def main():
    args=cli()
    if not np.isfinite(args.spo2_offset) or not np.isfinite(args.tolerance) or args.tolerance<=0:
        raise ValueError('SpO2 offset must be finite and time tolerance positive')
    if not (0<=args.zoom_min<1):raise ValueError('--zoom-min must lie in [0,1)')
    if args.bootstrap<300:raise ValueError('--bootstrap must be >=300')
    sources=select_sources(args.data_dir,args.manifest)
    offsets=load_offsets(args.offsets_csv)
    out=args.out_dir.expanduser().resolve();figs=out/'Figures';tables=out/'Tables'
    figs.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    records=[];audit=[];issues=[];channel_records=[];fused_source_records=[]
    for subject in SUBJECTS:
        source_path=sources[subject]
        with (ZipFile(source_path) if source_path.is_file() else FolderSource(source_path)) as src:
            files=list_run_members(src,subject)
            for run in RUNS:
                f=files[run]
                hemox=read_hemox(src,f['FUSED'][0],subject,run)
                windows=read_windows_channel_scores(src,f['WINDOWS'][0],subject,run)
                channel_records.append(windows)
                if len(windows)//2 != 147:
                    issues.append(dict(participant=subject,run=run,
                                       code='UNEXPECTED_SCHEDULED_WINDOW_COUNT',count=len(windows)//2))
                fused_source_records.append(source_fused_score_rows(hemox,subject,run))
                sp20=read_sp20(src,f['SP20'][0],subject,run)
                offset,verified,method,evidence=offsets.get((subject,run),
                    (0.,False,'assumed_first_record_zero_offset',''))
                quality=hemox.loc[hemox.quality_valid].copy()
                quality['hemox_aligned_t']=quality.hemox_t+offset
                quality=quality.loc[quality.hemox_aligned_t.between(0,300,inclusive='both')]
                pairs=pair_reference_first(sp20,quality,args.tolerance)
                for row in pairs:
                    row.update(participant=subject,run=run,
                        condition='Stationary' if run.startswith('S') else 'Controlled Motion',
                        offset_seconds=offset,alignment_verified=verified,alignment_method=method)
                    raw=row['native_spo2']
                    row['selected_spo2']=(float(np.clip(raw-args.spo2_offset,0,100))
                         if args.spo2_model=='m1' and np.isfinite(raw) and 0<=raw<=100 else raw)
                    row['selected_hr']=row['native_hr']
                records.extend(pairs)
                outside=int((hemox.fused_sqi.notna() & ~hemox.fused_sqi.between(0,1,inclusive='both')).sum())
                if outside:issues.append(dict(participant=subject,run=run,code='SQI_OUTSIDE_0_1',count=outside))
                audit.append(dict(participant=subject,run=run,source=str(source_path),
                    hemox_file=f['FUSED'][0],windows_file=f['WINDOWS'][0],sp20_file=f['SP20'][0],
                    scheduled_windows=len(windows)//2,
                    ch1_missing_sqi=int(windows[windows.channel.eq('CH1')].channel_sqi.isna().sum()),
                    ch2_missing_sqi=int(windows[windows.channel.eq('CH2')].channel_sqi.isna().sum()),
                    source_fused_rows=len(hemox),quality_valid_rows=int(hemox.quality_valid.sum()),
                    sp20_rows=len(sp20),time_matched_pairs=len(pairs),
                    fused_sqi_missing=int(hemox.fused_sqi.isna().sum()),
                    fused_sqi_outside_0_1=outside,offset_seconds=offset,
                    alignment_verified=verified,alignment_method=method,alignment_evidence=evidence))
                print(f'{subject}/{run}: valid fused={int(hemox.quality_valid.sum())}, time pairs={len(pairs)}, verified={verified}')
    paired=pd.DataFrame(records);audit=pd.DataFrame(audit)
    if len(audit)!=32 or audit[['participant','run']].duplicated().any():
        raise RuntimeError('Expected 32 unique runs (P03--P10 x S01/S02/NS01/NS02)')
    if paired.empty:raise ValueError('No reference/HEMOX candidate time pairs')
    evaluated={};exclusions=[]
    for typ in ('SpO2','HR'):
        data,ex=selected_measurements(paired,typ)
        if len(data)==0:raise ValueError(f'No valid SQI/reference pairs for {typ}')
        evaluated[typ]=data;exclusions.append(ex)
    window_scores=pd.concat(channel_records,ignore_index=True)
    fused_source=pd.concat(fused_source_records,ignore_index=True)
    if len(window_scores)!=2*int(audit.scheduled_windows.sum()):
        raise RuntimeError('Channel-window denominators are inconsistent')
    channel_table,fused_table,eligibility_summary=eligibility_tables(window_scores,fused_source,evaluated)
    interval_stats=calculate_bins(evaluated,args.bootstrap,args.seed)
    distribution=distribution_table(evaluated)
    all_errors=pd.concat(evaluated.values(),ignore_index=True)
    window_scores.to_csv(tables/'all_scheduled_window_channel_sqi_and_status.csv',index=False)
    fused_source.to_csv(tables/'all_native_fused_sqi_export_rows.csv',index=False)
    channel_table.to_csv(tables/'all_window_channel_sqi_acceptance_by_interval.csv',index=False)
    fused_table.to_csv(tables/'native_fused_sqi_selection_by_interval.csv',index=False)
    eligibility_summary.to_csv(tables/'low_sqi_coverage_and_missingness_QA.csv',index=False)
    paired.to_csv(tables/'native_reference_time_pairs_PROVISIONAL.csv',index=False)
    all_errors.to_csv(tables/'selected_model_sqi_absolute_errors_PROVISIONAL.csv',index=False)
    interval_stats.to_csv(tables/'agreement_by_sqi_3pp_5bpm_10bpm_PROVISIONAL.csv',index=False)
    distribution.to_csv(tables/'sqi_distribution_by_metric.csv',index=False)
    audit.to_csv(tables/'source_and_alignment_audit.csv',index=False)
    pd.DataFrame(exclusions).to_csv(tables/'metric_pair_exclusions_and_counts.csv',index=False)
    pd.DataFrame(issues,columns=['participant','run','code','count']).to_csv(tables/'QA_issues.csv',index=False)
    overall=[]
    for typ,df in evaluated.items():
        for limit in ((3.,) if typ=='SpO2' else (5.,10.)):
            passed=int(df.absolute_error.le(limit).sum())
            overall.append(dict(metric=typ,threshold=limit,unit='pp' if typ=='SpO2' else 'bpm',
                n_valid=len(df),n_agree=passed,agreement_pct=100*passed/len(df),
                min_sqi=float(df.fused_sqi.min()),n_sqi_below_0_75=int((df.fused_sqi<.75).sum())))
    summary=pd.DataFrame(overall)
    summary.to_csv(tables/'overall_agreement_criteria_PROVISIONAL.csv',index=False)
    font=style()
    figures=(make_scatter_figures(evaluated,figs,args.zoom_min)+
             make_agreement_figures(interval_stats,figs)+
             make_distribution_figures(distribution,figs)+
             make_eligibility_figures(window_scores,fused_source,evaluated,figs))
    status='All user-labeled offsets verified; independently review evidence.' if audit.alignment_verified.all() else \
           'PROVISIONAL: first-record zero-offset assumptions remain; independently verify device synchronization.'
    report=['# HEMOX SQI 0--1 and zoomed analysis','',
        f'Participants P03--P10; 32 five-minute recordings. {status}',
        f'Models: SpO2 {args.spo2_model.upper()} (frozen offset {args.spo2_offset:+.6f} pp if M1); HR M0.',
        f'Paired measurement time tolerance +/-{args.tolerance:g}s; match is SP20-first nearest unused quality-valid fused trend.',
        'Figures 1 full 0--1 and zoomed show the SAME paired observations; zoom clips their visual x-range ONLY, not the analysis.',
        'SQI bins: '+str(BIN_EDGES)+'; left-closed right-open except exact SQI 1.0 in final bin.',
        'Outcomes SpO2 <=3 pp; HR <=5 bpm (primary) and <=10 bpm (sensitivity). No estimates are interpolated or held.',
        'Fig. 3 distributes only valid reference pairs; Fig. 4 instead uses ALL scheduled WINDOWS CH SQI (accepted vs other status).',
        'Fig. 5 shows native FUSED SQI selection: source score -> native valid estimate -> valid reference-paired estimate. All categories are nested, not additive.',
        'Native FUSED_TRENDS quality-invalid rows can have fused_sqi=NaN; these missing values are counted, NOT scored as zero and NOT replaced with channel SQI.',
        'The WINDOW ch1/ch2 SQI are individual-channel scores, not equivalent to fused_sqi; compare distributions separately and note differing sample cadence.',
        'WINDOW status=accepted and spo2_eligible are separate export fields; Fig. 4 colors status, while both fields are preserved for QA.',
        'Fig. 4 denominator: 147 scheduled windows per run per channel when complete. Fig. 5 denominator: native exported FUSED_TRENDS rows, sampled differently.',
        'This script does not reconcile WINDOWS-versus-REPORT count discrepancies; review the separate dataset QA audit before using eligibility results as manuscript evidence.',
        'No figure uses held, imputed or reconstructed fused SQI for rejected observations. Zero n bins are empty, not 0% agreement.',
        'CIs use percentile 95% bootstrap resampling participants as entire clusters from all eight P03--P10; suppressed below 3 contributing participants.',
        'The zoomed panel cannot demonstrate behavior in sparse low-SQI bins; refer to full-range and n/P distribution.',
        '', '## Overall agreement',summary.to_string(index=False),
        '', '## Metric-specific validity',pd.DataFrame(exclusions).to_string(index=False),
        '', '## SQI-bin agreement with sample sizes and CIs',interval_stats.to_string(index=False),
        '', '## Low SQI / missing-score audit',eligibility_summary.to_string(index=False),
        '', '## All scheduled WINDOW channel-SQI / acceptance by interval',channel_table.to_string(index=False),
        '', '## Native FUSED SQI selection by interval',fused_table.to_string(index=False),
        '', '## Paired SQI distribution',distribution.to_string(index=False),
        '', '## Figures',*[f'- {x.name}' for x in figures],
        '', '## Source / offset audit',audit.to_string(index=False)]
    (out/'QA_report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    print('\n'+summary.to_string(index=False))
    print(f'\nGenerated {len(figures)} JPG figures; font={font}; output={out}')


if __name__=='__main__':
    main()