"""Generate only Heldout_SpO2_HR_AgreementBand_2Panel for HEMOX vs SP-20.

Requires pandas, numpy, matplotlib. Place beside extracted P01 ... P10 folders
or uniquely named participant ZIP archives. Supports --manifest with columns
participant,archive to select exact versions, or --data-dir path.

Default: all 10 required; P01/P02 development and P03-P10 held-out plotted
separately AND as a full N10 descriptive aggregate. Selected outputs: native
fused HR (M0), frozen SpO2 M1 offset (native SpO2 - 0.906754 pp, clipped).
Uses quality_valid, never the exported *_calibrated_* columns; independently
validates reference numerical ranges. Explorer-style provisional first-record
elapsed-time sync with nearest SP20-first, one-to-one tolerance +/-1 s.
An optional per-run offset CSV must have participant,run,offset_seconds,
method,evidence; offset is ADDED to HEMOX elapsed time. Verified requires
nonempty evidence. Independently unverified zero-offset plots are marked
PROVISIONAL; correlation is not an agreement or clinical accuracy measure.
"""
from __future__ import annotations
import argparse, io, re, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

RUNS=('S01','S02','NS01','NS02')
SUBJECTS=[f'P{i:02d}' for i in range(1,11)]
FROZEN_M1_OFFSET=0.906754
BLUE='#0072B2'; RED='#C62828'; GREY='#777F87'
ORANGE='#D55E00'; GREEN='#009E73'; LIGHTGREEN='#D9F2E8'
DPI=600
# User-adjustable agreement tolerances (absolute HEMOX minus SP-20 error).
BANDS={'SpO2':(3.0,None),'HR':(5.0,10.0)}

def args_parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--out-dir',type=Path,default=None)
    p.add_argument('--pairs-csv',type=Path,default=None,help='Optional prior exported native_and_selected_reference_pairs.csv; replot original real paired observations')
    p.add_argument('--spo2-band',type=float,default=3.0,help='Primary SpO2 agreement limit in percentage points (default 3)')
    p.add_argument('--hr-band',type=float,default=5.0,help='Primary HR agreement limit in bpm (default 5)')
    p.add_argument('--hr-sensitivity',type=float,default=10.0,help='Secondary HR sensitivity limit in bpm (default 10)')
    p.add_argument('--manifest',type=Path,default=None,help='Optional participant,archive CSV to select ZIP versions')
    p.add_argument('--offsets-csv',type=Path,default=None,help='Optional verified per-run time offset CSV')
    p.add_argument('--tolerance',type=float,default=1.0)
    p.add_argument('--spo2-offset',type=float,default=FROZEN_M1_OFFSET)
    p.add_argument('--spo2-model',choices=['m1','m0'],default='m1')
    p.add_argument('--allow-incomplete',action='store_true',help='Allow missing participants/runs, prominently documented in QA')
    p.add_argument('--show-calibration-on-heldout',action='store_true',help='Overlay P01-P02 lightly on validation plot without including them in its metrics')
    return p.parse_args()

def norm_sub(s):
    m=re.fullmatch(r'P\s*0*(\d{1,2})',str(s).strip(),re.I)
    return f'P{int(m.group(1)):02d}' if m and 1 <= int(m.group(1)) <= 10 else None

def run_of(name,subject):
    u=str(name).upper()
    m=re.search(r'P0*'+str(int(subject[1:]))+r'(NS0[12]|S0[12])(?=[_./\\]|$)',u)
    return m.group(1) if m else None

def choose_sources(root,manifest):
    explicit={}
    if manifest:
        x=pd.read_csv(manifest)
        if not {'participant','archive'}.issubset(x.columns):raise ValueError('Manifest requires participant,archive')
        for row in x.itertuples(index=False):
            sub=norm_sub(row.participant)
            if not sub or sub in explicit:raise ValueError(f'Invalid/duplicate manifest participant {row.participant}')
            path=Path(str(row.archive)).expanduser()
            if not path.is_absolute():path=manifest.resolve().parent/path
            if not path.is_file():raise FileNotFoundError(path)
            explicit[sub]=('zip',path)
    found={}
    for sub in SUBJECTS:
        if sub in explicit:found[sub]=explicit[sub];continue
        dirs=[p for p in root.iterdir() if p.is_dir() and p.name.upper()==sub]
        if len(dirs)>1:raise ValueError(f'Conflicting {sub} directories: {dirs}')
        if dirs:found[sub]=('directory',dirs[0]);continue
        zips=[p for p in root.iterdir() if p.is_file() and p.suffix.lower()=='.zip' and re.fullmatch(re.escape(sub)+r'(?:\(\d+\))?',p.stem,re.I)]
        exact=[p for p in zips if p.stem.upper()==sub]
        if len(exact)==1:found[sub]=('zip',exact[0])
        elif len(zips)==1:found[sub]=('zip',zips[0])
        elif len(zips)>1:raise ValueError(f'Ambiguous {sub} ZIPs: {zips}. Use --manifest.')
    return found

def index_files(kind,path,sub):
    if kind=='directory':names=[p for p in path.rglob('*.csv')]
    else:
        with zipfile.ZipFile(path) as z:names=[n for n in z.namelist() if n.lower().endswith('.csv')]
    out={run:{'he':[],'sp':[]} for run in RUNS}
    for n in names:
        base=Path(n).name.upper();run=run_of(n,sub)
        if run is None:continue
        if base.endswith('_FUSED_TRENDS.CSV') and 'CALIBRATED' not in base:out[run]['he'].append(n)
        if base.startswith('SP20'):out[run]['sp'].append(n)
    return out

def read_bytes(kind,source,item):
    if kind=='directory':return Path(item).read_bytes()
    with zipfile.ZipFile(source) as z:return z.read(item)

def read_he(data):
    r=pd.read_csv(io.BytesIO(data),low_memory=False)
    needed={'elapsed_s','fused_spo2_percent','fused_hr_bpm','quality_valid'}
    if not needed.issubset(r):raise ValueError(f'FUSED_TRENDS missing {needed-set(r.columns)}')
    t=pd.to_numeric(r.elapsed_s,errors='coerce')
    if t.isna().any() or not t.is_monotonic_increasing or t.duplicated().any():raise ValueError('HEMOX nonmonotone/duplicated elapsed_s')
    valid=r.quality_valid.astype(str).str.strip().str.lower().isin(('true','yes','1','valid'))
    return pd.DataFrame({'ht':t.to_numpy(float),'valid':valid.to_numpy(bool),
        'spo2':pd.to_numeric(r.fused_spo2_percent,errors='coerce').to_numpy(float),
        'hr':pd.to_numeric(r.fused_hr_bpm,errors='coerce').to_numpy(float)})

def read_sp(data):
    r=pd.read_csv(io.BytesIO(data),skiprows=1,skipinitialspace=True,dtype=str)
    r.columns=[str(c).strip().lstrip('#').strip() for c in r.columns]
    names={'Time','SpO2(%)','PR(bpm)'}
    if not names.issubset(r):raise ValueError(f'SP20 missing {names-set(r.columns)}')
    dt=pd.to_datetime(r['Time'].str.strip(),format='%m/%d/%Y %H:%M:%S',errors='coerce')
    if dt.isna().any() or not dt.is_monotonic_increasing or dt.duplicated().any():raise ValueError('SP20 invalid/duplicate time')
    return pd.DataFrame({'rt':(dt-dt.iloc[0]).dt.total_seconds().to_numpy(float),
       'timestamp':dt.astype(str),'spo2':pd.to_numeric(r['SpO2(%)'],errors='coerce'),
       'hr':pd.to_numeric(r['PR(bpm)'],errors='coerce')})

def time_offsets(path):
    if path is None:return {}
    r=pd.read_csv(path);need={'participant','run','offset_seconds','method','evidence'}
    if not need.issubset(r):raise ValueError(f'Offsets CSV must contain {need}')
    d={}
    for x in r.itertuples(index=False):
        k=(norm_sub(x.participant),str(x.run).upper().strip())
        if k[0] not in SUBJECTS or k[1] not in RUNS or k in d:raise ValueError(f'Invalid/duplicate offset {k}')
        off=float(x.offset_seconds)
        if not np.isfinite(off):raise ValueError(f'Nonfinite offset {k}')
        method=str(x.method).strip().lower();evidence=str(x.evidence).strip()
        if method=='verified' and (not evidence or evidence.lower() in ('nan','none')):raise ValueError(f'No evidence for {k}')
        d[k]=(off,method,evidence)
    return d

def pair(sp,he,tolerance,offset):
    # SP20-first one-to-one nearest elapsed time, inside common recording overlap.
    s=sp.reset_index(drop=True);h=he.reset_index(drop=True).copy()
    h['aligned']=h.ht+offset
    if s.empty or h.empty:return []
    low=max(s.rt.iloc[0],h.aligned.iloc[0]);high=min(s.rt.iloc[-1],h.aligned.iloc[-1])
    if high<low:return []
    ht=h.aligned.to_numpy(float);used=np.zeros(len(h),dtype=bool);rows=[]
    for x in s.itertuples(index=False):
        if not low<=x.rt<=high:continue
        at=int(np.searchsorted(ht,x.rt))
        # Search candidates near the reference point in sorted timeline.
        left=max(0,at-2);right=min(len(ht),at+3)
        choices=[i for i in range(left,right) if not used[i] and low<=ht[i]<=high]
        if not choices:continue
        i=min(choices,key=lambda j:(abs(ht[j]-x.rt),j))
        if abs(ht[i]-x.rt)>tolerance:continue
        used[i]=True;y=h.iloc[i]
        rows.append(dict(sp20_elapsed_s=x.rt,sp20_timestamp=x.timestamp,hemox_elapsed_s=y.ht,
           delta_t_s=y.aligned-x.rt,quality_valid=bool(y.valid),
           sp20_spo2=x.spo2,sp20_hr=x.hr,hemox_spo2_m0=y.spo2,hemox_hr_m0=y.hr))
    return rows

def setup():
    installed={f.name for f in font_manager.fontManager.ttflist}
    font=next((x for x in ['Arial Narrow','Arial','Liberation Sans','DejaVu Sans'] if x in installed),'DejaVu Sans')
    plt.rcParams.update({'font.family':font,'font.weight':'bold','axes.labelweight':'bold',
     'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'axes.linewidth':1.3,
     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':DPI})

def stats(x,y):
    if len(x)<3 or np.std(x)==0 or np.std(y)==0:return None
    slope,inter=np.polyfit(x,y,1)
    e=y-x
    return dict(n=len(x),r=float(np.corrcoef(x,y)[0,1]),slope=float(slope),intercept=float(inter),
        bias=float(e.mean()),mae=float(np.mean(abs(e))),rmse=float(np.sqrt(np.mean(e**2))))

def figure(df,metric,cohort,out,overlay_cal=False):
    ref,pred=('sp20_hr','selected_hr') if metric=='HR' else ('sp20_spo2','selected_spo2')
    x=df[ref].to_numpy(float);y=df[pred].to_numpy(float)
    good=np.isfinite(x)&np.isfinite(y);x=x[good];y=y[good]
    st=stats(x,y)
    if st is None:raise ValueError(f'Not enough varied data for {metric}/{cohort}')
    fig,ax=plt.subplots(figsize=(7.2,3.8))
    if overlay_cal and cohort=='Heldout':
        cal=df.attrs['calibration'];cc=cal[[ref,pred]].dropna()
        ax.scatter(cc[ref],cc[pred],s=12,color=GREY,alpha=.13,edgecolors='none',rasterized=True,label='P01–P02 (not fitted)')
    ax.scatter(x,y,s=17,facecolors='none',edgecolors=BLUE,linewidths=.55,alpha=.46,rasterized=True)
    # Identity (perfect agreement) and OLS fit are separate visual concepts.
    minimum=float(min(x.min(),y.min()));maximum=float(max(x.max(),y.max()))
    span=max(maximum-minimum,1 if metric=='HR' else .25)
    lower=minimum-.07*span;upper=maximum+.07*span
    xx=np.array([lower,upper]);ax.plot(xx,xx,color=GREY,ls='--',lw=1.1,zorder=1,label='Identity (y = x)')
    ax.plot(xx,st['slope']*xx+st['intercept'],color=RED,lw=2.0,zorder=4,label='Linear fit')
    ax.set_xlim(lower,upper);ax.set_ylim(lower,upper);ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel('SP-20 reference HR (bpm)' if metric=='HR' else r'SP-20 reference SpO$_2$ (%)')
    ax.set_ylabel('HEMOX fused HR (bpm)' if metric=='HR' else r'HEMOX fused SpO$_2$ (%)')
    ax.text(.03,.95,f'R = {st["r"]:.4f}   (n = {st["n"]:,})',transform=ax.transAxes,
       fontsize=11,weight='bold',va='top',bbox=dict(facecolor='white',edgecolor='none',alpha=.85))
    sign='+' if st['intercept']>=0 else '−'
    txt=f'Fit: y = {st["slope"]:.3f}x {sign} {abs(st["intercept"]):.3f}'
    ax.text(.03,.035,txt,transform=ax.transAxes,fontsize=10.5,weight='bold',
       bbox=dict(facecolor='white',edgecolor='none',alpha=.88))
    ax.grid(True,lw=.55,color='#CFD5DC',alpha=.5);ax.set_axisbelow(True)
    ax.legend(loc='lower right',fontsize=8.5,facecolor='white',edgecolor='#B0BAC3')
    fig.subplots_adjust(left=.16,right=.98,bottom=.22,top=.97)
    fig.savefig(out,dpi=DPI,facecolor='white',pil_kwargs={'quality':96,'subsampling':0});plt.close(fig)
    return {'cohort':cohort,'variable':metric,**st}

def combined_figure(df,cohort,out):
    """Two-panel manuscript companion, calculated from the same selected pairs."""
    fig,axes=plt.subplots(1,2,figsize=(10.2,4.2))
    for ax,metric,letter in zip(axes,('SpO2','HR'),'ab'):
        ref,pred,valid=(('sp20_spo2','selected_spo2','spo2_pair_valid') if metric=='SpO2' else
                       ('sp20_hr','selected_hr','hr_pair_valid'))
        subset=df[df[valid]];x=subset[ref].to_numpy(float);y=subset[pred].to_numpy(float)
        st=stats(x,y)
        if st is None:continue
        lo=min(float(x.min()),float(y.min()));hi=max(float(x.max()),float(y.max()))
        span=max(hi-lo,1 if metric=='HR' else .25);lo-=.07*span;hi+=.07*span
        v=np.array([lo,hi]);ax.scatter(x,y,s=13,facecolors='none',edgecolors=BLUE,
            lw=.5,alpha=.38,rasterized=True)
        ax.plot(v,v,'--',color=GREY,lw=1.0)
        ax.plot(v,st['slope']*v+st['intercept'],color=RED,lw=1.8)
        ax.set(xlim=(lo,hi),ylim=(lo,hi))
        ax.set_aspect('equal',adjustable='box')
        ax.set_xlabel('SP-20 SpO$_2$ (%)' if metric=='SpO2' else 'SP-20 HR (bpm)')
        ax.set_ylabel('HEMOX fused SpO$_2$ (%)' if metric=='SpO2' else 'HEMOX fused HR (bpm)')
        ax.text(.025,.97,f'({letter})  R = {st["r"]:.4f}; n = {st["n"]:,}',transform=ax.transAxes,
            va='top',fontsize=10,weight='bold',bbox=dict(facecolor='white',alpha=.82,edgecolor='none'))
        op='+' if st['intercept']>=0 else '−'
        ax.text(.025,.035,f'y = {st["slope"]:.3f}x {op} {abs(st["intercept"]):.3f}',
            transform=ax.transAxes,fontsize=9,weight='bold',
            bbox=dict(facecolor='white',alpha=.88,edgecolor='none'))
        ax.grid(True,color='#CFD5DC',lw=.5,alpha=.5);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.10,right=.99,bottom=.19,top=.97,wspace=.36)
    fig.savefig(out,dpi=DPI,facecolor='white',pil_kwargs={'quality':96,'subsampling':0});plt.close(fig)

def _limits(a, *, metric, ref_axis):
    """Independent x/y limits, intentionally NOT symmetrical and no equal aspect."""
    a=np.asarray(a,float);a=a[np.isfinite(a)]
    low=float(a.min()); high=float(a.max())
    span=max(high-low,1.0 if metric=='HR' else .3)
    margin=span*(.065 if ref_axis else .095)
    return (low-margin,high+margin)


def _thresholds(metric):
    return BANDS[metric]


def _metric_fields(metric):
    if metric=='SpO2':
        return 'sp20_spo2','selected_spo2','spo2_pair_valid', 'SpO$_2$ (%)', 'pp'
    return 'sp20_hr','selected_hr','hr_pair_valid','HR (bpm)','bpm'


def _band_statistics(x,y,metric):
    x=np.asarray(x,float); y=np.asarray(y,float)
    primary,secondary=_thresholds(metric)
    error=y-x
    base=stats(x,y)
    n=len(error)
    if not n or base is None:return None
    results={**base,'primary_threshold':primary,'primary_within':int((abs(error)<=primary).sum()),
        'primary_pct':100*float(np.mean(abs(error)<=primary)),
        'primary_outside':int((abs(error)>primary).sum()),
        'secondary_threshold':secondary,'secondary_within':np.nan,'secondary_pct':np.nan}
    if secondary is not None:
        results['secondary_within']=int((abs(error)<=secondary).sum())
        results['secondary_pct']=100*float(np.mean(abs(error)<=secondary))
    return results


def _draw_band(ax,x,y,metric,letter=None):
    """One actual paired sample per point; band statistics are descriptive."""
    primary,secondary=_thresholds(metric)
    st=_band_statistics(x,y,metric)
    if st is None:raise ValueError(f'Not enough variable data: {metric}')
    xlim=_limits(x,metric=metric,ref_axis=True)
    ylim=_limits(y,metric=metric,ref_axis=False)
    xx=np.linspace(xlim[0],xlim[1],400)
    # The true +/- error envelope is plotted relative to the identity, not OLS.
    ax.fill_between(xx,xx-primary,xx+primary,color=GREEN,alpha=.105,zorder=0)
    if secondary is not None:
        ax.plot(xx,xx-secondary,color=GREY,lw=.95,ls=':',zorder=1)
        ax.plot(xx,xx+secondary,color=GREY,lw=.95,ls=':',zorder=1,
                label=f'±{secondary:g} bpm (sensitivity)')
    ax.plot(xx,xx-primary,color=GREEN,lw=1.25,ls='--',zorder=2)
    ax.plot(xx,xx+primary,color=GREEN,lw=1.25,ls='--',zorder=2,
            label=f'±{primary:g} {"pp" if metric=="SpO2" else "bpm"} agreement band')
    ax.plot(xx,xx,color='#33383D',lw=1.25,ls='-',zorder=2,label='Identity: y = x')
    ax.plot(xx,st['slope']*xx+st['intercept'],color=RED,lw=1.65,zorder=3,label='Linear fit')
    within=abs(y-x)<=primary
    ax.scatter(x[within],y[within],s=19,marker='o',facecolor=BLUE,edgecolor='none',
               alpha=.28,rasterized=True,zorder=4,label=f'Within ±{primary:g}')
    if (~within).any():
        ax.scatter(x[~within],y[~within],s=28,marker='o',facecolor=ORANGE,
                   edgecolor='white',linewidth=.4,alpha=.92,rasterized=True,zorder=5,
                   label=f'Outside ±{primary:g}')
    ax.set_xlim(*xlim);ax.set_ylim(*ylim)
    # Deliberately NO ax.set_aspect('equal') and NO forced matching x/y limits.
    ax.set_xlabel('SP-20 reference '+('SpO$_2$ (%)' if metric=='SpO2' else 'HR (bpm)'))
    ax.set_ylabel('HEMOX fused '+('SpO$_2$ (%)' if metric=='SpO2' else 'HR (bpm)'))
    if letter:ax.text(.012,.975,f'({letter})',transform=ax.transAxes,ha='left',va='top',weight='bold',fontsize=12)
    stats_text=(f'R = {st["r"]:.4f}; n = {st["n"]:,}  |  '
                f'within ±{primary:g}: {st["primary_pct"]:.2f}% '
                f'({st["primary_within"]:,}/{st["n"]:,})')
    if secondary is not None:
        stats_text += f'\nwithin ±{secondary:g}: {st["secondary_pct"]:.2f}% ({st["secondary_within"]:,}/{st["n"]:,})'
    ax.text(.985,.035,stats_text,transform=ax.transAxes,ha='right',va='bottom',
            fontsize=9.0,weight='bold',bbox=dict(boxstyle='round,pad=.28',facecolor='white',edgecolor='#E5E8EB',alpha=.92),zorder=8)
    ax.grid(True,color='#CCD4DD',lw=.53,alpha=.53);ax.set_axisbelow(True)
    return st


def agreement_band_figure(df,metric,cohort,out):
    ref,pred,valid,_,_=_metric_fields(metric)
    sub=df.loc[df[valid],[ref,pred]].dropna()
    x=sub[ref].to_numpy(float); y=sub[pred].to_numpy(float)
    fig,ax=plt.subplots(figsize=(7.2,4.35))
    st=_draw_band(ax,x,y,metric)
    handles,labels=ax.get_legend_handles_labels()
    leg=fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.99),
                   ncol=3,fontsize=8.4,framealpha=.97,edgecolor='#ABB4BD',columnspacing=1.0)
    for t in leg.get_texts():t.set_fontweight('bold')
    fig.subplots_adjust(left=.145,right=.985,bottom=.18,top=.79)
    fig.savefig(out,dpi=DPI,facecolor='white',pil_kwargs={'quality':96,'subsampling':0});plt.close(fig)
    return {'cohort':cohort,'variable':metric,**st}


def heldout_rmse_agreement_figure(df,metric,out):
    """Additional RMSE-centered publication view, leaving ALL original plots intact.

    No identity or regression line; no R, n, or within-band statistic annotation.
    Agreement bands are still geometrically defined about the unplotted y=x line.
    The selected-model RMSE is independently computed from valid held-out pairs.
    """
    ref,pred,valid,_,unit=_metric_fields(metric)
    sub=df.loc[df[valid],[ref,pred]].dropna()
    x=sub[ref].to_numpy(float); y=sub[pred].to_numpy(float)
    if len(x)<3:
        raise ValueError(f'Insufficient held-out {metric} reference pairs')
    primary,secondary=_thresholds(metric)
    residual=y-x
    rmse=float(np.sqrt(np.mean(np.square(residual))))
    within=np.abs(residual)<=primary
    xlim=_limits(x,metric=metric,ref_axis=True)
    ylim=_limits(y,metric=metric,ref_axis=False)
    xx=np.linspace(xlim[0],xlim[1],400)

    fig,ax=plt.subplots(figsize=(7.2,4.35))
    ax.fill_between(xx,xx-primary,xx+primary,color=GREEN,alpha=.105,zorder=0)
    if secondary is not None:
        ax.plot(xx,xx-secondary,color=GREY,lw=.95,ls=':',zorder=1)
        ax.plot(xx,xx+secondary,color=GREY,lw=.95,ls=':',zorder=1,
                label=f'±{secondary:g} bpm (sensitivity)')
    ax.plot(xx,xx-primary,color=GREEN,lw=1.25,ls='--',zorder=2)
    ax.plot(xx,xx+primary,color=GREEN,lw=1.25,ls='--',zorder=2,
            label=f'±{primary:g} {unit} agreement band')
    # Do NOT plot either identity y=x or a fitted regression in this variant.
    ax.scatter(x[within],y[within],s=19,marker='o',facecolor=BLUE,
               edgecolor='none',alpha=.28,rasterized=True,zorder=4,
               label=f'Within ±{primary:g}')
    if (~within).any():
        ax.scatter(x[~within],y[~within],s=28,marker='o',facecolor=ORANGE,
                   edgecolor='white',linewidth=.4,alpha=.92,rasterized=True,
                   zorder=5,label=f'Outside ±{primary:g}')
    ax.set_xlim(*xlim);ax.set_ylim(*ylim)  # independent axes, unchanged from originals
    ax.set_xlabel('SP-20 reference '+('SpO$_2$ (%)' if metric=='SpO2' else 'HR (bpm)'))
    ax.set_ylabel('HEMOX fused '+('SpO$_2$ (%)' if metric=='SpO2' else 'HR (bpm)'))
    ax.grid(True,color='#CCD4DD',lw=.53,alpha=.53);ax.set_axisbelow(True)

    handles,labels=ax.get_legend_handles_labels()
    leg=fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.99),
                   ncol=3,fontsize=8.4,framealpha=.97,edgecolor='#ABB4BD',
                   columnspacing=1.0)
    for txt in leg.get_texts():txt.set_fontweight('bold')
    # The ONLY statistical annotation in these new figures is RMSE.
    fig.text(.985,.855,f'RMSE = {rmse:.3f} {unit}',ha='right',va='center',
             fontsize=12,weight='bold',color='#253348',
             bbox=dict(boxstyle='round,pad=.28',facecolor='white',
                       edgecolor='#AEB8C2',alpha=.98))
    fig.subplots_adjust(left=.145,right=.985,bottom=.18,top=.76)
    fig.savefig(out,dpi=DPI,facecolor='white',
                pil_kwargs={'quality':96,'subsampling':0})
    plt.close(fig)
    return dict(variable=metric,rmse=rmse,n=len(x),output=str(out))


def agreement_band_combined(df,cohort,out):
    fig,axes=plt.subplots(1,2,figsize=(10.6,4.5))
    results=[]
    for ax,metric,letter in zip(axes,['SpO2','HR'],'ab'):
        ref,pred,valid,_,_=_metric_fields(metric)
        sub=df.loc[df[valid],[ref,pred]].dropna()
        st=_draw_band(ax,sub[ref].to_numpy(float),sub[pred].to_numpy(float),metric,letter=letter)
        results.append({'cohort':cohort,'variable':metric,**st})
        ax.xaxis.label.set_size(10.5);ax.yaxis.label.set_size(10.5)
        ax.tick_params(axis='both',labelsize=9.6)
        for txt in ax.texts:
            if txt.get_text().startswith('R ='):txt.set_fontsize(7.8)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend=[Line2D([0],[0],color=BLUE,marker='o',ls='none',markersize=5,label='Within primary limit'),
            Line2D([0],[0],color=ORANGE,marker='o',ls='none',markersize=5,label='Outside primary limit'),
            Line2D([0],[0],color='#33383D',lw=1.2,label='Identity'),
            Line2D([0],[0],color=RED,lw=1.6,label='Linear fit'),
            Line2D([0],[0],color=GREEN,ls='--',lw=1.2,label=f'Primary ±{BANDS["SpO2"][0]:g} pp / ±{BANDS["HR"][0]:g} bpm'),
            Line2D([0],[0],color=GREY,ls=':',lw=1.1,label=f'HR ±{BANDS["HR"][1]:g} bpm')]
    fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.994),ncol=3,
               fontsize=9.1,framealpha=.96,edgecolor='#ADB7C0')
    fig.subplots_adjust(left=.09,right=.99,bottom=.18,top=.75,wspace=.33)
    fig.savefig(out,dpi=DPI,facecolor='white',pil_kwargs={'quality':96,'subsampling':0});plt.close(fig)
    return results


def _per_participant_agreement(df):
    rows=[]
    for participant,group in df.groupby('participant',sort=True):
        for metric in ['SpO2','HR']:
            ref,pred,valid,_,_=_metric_fields(metric)
            sub=group.loc[group[valid],[ref,pred]].dropna()
            x=sub[ref].to_numpy(float);y=sub[pred].to_numpy(float)
            st=_band_statistics(x,y,metric) if len(x)>=3 else None
            if st:rows.append({'participant':participant,'variable':metric,**st})
    return pd.DataFrame(rows)

def main():
    a=args_parser()
    if min(a.spo2_band,a.hr_band,a.hr_sensitivity)<=0 or a.hr_sensitivity<a.hr_band:
        raise ValueError('Agreement limits must be positive and secondary HR limit >= primary HR limit')
    BANDS['SpO2']=(a.spo2_band,None)
    BANDS['HR']=(a.hr_band,a.hr_sensitivity)
    root=a.data_dir.expanduser().resolve();out=(a.out_dir or root/'Study_HEMOX_Correlation_outputs').resolve()
    if not root.is_dir() or not (a.tolerance>0) or not np.isfinite(a.spo2_offset):
        raise ValueError('Invalid input directory/tolerance/offset')
    out.mkdir(parents=True,exist_ok=True);figdir=out/'Figures';figdir.mkdir(exist_ok=True)
    setup()

    missing=[]
    if a.pairs_csv is not None:
        d=pd.read_csv(a.pairs_csv,low_memory=False)
        needed={'participant','run','quality_valid','sp20_spo2','sp20_hr','hemox_spo2_m0','hemox_hr_m0'}
        if not needed.issubset(d):
            raise ValueError(f'Prior pairs file lacks {needed-set(d.columns)}')
        if d.empty:
            raise ValueError('Provided paired measurement CSV is empty')
        d['participant']=d['participant'].map(norm_sub)
        if d.participant.isna().any():
            raise ValueError('Invalid participant ID in prior paired measurements')
        d['quality_valid']=d.quality_valid.astype(str).str.lower().isin(['true','yes','1','valid'])
        missing=[s for s in SUBJECTS if s not in set(d.participant)]
        if missing and not a.allow_incomplete:
            raise FileNotFoundError(f'Prior pair CSV lacks participants {missing}; pass --allow-incomplete only intentionally')
        d['split']=np.where(d.participant.isin(SUBJECTS[:2]),'Calibration','Heldout')
        qa=d.groupby(['participant','run'],as_index=False).size().rename(columns={'size':'reference_pairs'})
        qa['status']='reused_previously_exported_pairs'; qa['source']=str(a.pairs_csv)
        qa.to_csv(out/'run_level_QA.csv',index=False)
        incomplete=bool(missing or len(qa)!=40)
        if incomplete and not a.allow_incomplete:
            raise RuntimeError('Incomplete participant/run set in prior pair CSV. Use --allow-incomplete only intentionally.')
        provenance=str(a.pairs_csv)
    else:
        sources=choose_sources(root,a.manifest); offsets=time_offsets(a.offsets_csv)
        missing=[x for x in SUBJECTS if x not in sources]
        if missing and not a.allow_incomplete:
            raise FileNotFoundError(f'Missing {missing}. Supply P01–P10 directories/ZIPs, or explicitly use --allow-incomplete.')
        allpairs=[]; qa=[]
        for sub,(kind,src) in sources.items():
            index=index_files(kind,src,sub)
            for run in RUNS:
                files=index[run]
                q={'participant':sub,'run':run,'source':str(src),'kind':kind,'hemox_files':len(files['he']),'sp20_files':len(files['sp'])}
                if len(files['he'])!=1 or len(files['sp'])!=1:
                    q.update(status='missing_or_duplicate_source',reference_pairs=0); qa.append(q); continue
                try:
                    he=read_he(read_bytes(kind,src,files['he'][0]))
                    sp=read_sp(read_bytes(kind,src,files['sp'][0]))
                    he=he.loc[he.ht.between(0,300) & he.valid].reset_index(drop=True)
                    sp=sp.loc[sp.rt.between(0,300)].reset_index(drop=True)
                    off,method,evidence=offsets.get((sub,run),(0.0,'provisional_first_record',''))
                    he=he.loc[(he.ht+off>=0) & (he.ht+off<=300)].reset_index(drop=True)
                    pairs=pair(sp,he,a.tolerance,off)
                    for row in pairs:
                        row.update(participant=sub,run=run,split='Calibration' if sub in ('P01','P02') else 'Heldout',
                                   offset_seconds=off,alignment_method=method)
                    allpairs += pairs
                    q.update(status='ok',hemox_rows=len(he),sp20_rows=len(sp),reference_pairs=len(pairs),
                             alignment_method=method,offset_seconds=off,quality_valid_pairs=sum(p['quality_valid'] for p in pairs))
                except Exception as exc:
                    q.update(status=f'ERROR: {exc}',reference_pairs=0)
                qa.append(q)
        qa=pd.DataFrame(qa); qa.to_csv(out/'run_level_QA.csv',index=False)
        incomplete=bool(missing or qa.status.ne('ok').any() or set(qa.participant)!=set(SUBJECTS))
        if incomplete and not a.allow_incomplete:
            raise RuntimeError('Incomplete sources/runs; see run_level_QA.csv. Use --allow-incomplete only intentionally.')
        if not allpairs:
            raise ValueError('No timestamp-paired measurements; see QA CSV')
        d=pd.DataFrame(allpairs)
        provenance='raw HEMOX FUSED_TRENDS and SP-20 exports'

    d['selected_spo2']=np.clip(d.hemox_spo2_m0-a.spo2_offset,0,100) if a.spo2_model=='m1' else d.hemox_spo2_m0
    d['selected_hr']=d.hemox_hr_m0
    spgood=d.quality_valid & d.sp20_spo2.between(0,100) & d.selected_spo2.between(0,100)
    hrgood=d.quality_valid & d.sp20_hr.gt(0) & d.selected_hr.gt(0) & np.isfinite(d.sp20_hr) & np.isfinite(d.selected_hr)
    d['spo2_pair_valid']=spgood; d['hr_pair_valid']=hrgood
    d.to_csv(out/'native_and_selected_reference_pairs.csv',index=False)

    heldout=d[d.participant.isin(SUBJECTS[2:])].copy()
    if heldout.empty:
        raise ValueError('No held-out (P03–P10) paired observations found')
    fig_path = figdir/'Heldout_SpO2_HR_AgreementBand_2Panel.jpg'
    results = agreement_band_combined(heldout,'Heldout',fig_path)
    pd.DataFrame(results).to_csv(out/'heldout_agreement_band_metrics.csv',index=False,float_format='%.8f')
    _per_participant_agreement(heldout).to_csv(out/'participant_agreement_band_metrics.csv',index=False,float_format='%.8f')
    note=('Generated only Heldout_SpO2_HR_AgreementBand_2Panel. '          'P03–P10 held-out cohort; SpO2 model {model} (offset {offset:.6f} pp if M1), HR M0. '          'Primary agreement band: SpO2 ±{sband:g} pp, HR ±{hband:g} bpm; HR sensitivity ±{hsens:g} bpm. '          'Source: {src}. Missing participants: {missing}; incomplete: {incomplete}. '          'Band percentages use eligible finite paired observations.').format(
              model=a.spo2_model,offset=a.spo2_offset,sband=a.spo2_band,hband=a.hr_band,
              hsens=a.hr_sensitivity,src=provenance,missing=missing,incomplete=bool(incomplete))
    (out/'QA_readme.txt').write_text(note+'\n',encoding='utf-8')
    print(note)
    print('Saved:', fig_path)
    print(pd.DataFrame(results).to_string(index=False))

if __name__=='__main__':main()
