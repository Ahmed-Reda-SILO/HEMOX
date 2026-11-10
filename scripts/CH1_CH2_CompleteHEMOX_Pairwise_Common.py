
from __future__ import annotations
import argparse, re, warnings
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

PARTICIPANTS=tuple(f'P{i:02d}' for i in range(1,11))
DEV=('P01','P02'); HELDOUT=tuple(f'P{i:02d}' for i in range(3,11))
RUNS=('S01','S02','NS01','NS02')
FS=50.0; WINDOW_S=8.0; DURATION_S=300.0; TOL_S=1.0
IMU_SUPPORT_S=.5; BASE_FRAC=.40; BASE_K=5.0
ACC_FLOOR,ACC_CEIL=.002,.010; GYR_FLOOR,GYR_CEIL=.15,.75
CAL_SETTLE,CAL_SEG,CAL_SEARCH=2.0,5.0,30.0; CAL_ON_K,CAL_OFF_K=6.0,3.0
ACC_MAX,GYR_MAX=.030,3.5; RELEASE_RATIO=.80
ONSET_S,RELEASE_S,PAD_S=.15,.30,.20; MOTION_REJECT=26.0
SQI_FALLBACK=.70; SPO2_CONS=3.0; HR_CONS=5.0
SPO2_WITHIN=3.0; HR_WITHIN=5.0; BOOT_N=2000; SEED=1986; DPI=600
BLUE='#0072B2'; RED='#C62828'; GREY='#777F87'; GREEN='#009E73'; ORANGE='#D55E00'; GRID='#D9DEE3'; BLACK='#202124'
METHODS=('CH1','CH2','Simple average','SQI-only fusion','Motion-only fusion','Fixed thresholds','Adaptive thresholds','Complete HEMOX')

def cli():
    h=Path(__file__).resolve().parent
    p=argparse.ArgumentParser()
    p.add_argument('--data-dir',type=Path,default=h)
    p.add_argument('--out-dir',type=Path,default=h/'Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common_output')
    p.add_argument('--duration',type=float,default=DURATION_S)
    p.add_argument('--tolerance',type=float,default=TOL_S)
    p.add_argument('--offsets-csv',type=Path,default=None)
    p.add_argument('--bootstrap',type=int,default=BOOT_N)
    p.add_argument('--seed',type=int,default=SEED)
    p.add_argument('--no-method-specific-m1',action='store_true')
    return p.parse_args()

def set_style():
    fs={f.name for f in font_manager.fontManager.ttflist}
    fam=next((f for f in ('Arial Narrow','Arial','Liberation Sans','DejaVu Sans') if f in fs),'DejaVu Sans')
    plt.rcParams.update({'font.family':fam,'font.weight':'bold','axes.labelweight':'bold','axes.labelsize':10,'axes.linewidth':1.0,'axes.spines.top':False,'axes.spines.right':False,'xtick.labelsize':8,'ytick.labelsize':8.5,'figure.facecolor':'white','savefig.facecolor':'white'})

def safe(s): return re.sub(r'[^A-Za-z0-9]+','_',s).strip('_').lower()
def run_token(text,subject=None):
    s=str(text).upper()
    if subject:
        m=re.search(rf'P0*{int(subject[1:])}(NS01|NS02|S01|S02)',s)
        if m:return m.group(1)
    m=re.search(r'(NS01|NS02|S01|S02)',s); return m.group(1) if m else None

class Source:
    def __init__(self,root,subject):
        self.subject=subject; self.z=None; root=root.resolve(); folder=root/subject; zp=root/f'{subject}.zip'
        if folder.is_dir():
            self.kind='dir'; self.base=folder; self.names=[str(p.relative_to(folder)).replace('\\','/') for p in folder.rglob('*') if p.is_file()]
        elif zp.is_file():
            self.kind='zip'; self.base=zp; self.z=ZipFile(zp); self.names=[n for n in self.z.namelist() if not n.endswith('/')]
        else: raise FileNotFoundError(f'No {subject}/ or {subject}.zip')
    def close(self):
        if self.z:self.z.close()
    def find(self,run,kind):
        n=[x for x in self.names if run_token(Path(x).name,self.subject)==run]
        if kind=='WINDOWS': m=[x for x in n if x.upper().endswith('_WINDOWS.CSV')]
        elif kind=='RAW': m=[x for x in n if x.upper().endswith('_RAW.CSV')]
        elif kind=='SP20': m=[x for x in n if Path(x).name.upper().startswith('SP20') and x.upper().endswith('.CSV')]
        else: raise ValueError(kind)
        m=sorted(m,key=lambda x:(1 if re.search(r'\(\d+\)',Path(x).name) else 0,len(x)))
        if not m: raise FileNotFoundError(f'{self.subject}/{run}: missing {kind}')
        return m[0]
    def csv(self,m,**kw):
        if self.kind=='dir': return pd.read_csv(self.base/m,**kw)
        with self.z.open(m) as f:return pd.read_csv(BytesIO(f.read()),**kw)

def num(df,c,required=True):
    if c not in df:
        if required: raise ValueError(f'Missing {c}')
        return np.full(len(df),np.nan)
    return pd.to_numeric(df[c],errors='coerce').to_numpy(float)
def boolcol(df,c,default=False):
    if c not in df:return np.full(len(df),default,bool)
    return df[c].astype(str).str.strip().str.lower().isin(('true','1','yes','valid','accepted')).to_numpy(bool)
def med(x):
    x=np.sort(np.asarray(x,float)[np.isfinite(x)]); n=len(x)
    if not n:return np.nan
    return float(x[n//2] if n%2 else (x[n//2-1]+x[n//2])/2)
def quant(x,q):
    x=np.sort(np.asarray(x,float)[np.isfinite(x)]); n=len(x)
    if not n:return np.nan
    p=(n-1)*q;i=int(np.floor(p));f=p-i;j=min(i+1,n-1);return float(x[i]+(x[j]-x[i])*f)
def rsig(x):
    x=np.asarray(x,float);x=x[np.isfinite(x)]
    if not len(x):return np.nan
    m=med(x);return 1.4826*med(np.abs(x-m))
def centered_std(x,w):
    x=np.asarray(x,float);n=len(x);r=w//2;ok=np.isfinite(x);xv=np.where(ok,x,0.)
    cc=np.r_[0,np.cumsum(ok.astype(int))];cs=np.r_[0.,np.cumsum(xv)];cq=np.r_[0.,np.cumsum(xv*xv)]
    idx=np.arange(n);a=np.maximum(0,idx-r);b=np.minimum(n,idx+r+1);cnt=cc[b]-cc[a];sm=cs[b]-cs[a];sq=cq[b]-cq[a]
    out=np.full(n,np.nan);g=cnt>0;mu=np.zeros(n);mu[g]=sm[g]/cnt[g];v=np.zeros(n);v[g]=sq[g]/cnt[g]-mu[g]**2;out[g]=np.sqrt(np.maximum(0,v[g]));return out

def centered_rms(x,w):
    x=np.asarray(x,float);n=len(x);r=w//2;ok=np.isfinite(x);xv=np.where(ok,x,0.)
    cc=np.r_[0,np.cumsum(ok.astype(int))];cq=np.r_[0.,np.cumsum(xv*xv)]
    idx=np.arange(n);a=np.maximum(0,idx-r);b=np.minimum(n,idx+r+1);cnt=cc[b]-cc[a];sq=cq[b]-cq[a]
    out=np.full(n,np.nan);g=cnt>0;out[g]=np.sqrt(sq[g]/cnt[g]);return out

def imu_activity(raw,ch,prefix=None):
    r=raw if prefix is None else raw.iloc[:prefix]
    ax,ay,az=[num(r,f'ch{ch}_acc_{a}_g') for a in 'xyz'];gx,gy,gz=[num(r,f'ch{ch}_gyro_{a}_dps') for a in 'xyz']
    amag=np.sqrt(ax*ax+ay*ay+az*az);bx,by,bz=med(gx),med(gy),med(gz);gmag=np.sqrt((gx-bx)**2+(gy-by)**2+(gz-bz)**2);w=max(5,int(round(IMU_SUPPORT_S*FS)))
    return centered_std(amag,w),centered_rms(gmag,w)

def adaptive_limit(x,lo,hi):
    x=np.sort(np.asarray(x,float)[np.isfinite(x)])
    if not len(x):return np.nan,lo
    k=max(1,int(np.ceil(BASE_FRAC*len(x))));b=x[:k];m=med(b);s=rsig(b);s=0 if not np.isfinite(s) else s
    return m,float(np.clip(m+BASE_K*s,lo,hi))

def cstats(x,a,b):
    z=np.asarray(x[a:b],float);z=z[np.isfinite(z)]
    if len(z)<max(5,(b-a)//2):return None
    return {'median':med(z),'sigma':rsig(z),'p95':quant(z,.95),'p99':quant(z,.99)}

def initial_limits(raw,ch):
    n=min(len(raw),int(round(WINDOW_S*FS)));acc,gyr=imu_activity(raw,ch,n);ab,aa=adaptive_limit(acc,ACC_FLOOR,ACC_CEIL);gb,ga=adaptive_limit(gyr,GYR_FLOOR,GYR_CEIL)
    settle=int(CAL_SETTLE*FS);seg=int(CAL_SEG*FS);end=min(n,int(CAL_SEARCH*FS));best=None;score=np.inf
    for s in range(settle,max(settle,end-seg+1),int(FS)):
        if s+seg>end:break
        A,G=cstats(acc,s,s+seg),cstats(gyr,s,s+seg)
        if A is None or G is None:continue
        aon=max(A['p99'],A['median']+CAL_ON_K*A['sigma']);gon=max(G['p99'],G['median']+CAL_ON_K*G['sigma'])
        if aon>ACC_MAX or gon>GYR_MAX:continue
        z=aon/ACC_FLOOR+gon/GYR_FLOOR
        if z<score:score=z;best=(A,G)
    if best is None:return {'accel_on':aa,'gyro_on':ga,'accel_off':max(RELEASE_RATIO*aa,ab),'gyro_off':max(RELEASE_RATIO*ga,gb)}
    A,G=best;aon=max(aa,A['p99'],A['median']+CAL_ON_K*A['sigma']);gon=max(ga,G['p99'],G['median']+CAL_ON_K*G['sigma'])
    aoff=min(aon,max(RELEASE_RATIO*aon,A['p95'],A['median']+CAL_OFF_K*A['sigma']));goff=min(gon,max(RELEASE_RATIO*gon,G['p95'],G['median']+CAL_OFF_K*G['sigma']))
    return {'accel_on':float(aon),'gyro_on':float(gon),'accel_off':float(aoff),'gyro_off':float(goff)}

def freeze_fixed(devraw):
    rows=[]
    for s,r,raw in devraw:
        for ch in (1,2):rows.append({'subject':s,'run':r,'channel':f'CH{ch}',**initial_limits(raw,ch)})
    d=pd.DataFrame(rows);fixed={};summary=[]
    for ch in ('CH1','CH2'):
        q=d[d.channel.eq(ch)];vals={k:float(np.nanmedian(q[k])) for k in ('accel_on','gyro_on','accel_off','gyro_off')};fixed[ch]=vals;summary.append({'channel':ch,'development_runs':len(q),**vals})
    return fixed,d,pd.DataFrame(summary)

def fixed_mask(acc,gyr,lim):
    n=len(acc);above=(acc>lim['accel_on'])|(gyr>lim['gyro_on']);below=(acc<lim['accel_off'])&(gyr<lim['gyro_off']);on=max(1,int(round(ONSET_S*FS)));off=max(1,int(round(RELEASE_S*FS)))
    m=np.zeros(n,bool);active=False;oc=rc=0
    for i in range(n):
        if not active:
            oc=oc+1 if above[i] else 0
            if oc>=on:active=True;m[max(0,i-on+1):i+1]=True;rc=0
        else:
            m[i]=True;rc=rc+1 if below[i] else 0
            if rc>=off:active=False;oc=rc=0
    pad=max(0,int(round(PAD_S*FS)))
    if pad and m.any():
        idx=np.flatnonzero(m);starts=idx[np.r_[True,np.diff(idx)>1]];ends=idx[np.r_[np.diff(idx)>1,True]];diff=np.zeros(n+1,int)
        for a,b in zip(starts,ends):a=max(0,a-pad);b=min(n-1,b+pad);diff[a]+=1;diff[b+1]-=1
        m=np.cumsum(diff[:-1])>0
    return m

def add_fixed_burden(raw,w,fixed):
    t=num(raw,'common_time_us')/1e6;out=w.copy()
    for ch in (1,2):
        acc,gyr=imu_activity(raw,ch);mask=fixed_mask(acc,gyr,fixed[f'CH{ch}']);b=[]
        for r in w.itertuples(index=False):
            sel=(t>=float(r.start_s))&(t<float(r.end_s));b.append(100*np.mean(mask[sel]) if np.any(sel) else np.nan)
        out[f'ch{ch}_fixed_motion_percent']=b
    return out
# =============================================================================
# WINDOW PREPARATION + EIGHT ABLATION METHODS
# =============================================================================

def prep_windows(w):
    need=('start_s','end_s','ch1_hr_bpm','ch2_hr_bpm','ch1_spo2_percent','ch2_spo2_percent','ch1_sqi','ch2_sqi','ch1_motion_percent','ch2_motion_percent')
    miss=set(need)-set(w.columns)
    if miss:raise ValueError(f'WINDOWS missing {sorted(miss)}')
    w=w.sort_values('end_s',kind='stable').copy()
    for c in need:w[c]=pd.to_numeric(w[c],errors='coerce')
    for ch in (1,2):
        ac=f'ch{ch}_accept_threshold';w[ac]=pd.to_numeric(w[ac],errors='coerce') if ac in w else SQI_FALLBACK
        clipped=f'ch{ch}_clipped';stale=f'ch{ch}_stale_fraction'
        w[f'ch{ch}_base_not_clipped']=~boolcol(w,clipped,False) if clipped in w else True
        if stale in w:
            z=pd.to_numeric(w[stale],errors='coerce');w[f'ch{ch}_base_fresh']=~(z>=.90)
        else:w[f'ch{ch}_base_fresh']=True
    return w

def average(x1,x2,e1,e2):
    o=np.full(len(x1),np.nan);o[e1&~e2]=x1[e1&~e2];o[e2&~e1]=x2[e2&~e1];b=e1&e2;o[b]=(x1[b]+x2[b])/2;return o

def fuse(x1,x2,e1,e2,w1,w2,consistency=None):
    o=np.full(len(x1),np.nan);o[e1&~e2]=x1[e1&~e2];o[e2&~e1]=x2[e2&~e1]
    for i in np.flatnonzero(e1&e2):
        if consistency is not None and abs(x1[i]-x2[i])>consistency:
            o[i]=x1[i] if w1[i]>w2[i] else (x2[i] if w2[i]>w1[i] else (x1[i]+x2[i])/2)
        else:
            den=w1[i]+w2[i];o[i]=(w1[i]*x1[i]+w2[i]*x2[i])/den if den>0 else (x1[i]+x2[i])/2
    return o

def predictions(w):
    out=pd.DataFrame({'hemox_t':pd.to_numeric(w.end_s,errors='coerce')})
    for var in ('hr','spo2'):
        c1=f'ch1_hr_bpm' if var=='hr' else 'ch1_spo2_percent';c2=f'ch2_hr_bpm' if var=='hr' else 'ch2_spo2_percent'
        x1=pd.to_numeric(w[c1],errors='coerce').to_numpy(float);x2=pd.to_numeric(w[c2],errors='coerce').to_numpy(float)
        q1,q2=[pd.to_numeric(w[f'ch{ch}_sqi'],errors='coerce').to_numpy(float) for ch in (1,2)]
        t1,t2=[pd.to_numeric(w[f'ch{ch}_accept_threshold'],errors='coerce').to_numpy(float) for ch in (1,2)]
        ma1,ma2=[pd.to_numeric(w[f'ch{ch}_motion_percent'],errors='coerce').to_numpy(float) for ch in (1,2)]
        mf1,mf2=[pd.to_numeric(w[f'ch{ch}_fixed_motion_percent'],errors='coerce').to_numpy(float) for ch in (1,2)]
        b1=np.isfinite(x1)&w.ch1_base_not_clipped.to_numpy(bool)&w.ch1_base_fresh.to_numpy(bool);b2=np.isfinite(x2)&w.ch2_base_not_clipped.to_numpy(bool)&w.ch2_base_fresh.to_numpy(bool)
        s1=b1&np.isfinite(q1)&np.isfinite(t1)&(q1>=t1);s2=b2&np.isfinite(q2)&np.isfinite(t2)&(q2>=t2)
        am1=b1&np.isfinite(ma1)&(ma1<=MOTION_REJECT);am2=b2&np.isfinite(ma2)&(ma2<=MOTION_REJECT)
        fm1=b1&np.isfinite(mf1)&(mf1<=MOTION_REJECT);fm2=b2&np.isfinite(mf2)&(mf2<=MOTION_REJECT)
        qw1=np.maximum(0,np.nan_to_num(q1,nan=0));qw2=np.maximum(0,np.nan_to_num(q2,nan=0))
        aw1=np.maximum(0,1-np.nan_to_num(ma1,nan=100)/100);aw2=np.maximum(0,1-np.nan_to_num(ma2,nan=100)/100)
        fw1=np.maximum(0,1-np.nan_to_num(mf1,nan=100)/100);fw2=np.maximum(0,1-np.nan_to_num(mf2,nan=100)/100)
        con=HR_CONS if var=='hr' else SPO2_CONS
        P={
            'CH1':np.where(b1,x1,np.nan),
            'CH2':np.where(b2,x2,np.nan),
            'Simple average':average(x1,x2,b1,b2),
            'SQI-only fusion':fuse(x1,x2,s1,s2,qw1,qw2),
            'Motion-only fusion':fuse(x1,x2,am1,am2,aw1,aw2),
            'Fixed thresholds':fuse(x1,x2,s1&fm1,s2&fm2,qw1*fw1,qw2*fw2),
            'Adaptive thresholds':fuse(x1,x2,s1&am1,s2&am2,qw1*aw1,qw2*aw2),
            'Complete HEMOX':fuse(x1,x2,s1&am1,s2&am2,qw1*aw1,qw2*aw2,con),
        }
        for method,val in P.items():out[f'{var}__{safe(method)}']=val
    return out

# =============================================================================
# REFERENCE READING / PAIRING
# =============================================================================

def read_sp20(src,member):
    r=src.csv(member,skiprows=1,skipinitialspace=True,dtype=str);r.columns=[str(c).strip().lstrip('#').strip() for c in r.columns];need=('Time','SpO2(%)','PR(bpm)')
    if not set(need).issubset(r.columns):
        r=src.csv(member,dtype=str);r.columns=[str(c).strip().lstrip('#').strip() for c in r.columns]
    if not set(need).issubset(r.columns):raise ValueError(f'SP20 missing {set(need)-set(r.columns)}')
    raw=r.Time.astype(str).str.strip();tm=pd.to_datetime(raw,errors='coerce',format='%m/%d/%Y %H:%M:%S')
    if tm.isna().any():
        try:tm=tm.fillna(pd.to_datetime(raw,errors='coerce',format='mixed'))
        except Exception:pass
    if tm.isna().any() or not tm.is_monotonic_increasing:raise ValueError('Invalid SP20 timestamps')
    return pd.DataFrame({'ref_t':(tm-tm.iloc[0]).dt.total_seconds(),'sp20_timestamp':tm.astype(str),'ref_spo2':pd.to_numeric(r['SpO2(%)'],errors='coerce'),'ref_hr':pd.to_numeric(r['PR(bpm)'],errors='coerce')})

def load_offsets(path):
    if path is None:return {}
    d=pd.read_csv(path);req={'subject','run','offset_seconds','method','evidence'}
    if not req.issubset(d.columns):raise ValueError(f'offsets CSV missing {req-set(d.columns)}')
    return {(str(r.subject).upper(),str(r.run).upper()):{'offset':float(r.offset_seconds),'method':str(r.method),'evidence':str(r.evidence)} for r in d.itertuples(index=False)}

def pair(sp,he,offset,tol,duration):
    sp=sp[(sp.ref_t>=0)&(sp.ref_t<=duration)].copy();he=he[(he.hemox_t>=0)&(he.hemox_t<=duration)].copy();he['aligned']=he.hemox_t+offset
    R=list(sp.itertuples(index=False));H=list(he.itertuples(index=False));used=set();rows=[]
    if not R or not H:return pd.DataFrame()
    lo=max(R[0].ref_t,H[0].aligned);hi=min(R[-1].ref_t,H[-1].aligned)
    for sr in R:
        if not lo<=sr.ref_t<=hi:continue
        bi=None;bd=np.inf
        for i,hr in enumerate(H):
            if i in used or not lo<=hr.aligned<=hi:continue
            dt=abs(hr.aligned-sr.ref_t)
            if dt<bd:bi,bd=i,dt
            if hr.aligned>sr.ref_t+tol and bi is not None:break
        if bi is None or bd>tol:continue
        used.add(bi);hr=H[bi];row={'sp20_elapsed_s':sr.ref_t,'sp20_timestamp':sr.sp20_timestamp,'hemox_elapsed_s':hr.hemox_t,'delta_t_s':hr.aligned-sr.ref_t,'sp20_spo2':sr.ref_spo2,'sp20_hr':sr.ref_hr}
        for c in he.columns:
            if c not in ('hemox_t','aligned'):row[c]=getattr(hr,c)
        rows.append(row)
    return pd.DataFrame(rows)
# =============================================================================
# BUILD DATASET
# =============================================================================

def build(root,duration,tol,off):
    src={}
    for s in PARTICIPANTS:
        try:src[s]=Source(root,s)
        except Exception as e:warnings.warn(f'{s}: {e}')
    devraw=[];cache={}
    for s in DEV:
        if s not in src:continue
        for run in RUNS:
            try:
                m=src[s].find(run,'RAW');r=src[s].csv(m,low_memory=False);cache[(s,run)]=r;devraw.append((s,run,r))
            except Exception as e:warnings.warn(f'DEV {s}/{run}: {e}')
    if not devraw:raise RuntimeError('No P01-P02 RAW data; cannot freeze fixed thresholds')
    fixed,fd,fs=freeze_fixed(devraw)
    paired=[];audit=[]
    for s in PARTICIPANTS:
        if s not in src:continue
        for run in RUNS:
            a={'subject':s,'run':run,'status':'SKIPPED','reason':'','paired_rows':0}
            try:
                wm=src[s].find(run,'WINDOWS');rm=src[s].find(run,'RAW');sm=src[s].find(run,'SP20')
                w=prep_windows(src[s].csv(wm,low_memory=False));r=cache.get((s,run));r=src[s].csv(rm,low_memory=False) if r is None else r;sp=read_sp20(src[s],sm)
                w=add_fixed_burden(r,w,fixed);he=predictions(w)
                oi=off.get((s,run),{'offset':0.,'method':'provisional_zero_offset','evidence':''});q=pair(sp,he,oi['offset'],tol,duration)
                if q.empty:raise RuntimeError('No matched HEMOX-SP20 pairs')
                q.insert(0,'run',run);q.insert(0,'subject',s);q['cohort']='Development' if s in DEV else 'Held-out';q['run_type']='Stationary' if run in ('S01','S02') else 'Motion';q['alignment_method']=oi['method'];q['alignment_offset_s']=oi['offset'];paired.append(q)
                a.update(status='OK',paired_rows=len(q),windows_rows=len(w),raw_rows=len(r),sp20_rows=len(sp),windows_file=wm,raw_file=rm,sp20_file=sm);print(f'[OK] {s}/{run}: windows={len(w)}, pairs={len(q)}')
            except Exception as e:a['reason']=str(e);print(f'[SKIP] {s}/{run}: {e}')
            audit.append(a)
    for x in src.values():x.close()
    if not paired:raise RuntimeError('No usable recordings')
    return pd.concat(paired,ignore_index=True),pd.DataFrame(audit),fd,fs

# =============================================================================
# CALIBRATION / METRICS / BOOTSTRAP
# =============================================================================

def fit_m1(p,do=True):
    d=p[p.subject.isin(DEV)];rows=[];out={}
    for m in METHODS:
        x=pd.to_numeric(d[f'spo2__{safe(m)}'],errors='coerce').to_numpy(float);y=pd.to_numeric(d.sp20_spo2,errors='coerce').to_numpy(float);ok=np.isfinite(x)&np.isfinite(y);off=float(np.mean(x[ok]-y[ok])) if do and ok.any() else 0.0;out[m]=off;rows.append({'method':m,'development_pairs':int(ok.sum()),'offset_pp_subtracted':off,'fit':'P01-P02 M1' if do else 'none'})
    return out,pd.DataFrame(rows)

def apply_m1(p,offs):
    p=p.copy()
    for m in METHODS:
        s=safe(m);p[f'spo2_cal__{s}']=np.clip(pd.to_numeric(p[f'spo2__{s}'],errors='coerce')-offs[m],0,100);p[f'hr_cal__{s}']=pd.to_numeric(p[f'hr__{s}'],errors='coerce')
    return p

def metric(pred,ref,thr):
    a=np.asarray(pred,float);b=np.asarray(ref,float);ok=np.isfinite(a)&np.isfinite(b);a=a[ok];b=b[ok]
    if not len(a):return {'valid_pairs':0,'bias':np.nan,'absolute_bias':np.nan,'sd_error':np.nan,'loa_low':np.nan,'loa_high':np.nan,'mae':np.nan,'rmse':np.nan,'within_threshold_pct':np.nan}
    e=a-b;bias=float(e.mean());sd=float(e.std(ddof=1)) if len(e)>1 else np.nan
    return {'valid_pairs':int(len(e)),'bias':bias,'absolute_bias':abs(bias),'sd_error':sd,'loa_low':bias-1.96*sd if np.isfinite(sd) else np.nan,'loa_high':bias+1.96*sd if np.isfinite(sd) else np.nan,'mae':float(np.mean(np.abs(e))),'rmse':float(np.sqrt(np.mean(e*e))),'within_threshold_pct':100*float(np.mean(np.abs(e)<=thr))}

def _subject_sse_counts(df, pred, ref):
    """Return subject-level SSE and valid-pair counts for exact cluster bootstrap."""
    rows = []
    for subject, g in df.groupby('subject', sort=True):
        a = pd.to_numeric(g[pred], errors='coerce').to_numpy(float)
        b = pd.to_numeric(g[ref], errors='coerce').to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.any():
            e = a[ok] - b[ok]
            rows.append((subject, float(np.sum(e * e)), int(ok.sum())))
        else:
            rows.append((subject, 0.0, 0))
    return rows


def boot_rmse(df,pred,ref,thr,n,rng):
    """
    Exact participant-cluster bootstrap for RMSE, implemented from subject-level
    SSE/count sufficient statistics. This is mathematically equivalent to
    rebuilding the resampled dataframe, but is orders of magnitude faster.
    """
    stats = _subject_sse_counts(df, pred, ref)
    if len(stats) < 2 or n <= 0:
        return np.nan, np.nan

    sse = np.asarray([z[1] for z in stats], float)
    cnt = np.asarray([z[2] for z in stats], float)
    k = len(stats)

    # Each row contains multiplicities of the k participants in one bootstrap.
    draws = rng.multinomial(k, np.full(k, 1.0 / k), size=n)
    total_sse = draws @ sse
    total_cnt = draws @ cnt

    good = total_cnt > 0
    values = np.full(n, np.nan, float)
    values[good] = np.sqrt(total_sse[good] / total_cnt[good])
    values = values[np.isfinite(values)]

    return tuple(np.quantile(values, [.025, .975])) if len(values) else (np.nan, np.nan)

def evaluate(p,cohort,nboot,seed):
    print(f'[BOOT] {cohort}: participant-cluster RMSE CIs, B={nboot}')
    d=p[p.cohort.eq(cohort)].copy();rng=np.random.default_rng(seed);rows=[]
    strata=(('All',d),('Stationary',d[d.run_type.eq('Stationary')]),('Motion',d[d.run_type.eq('Motion')]))
    for st,q in strata:
        if q.empty:continue
        for var,ref,prefix,thr in (('SpO2','sp20_spo2','spo2_cal__',SPO2_WITHIN),('HR','sp20_hr','hr_cal__',HR_WITHIN)):
            total=int(np.isfinite(pd.to_numeric(q[ref],errors='coerce')).sum())
            for m in METHODS:
                c=prefix+safe(m);z=metric(q[c],q[ref],thr);lo,hi=boot_rmse(q,c,ref,thr,nboot,rng);rows.append({'cohort':cohort,'stratum':st,'variable':var,'method':m,'participants':int(q.subject.nunique()),'runs':int(q[['subject','run']].drop_duplicates().shape[0]),'reference_pairs':total,'availability_pct':100*z['valid_pairs']/total if total else np.nan,'rmse_ci95_low':lo,'rmse_ci95_high':hi,**z})
    return pd.DataFrame(rows)

def pairwise(p,nboot,seed):
    """
    Common-sample ablation comparison versus Complete HEMOX.
    Participant-cluster bootstrap is computed from per-subject SSE sufficient
    statistics, avoiding repeated pandas concatenation.
    """
    d=p[p.cohort.eq('Held-out')].copy()
    rng=np.random.default_rng(seed)
    rows=[]

    strata=(
        ('All',d),
        ('Stationary',d[d.run_type.eq('Stationary')]),
        ('Motion',d[d.run_type.eq('Motion')]),
    )

    for st,q in strata:
        print(f'[BOOT] Pairwise common-sample comparison: {st}')
        for var,ref,prefix,thr in (
            ('SpO2','sp20_spo2','spo2_cal__',SPO2_WITHIN),
            ('HR','sp20_hr','hr_cal__',HR_WITHIN),
        ):
            cc=prefix+safe('Complete HEMOX')

            for m in METHODS[:-1]:
                mc=prefix+safe(m)
                mask=(
                    np.isfinite(pd.to_numeric(q[mc],errors='coerce'))
                    & np.isfinite(pd.to_numeric(q[cc],errors='coerce'))
                    & np.isfinite(pd.to_numeric(q[ref],errors='coerce'))
                )
                g=q.loc[mask].copy()
                if g.empty:
                    continue

                M=metric(g[mc],g[ref],thr)
                C=metric(g[cc],g[ref],thr)

                subs=sorted(g.subject.unique())
                k=len(subs)

                method_sse=[]
                complete_sse=[]
                counts=[]

                for s in subs:
                    gs=g[g.subject.eq(s)]
                    a=pd.to_numeric(gs[mc],errors='coerce').to_numpy(float)
                    c=pd.to_numeric(gs[cc],errors='coerce').to_numpy(float)
                    r=pd.to_numeric(gs[ref],errors='coerce').to_numpy(float)
                    ok=np.isfinite(a)&np.isfinite(c)&np.isfinite(r)

                    if ok.any():
                        em=a[ok]-r[ok]
                        ec=c[ok]-r[ok]
                        method_sse.append(float(np.sum(em*em)))
                        complete_sse.append(float(np.sum(ec*ec)))
                        counts.append(int(ok.sum()))
                    else:
                        method_sse.append(0.0)
                        complete_sse.append(0.0)
                        counts.append(0)

                method_sse=np.asarray(method_sse,float)
                complete_sse=np.asarray(complete_sse,float)
                counts=np.asarray(counts,float)

                if k>=2 and nboot>0:
                    draws=rng.multinomial(
                        k,
                        np.full(k,1.0/k),
                        size=nboot,
                    )
                    nobs=draws@counts
                    good=nobs>0
                    rm_m=np.full(nboot,np.nan,float)
                    rm_c=np.full(nboot,np.nan,float)
                    rm_m[good]=np.sqrt((draws[good]@method_sse)/nobs[good])
                    rm_c[good]=np.sqrt((draws[good]@complete_sse)/nobs[good])
                    z=(rm_m-rm_c)
                    z=z[np.isfinite(z)]
                    lo,hi=(
                        np.quantile(z,[.025,.975])
                        if len(z)
                        else (np.nan,np.nan)
                    )
                else:
                    lo,hi=np.nan,np.nan

                rows.append({
                    'stratum':st,
                    'variable':var,
                    'method':m,
                    'reference_method':'Complete HEMOX',
                    'common_pairs':len(g),
                    'method_rmse':M['rmse'],
                    'complete_rmse':C['rmse'],
                    'delta_rmse_method_minus_complete':M['rmse']-C['rmse'],
                    'delta_rmse_ci95_low':lo,
                    'delta_rmse_ci95_high':hi,
                })

    return pd.DataFrame(rows)
# =============================================================================
# FIGURES / TABLE
# =============================================================================

def plot_rmse(met,out):
    q=met[(met.cohort.eq('Held-out'))&(met.stratum.eq('All'))]
    labels=['CH1','CH2','Simple\navg.','SQI-only','Motion-only','Fixed\nthr.','Adaptive\nthr.','Complete\nHEMOX']
    colors=[GREY,GREY,GREY,BLUE,ORANGE,RED,BLUE,GREEN];x=np.arange(len(METHODS));fig,ax=plt.subplots(2,1,figsize=(7.2,4.8),sharex=True)
    for a,var,yl in ((ax[0],'SpO2',r'SpO$_2$ RMSE (pp)'),(ax[1],'HR','HR RMSE (bpm)')):
        d=q[q.variable.eq(var)].set_index('method').reindex(METHODS);y=d.rmse.to_numpy(float);lo=d.rmse_ci95_low.to_numpy(float);hi=d.rmse_ci95_high.to_numpy(float);err=np.vstack([np.maximum(0,y-lo),np.maximum(0,hi-y)])
        a.bar(x,y,.68,color=colors,edgecolor='none',zorder=3);a.errorbar(x,y,yerr=err,fmt='none',ecolor='#30343B',capsize=2.5,lw=.9,zorder=4);a.set_ylabel(yl);a.grid(axis='y',color=GRID,lw=.6);a.set_axisbelow(True)
    ax[0].text(.01,1.02,'(a)',transform=ax[0].transAxes);ax[1].text(.01,1.02,'(b)',transform=ax[1].transAxes);ax[1].set_xticks(x);ax[1].set_xticklabels(labels);ax[1].set_xlabel('Ablation method');fig.subplots_adjust(left=.12,right=.99,top=.98,bottom=.22,hspace=.30)
    fig.savefig(out/'Fig_Ablation_RMSE.png',dpi=DPI,bbox_inches='tight');fig.savefig(out/'Fig_Ablation_RMSE.svg',bbox_inches='tight');plt.close(fig)

def plot_availability(met,out):
    q=met[(met.cohort.eq('Held-out'))&(met.stratum.eq('All'))]
    labels=['CH1','CH2','Simple\navg.','SQI-only','Motion-only','Fixed\nthr.','Adaptive\nthr.','Complete\nHEMOX']
    colors=[GREY,GREY,GREY,BLUE,ORANGE,RED,BLUE,GREEN];x=np.arange(len(METHODS));fig,ax=plt.subplots(2,1,figsize=(7.2,4.8),sharex=True)
    for a,var,yl in ((ax[0],'SpO2',r'SpO$_2$ availability (%)'),(ax[1],'HR','HR availability (%)')):
        d=q[q.variable.eq(var)].set_index('method').reindex(METHODS);a.bar(x,d.availability_pct,.68,color=colors,edgecolor='none',zorder=3);a.set_ylim(0,105);a.set_ylabel(yl);a.grid(axis='y',color=GRID,lw=.6);a.set_axisbelow(True)
    ax[0].text(.01,1.02,'(a)',transform=ax[0].transAxes);ax[1].text(.01,1.02,'(b)',transform=ax[1].transAxes);ax[1].set_xticks(x);ax[1].set_xticklabels(labels);ax[1].set_xlabel('Ablation method');fig.subplots_adjust(left=.12,right=.99,top=.98,bottom=.22,hspace=.30)
    fig.savefig(out/'Fig_Ablation_Availability.png',dpi=DPI,bbox_inches='tight');fig.savefig(out/'Fig_Ablation_Availability.svg',bbox_inches='tight');plt.close(fig)


def ch1_ch2_complete_summary(pw):
    """
    Reviewer-facing direct comparison.

    Uses only Held-out / All pairwise-common rows already calculated by
    pairwise(). Positive RMSE reduction means Complete HEMOX has lower RMSE.
    """
    q = pw[
        pw.stratum.eq('All')
        & pw.method.isin(('CH1','CH2'))
        & pw.variable.isin(('SpO2','HR'))
    ].copy()

    if len(q) != 4:
        raise RuntimeError(
            f'Expected four CH1/CH2 pairwise rows, found {len(q)}.'
        )

    method_rmse = pd.to_numeric(q.method_rmse, errors='coerce')
    complete_rmse = pd.to_numeric(q.complete_rmse, errors='coerce')

    q['rmse_reduction_pct'] = np.where(
        method_rmse > 0,
        100.0 * (method_rmse - complete_rmse) / method_rmse,
        np.nan,
    )
    q['rmse_ratio_baseline_over_complete'] = np.where(
        complete_rmse > 0,
        method_rmse / complete_rmse,
        np.nan,
    )
    q['ci_excludes_zero'] = (
        (pd.to_numeric(q.delta_rmse_ci95_low, errors='coerce') > 0)
        | (pd.to_numeric(q.delta_rmse_ci95_high, errors='coerce') < 0)
    )

    cols = [
        'variable','method','reference_method','common_pairs',
        'method_rmse','complete_rmse',
        'delta_rmse_method_minus_complete',
        'rmse_reduction_pct',
        'rmse_ratio_baseline_over_complete',
        'delta_rmse_ci95_low','delta_rmse_ci95_high',
        'ci_excludes_zero',
    ]
    return q[cols].sort_values(['variable','method']).reset_index(drop=True)


def plot_ch1_ch2_pairwise(summary, out):
    """
    Correct publication figure for direct improvement claims.

    Each x-position is one matched comparison:
      CH1 vs Complete HEMOX
      CH2 vs Complete HEMOX

    Complete HEMOX can have a slightly different RMSE in the two SpO2
    comparisons because the pairwise-common observation sets can differ.
    """
    methods = ('CH1','CH2')
    x = np.arange(2, dtype=float)
    width = 0.30

    fig, axes = plt.subplots(
        2, 1,
        figsize=(7.2, 4.8),
        sharex=True,
    )

    for ax, var, ylabel, panel in (
        (axes[0], 'SpO2', r'SpO$_2$ RMSE (pp)', '(a)'),
        (axes[1], 'HR', 'HR RMSE (bpm)', '(b)'),
    ):
        q = (
            summary[summary.variable.eq(var)]
            .set_index('method')
            .reindex(methods)
        )

        base = q.method_rmse.to_numpy(float)
        comp = q.complete_rmse.to_numpy(float)
        red = q.rmse_reduction_pct.to_numpy(float)
        ratio = q.rmse_ratio_baseline_over_complete.to_numpy(float)

        b1 = ax.bar(
            x - width/2, base, width,
            color=GREY, edgecolor='none', zorder=3,
        )
        b2 = ax.bar(
            x + width/2, comp, width,
            color=GREEN, edgecolor='none', zorder=3,
        )

        ymax = float(np.nanmax(np.r_[base, comp]))

        for bars, vals in ((b1, base), (b2, comp)):
            for bar, val in zip(bars, vals):
                ax.text(
                    bar.get_x()+bar.get_width()/2,
                    val + 0.025*ymax,
                    f'{val:.3f}',
                    ha='center', va='bottom',
                    fontsize=8.0, color=BLACK,
                    fontweight='bold', zorder=5,
                )

        ax.set_ylim(0, ymax*1.62)

        for i, method in enumerate(methods):
            if not np.isfinite(red[i]):
                headline = 'Unavailable'
            elif abs(red[i]) < 0.1:
                headline = 'Comparable RMSE'
            elif red[i] > 0:
                headline = f'{red[i]:.1f}% reduction'
            else:
                headline = f'{abs(red[i]):.1f}% increase'

            color = RED if method == 'CH1' else BLUE
            ax.text(
                x[i], 0.965,
                f'{headline}\n'
                f'RMSE ratio={ratio[i]:.2f}x',
                transform=ax.get_xaxis_transform(),
                ha='center', va='top',
                fontsize=7.4, color=color,
                bbox=dict(
                    boxstyle='round,pad=0.25',
                    facecolor='white',
                    edgecolor=color,
                    linewidth=0.85,
                ),
                zorder=10,
            )

        ax.set_ylabel(ylabel)
        ax.grid(axis='y', color=GRID, lw=.6, zorder=0)
        ax.set_axisbelow(True)
        ax.text(
            .01, 1.02, panel,
            transform=ax.transAxes,
            fontsize=10.2,
            fontweight='bold',
        )

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(
        ['CH1 vs Complete HEMOX', 'CH2 vs Complete HEMOX']
    )
    axes[1].set_xlabel('Window-level pairwise comparison')

    handles = [
        plt.Rectangle((0,0),1,1,color=GREY),
        plt.Rectangle((0,0),1,1,color=GREEN),
    ]
    fig.legend(
        handles,
        ['Individual channel','Complete HEMOX'],
        loc='upper center',
        bbox_to_anchor=(.5,.995),
        ncol=2,
        frameon=False,
        fontsize=8,
    )

    fig.text(
        0.5, 0.025,
        'Primary agreement: N=2,202 (SpO2) and N=2,255 (HR), '
        'using native quality-valid FUSED_TRENDS.',
        ha='center', va='bottom',
        fontsize=7.0, color=GREY, fontweight='bold'
    )

    fig.subplots_adjust(
        left=.13, right=.985, top=.91, bottom=.20, hspace=.32
    )

    fig.savefig(
        out/'Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common.png',
        dpi=DPI, bbox_inches='tight',
    )
    fig.savefig(
        out/'Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common.svg',
        bbox_inches='tight',
    )
    plt.close(fig)


def write_population_note(out, audit):
    """
    Document the denominator distinction so the output folder is self-auditing.
    """
    held = audit[
        audit.subject.isin(HELDOUT)
        & audit.status.eq('OK')
    ].copy()

    methods = sorted(set(
        str(x) for x in held.get('alignment_method', pd.Series(dtype=str))
        if str(x).strip()
    ))

    lines = [
        'HEMOX R2 analysis-population note',
        '=================================',
        '',
        'Primary agreement analysis:',
        '  native quality-valid FUSED_TRENDS',
        '  SpO2 N = 2,202; HR N = 2,255',
        '',
        'This ablation analysis:',
        '  WINDOW/PARAMETER_TRENDS-level reconstruction and SP-20 pairing',
        '  method-specific availability is reported separately',
        '  direct method-vs-Complete claims use pairwise-common observations',
        '',
        'Do not substitute the ablation denominators for the primary',
        'FUSED_TRENDS Bland-Altman denominators, and do not force a post-hoc',
        'one-to-one mapping between the two exported representations.',
        '',
        'Alignment labels present in this run:',
        '  ' + (', '.join(methods) if methods else 'none recorded'),
        '',
    ]
    (out/'analysis_population_note.txt').write_text(
        '\n'.join(lines),
        encoding='utf-8',
    )


def latex_table(met, path):
    """
    Write the held-out ablation table directly as LaTeX.
    This avoids pandas.DataFrame.to_latex(), which requires the optional
    Jinja2 dependency in recent pandas versions.
    """
    q = met[
        (met.cohort.eq('Held-out'))
        & (met.stratum.eq('All'))
    ]

    s = (
        q[q.variable.eq('SpO2')]
        .set_index('method')
        .reindex(METHODS)
    )
    h = (
        q[q.variable.eq('HR')]
        .set_index('method')
        .reindex(METHODS)
    )

    def fmt(v):
        try:
            v = float(v)
        except Exception:
            return '--'
        return f'{v:.3f}' if np.isfinite(v) else '--'

    def esc(s):
        return (
            str(s)
            .replace('\\', r'\textbackslash{}')
            .replace('&', r'\&')
            .replace('%', r'\%')
            .replace('_', r'\_')
            .replace('#', r'\#')
        )

    lines = [
        r'\begin{tabular}{lrrrrrr}',
        r'\toprule',
        r'Method & SpO$_2$ RMSE (pp) & SpO$_2$ MAE (pp) & SpO$_2$ avail. (\%) & '
        r'HR RMSE (bpm) & HR MAE (bpm) & HR avail. (\%) \\',
        r'\midrule',
    ]

    for method in METHODS:
        lines.append(
            f"{esc(method)} & "
            f"{fmt(s.loc[method, 'rmse'])} & "
            f"{fmt(s.loc[method, 'mae'])} & "
            f"{fmt(s.loc[method, 'availability_pct'])} & "
            f"{fmt(h.loc[method, 'rmse'])} & "
            f"{fmt(h.loc[method, 'mae'])} & "
            f"{fmt(h.loc[method, 'availability_pct'])} \\\\"
        )

    lines += [
        r'\bottomrule',
        r'\end{tabular}',
        '',
    ]

    path.write_text('\n'.join(lines), encoding='utf-8')


# =============================================================================
# MAIN
# =============================================================================

def main():
    a=cli()
    set_style()

    root=a.data_dir.expanduser().resolve()
    out=a.out_dir.expanduser().resolve()
    out.mkdir(parents=True,exist_ok=True)

    print('='*78)
    print('HEMOX R2: CH1/CH2 vs COMPLETE HEMOX — PAIRWISE-COMMON FIGURE ONLY')
    print('='*78)
    print('P01-P02: development/calibration only')
    print('P03-P10: held-out pairwise-common comparison')
    print('SpO2: method-specific M1 fitted on P01-P02; HR: M0')
    print()

    offsets=load_offsets(a.offsets_csv)
    if not offsets:
        print(
            '[WARNING] No --offsets-csv supplied. '
            'HEMOX-SP20 window pairing will use zero offset and the output '
            'must be interpreted under that alignment assumption.'
        )
    else:
        print(f'[INFO] Loaded {len(offsets)} documented run offsets.')

    # These steps are required to reproduce the original figure correctly.
    P,audit,fd,fs=build(
        root,
        a.duration,
        a.tolerance,
        offsets,
    )

    print(
        '[STEP] Fitting P01-P02 method-specific SpO2 M1 offsets '
        'required by the pairwise-common comparison...'
    )
    m1,cal=fit_m1(P,not a.no_method_specific_m1)
    P=apply_m1(P,m1)

    print('[STEP] Computing held-out pairwise-common CH1/CH2 comparisons...')
    pw=pairwise(P,a.bootstrap,a.seed+2)
    focus=ch1_ch2_complete_summary(pw)

    print('[STEP] Generating Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common only...')
    plot_ch1_ch2_pairwise(focus,out)

    png=out/'Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common.png'
    svg=out/'Fig_CH1_CH2_CompleteHEMOX_Pairwise_Common.svg'

    print()
    print('PAIRWISE-COMMON VALUES USED IN THE FIGURE')
    print(
        focus[
            [
                'variable','method','common_pairs',
                'method_rmse','complete_rmse',
                'rmse_reduction_pct',
                'delta_rmse_ci95_low','delta_rmse_ci95_high',
                'ci_excludes_zero',
            ]
        ].to_string(index=False)
    )
    print()
    print('[DONE] Figure-only analysis completed.')
    print(f'PNG: {png}')
    print(f'SVG: {svg}')
    return png,svg


if __name__=='__main__':
    main()
