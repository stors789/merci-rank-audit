"""Fixed-label FWL effects and explicitly named conditional uncertainties."""
import numpy as np
import pandas as pd
from scipy.stats import t as student_t, false_discovery_control


def projection(x,y):
    u,s,_=np.linalg.svd(x,full_matrices=False)
    rank=int((s>s[0]*max(x.shape)*np.finfo(float).eps).sum())
    q=u[:,:rank]
    return q,y-q@(q.T@y)


def fit(y,x,label,names):
    y=np.asarray(y,float);x=np.asarray(x,float);label=np.asarray(label,float)
    n=len(label);n1=int(label.sum());n0=n-n1
    base=dict(N=n,selected=n1,unselected=n0)
    if min(n1,n0)<5:
        return pd.DataFrame([dict(**base,nuclear_feature=name,status='not_estimable_empty_or_small_group') for name in names])
    q,yr=projection(x,y);zr=label-q@(q.T@label);den=float(zr@zr);df=n-q.shape[1]-1
    if den<1e-10 or df<1:
        return pd.DataFrame([dict(**base,nuclear_feature=name,status='not_estimable_no_residual_label_information') for name in names])
    beta=(zr@yr)/den;error=yr-zr[:,None]*beta
    h=(q*q).sum(axis=1)+zr*zr/den
    if np.max(h)>=1-1e-10:
        return pd.DataFrame([dict(**base,nuclear_feature=name,status='not_estimable_unit_leverage') for name in names])
    se=np.sqrt(((zr/(1-h))[:,None]**2*error**2).sum(axis=0))/den
    classic_se=np.sqrt((error**2).sum(axis=0)/df/den)
    old_cr1_se=np.sqrt((n/df)*(zr[:,None]**2*error**2).sum(axis=0))/den
    sd=y.std(axis=0,ddof=1);label_sd=label.std(ddof=1);limit=student_t.ppf(.975,df)
    rows=[]
    for j,name in enumerate(names):
        valid=bool(sd[j]>1e-14 and se[j]>0 and np.isfinite(se[j]))
        rows.append(dict(**base,nuclear_feature=name,status='estimable' if valid else 'not_estimable_constant_outcome',
            beta=float(beta[j]),beta_per_outcome_SD=float(beta[j]/sd[j]) if sd[j]>0 else np.nan,
            standardized_coefficient=float(beta[j]*label_sd/sd[j]) if sd[j]>0 else np.nan,
            HC3_SE=float(se[j]),CI_low=float(beta[j]-limit*se[j]) if valid else np.nan,
            CI_high=float(beta[j]+limit*se[j]) if valid else np.nan,
            p_value=float(2*student_t.sf(abs(beta[j]/se[j]),df)) if valid else np.nan,
            classical_SE=float(classic_se[j]),classical_t=float(beta[j]/classic_se[j]) if classic_se[j]>0 else np.nan,
            classical_p=float(2*student_t.sf(abs(beta[j]/classic_se[j]),df)) if classic_se[j]>0 else np.nan,
            legacy_CR1_SE=float(old_cr1_se[j]),legacy_CR1_p=float(2*student_t.sf(abs(beta[j]/old_cr1_se[j]),n-1)) if old_cr1_se[j]>0 else np.nan,
            outcome_SD=float(sd[j]),df=df,nuisance_rank=q.shape[1],residual_label_sum_squares=den,
            unadjusted_mean_difference=float(y[label==1,j].mean()-y[label==0,j].mean()),
            residual_mean_contrast=float(yr[label==1,j].mean()-yr[label==0,j].mean())))
    return pd.DataFrame(rows)


def bh(table,planned_size=None,column='p_value'):
    result=table.copy()
    if column not in result:result[column]=np.nan
    p=result[column].to_numpy(float);n=len(p) if planned_size is None else planned_size
    if n<len(p):raise ValueError('Family smaller than observed tests')
    filled=np.r_[np.where(np.isfinite(p),p,1.),np.ones(n-len(p))]
    q=false_discovery_control(filled,method='bh')[:len(p)];q[~np.isfinite(p)]=np.nan
    result['q_value']=q;result['planned_family_size']=n
    return result


def beta_batch(y,x,labels,block=256):
    q,yr=projection(x,y);out=np.full((len(labels),y.shape[1]),np.nan)
    for start in range(0,len(labels),block):
        z=labels[start:start+block].astype(float);zr=z-(z@q)@q.T;den=(zr*zr).sum(axis=1)
        valid=(z.sum(axis=1)>=5)&((len(x)-z.sum(axis=1))>=5)&(den>1e-10)
        values=np.full((len(z),y.shape[1]),np.nan);values[valid]=(zr[valid]@yr)/den[valid,None]
        out[start:start+len(z)]=values
    return out
