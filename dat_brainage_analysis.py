# -*- coding: utf-8 -*-
"""PD特异新分析: 区域皮层脑龄 vs 多巴胺能失神经(DaTSCAN SBR)
   A) 合理性检验: PD vs HC 壳核SBR
   B) 耦合: 区域brain-PAD 与 SBR 是否相关(皮层老化是否继发于多巴胺缺失)
   C) 核心增量效度: 控制SBR后 brain-PAD 是否仍解释 UPDRS-III / MoCA
   D) 偏侧化: SBR不对称 vs 发病侧 vs 区域脑龄"""
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
import warnings; warnings.filterwarnings("ignore")
BASE=r"F:/PD_brainage_data"
REG=["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN={"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attn","vnt":"Ventral attn","dft":"Default",
    "slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}

sb=pd.read_csv(f'{BASE}/DaTSCAN/Xing_Core_Lab_-_Quant_SBR_24Jul2026.csv')
order={'SC':0,'BL':1,'V04':2,'V06':3,'V10':4}
sb['ord']=sb.EVENT_ID.map(order).fillna(9)
sb=sb.sort_values('ord').drop_duplicates('PATNO')
cols=['PATNO','EVENT_ID','STRIATUM_REF_CWM','CAUDATE_REF_CWM','PUTAMEN_REF_CWM',
      'PUTAMEN_L_REF_CWM','PUTAMEN_R_REF_CWM','CAUDATE_L_REF_CWM','CAUDATE_R_REF_CWM']
sb=sb[[c for c in cols if c in sb.columns]].rename(columns={'EVENT_ID':'DAT_EVENT'})
df=pd.read_csv(f'{BASE}/analysis_with_thickness.csv').merge(sb,on='PATNO',how='inner')
df=df.rename(columns={'PUTAMEN_REF_CWM':'putamen','STRIATUM_REF_CWM':'striatum','CAUDATE_REF_CWM':'caudate'})
print(f"合并后 n={len(df)} (PD={(df.group=='PD').sum()}, HC={(df.group=='HC').sum()}); SBR访视 {df.DAT_EVENT.value_counts().to_dict()}")
pdf=df[df.group=='PD'].copy()

print("\n########## A) 合理性检验: PD vs HC 壳核SBR ##########")
for c in ['putamen','caudate','striatum']:
    a=df[df.group=='PD'][c].dropna(); b=df[df.group=='HC'][c].dropna()
    t,p=stats.ttest_ind(a,b,equal_var=False)
    print(f"  {c:9s} PD {a.mean():.2f} vs HC {b.mean():.2f}  t={t:.1f} p={p:.3g}")

def lme(data,outcome,terms,group='batch'):
    need=[outcome,'age_at_visit','SEX','EDUCYRS',group]+[t for t in terms if t in data.columns]
    d=data.dropna(subset=need).copy()
    f=f"{outcome} ~ {' + '.join(terms)} + age_at_visit + C(SEX) + EDUCYRS"
    m=smf.mixedlm(f,d,groups=d[group]).fit()
    return m,len(d)

print("\n########## B) 耦合: 区域brain-PAD ~ 壳核SBR (PD内, 校正年龄/性别/教育+扫描仪) ##########")
rows=[]
for r in REG:
    m,n=lme(pdf,f'{r}_PAD',['putamen'])
    rows.append([EN[r],m.params['putamen'],m.pvalues['putamen'],n])
B=pd.DataFrame(rows,columns=['network','beta_SBR','p','n']); B['fdr']=multipletests(B['p'],method='fdr_bh')[1]
print(B.round(4).to_string(index=False))
print("  FDR<0.05:",B[B.fdr<0.05]['network'].tolist() or "无 -> 皮层脑龄与多巴胺能失神经基本解耦")

print("\n########## C) 核心: 控制SBR后 brain-PAD 是否仍解释临床 ##########")
SIG={'updrs3_score':['vsl','vnt','slt','lng'],'moca':['fpn','vnt']}
out=[]
for oc,nets in SIG.items():
    print(f"\n=== {oc} ===")
    for r in nets:
        m1,n1=lme(pdf,oc,[f'{r}_z'])                 # 仅brain-PAD
        m2,n2=lme(pdf,oc,[f'{r}_z','putamen'])       # +SBR
        b1,p1=m1.params[f'{r}_z'],m1.pvalues[f'{r}_z']
        b2,p2=m2.params[f'{r}_z'],m2.pvalues[f'{r}_z']
        bs,ps=m2.params['putamen'],m2.pvalues['putamen']
        att=100*(1-abs(b2)/abs(b1)) if b1 else np.nan
        print(f"  {EN[r]:14s} PAD单独 b={b1:+.3f}(p={p1:.3f}) | +SBR后 PAD b={b2:+.3f}(p={p2:.3f}) 衰减{att:.0f}% | SBR b={bs:+.3f}(p={ps:.3f})")
        out.append([oc,EN[r],b1,p1,b2,p2,bs,ps,att,n2])
C=pd.DataFrame(out,columns=['outcome','network','b_PADonly','p_PADonly','b_PADadjDAT','p_PADadjDAT','b_SBR','p_SBR','atten_pct','n'])
C.to_csv(f'{BASE}/DAT_incremental_validity.csv',index=False,encoding='utf-8-sig')

print("\n########## D) 偏侧化: SBR不对称 vs 发病侧 ##########")
pdf['put_AI']=(pdf['PUTAMEN_R_REF_CWM']-pdf['PUTAMEN_L_REF_CWM'])/(pdf['PUTAMEN_R_REF_CWM']+pdf['PUTAMEN_L_REF_CWM'])
for ds,lab in [(1.0,'DOMSIDE=1'),(2.0,'DOMSIDE=2'),(3.0,'DOMSIDE=3')]:
    s=pdf[pdf.DOMSIDE==ds]['put_AI'].dropna()
    if len(s)>3: print(f"  {lab}: n={len(s)} 壳核不对称指数(R-L)/(R+L) 均值={s.mean():+.4f}")
g1=pdf[pdf.DOMSIDE==1.0]['put_AI'].dropna(); g2=pdf[pdf.DOMSIDE==2.0]['put_AI'].dropna()
if len(g1)>3 and len(g2)>3:
    t,p=stats.ttest_ind(g1,g2,equal_var=False)
    print(f"  DOMSIDE 1 vs 2 的不对称差异: t={t:.2f} p={p:.3g}  (应显著=内部验证通过)")
print("\nDONE")
