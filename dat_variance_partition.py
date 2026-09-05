# -*- coding: utf-8 -*-
"""方差分解: 临床严重度的变异中, 区域脑龄 vs 多巴胺能失神经(SBR) 各解释多少。
   commonality analysis: unique_PAD / unique_SBR / shared, + 嵌套模型F检验。
   base = age + sex + edu + scanner(batch) 固定效应; 所有模型同一子样本。"""
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import warnings; warnings.filterwarnings("ignore")
BASE=r"F:/PD_brainage_data"
REG={"vsl":"Visual","vnt":"Ventral attn","slt":"Salience","lng":"Language","fpn":"Frontoparietal"}

sb=pd.read_csv(f'{BASE}/DaTSCAN/Xing_Core_Lab_-_Quant_SBR_24Jul2026.csv')
order={'SC':0,'BL':1,'V04':2,'V06':3,'V10':4}
sb['ord']=sb.EVENT_ID.map(order).fillna(9)
sb=sb.sort_values('ord').drop_duplicates('PATNO')[['PATNO','PUTAMEN_REF_CWM']].rename(columns={'PUTAMEN_REF_CWM':'putamen'})
df=pd.read_csv(f'{BASE}/analysis_with_thickness.csv').merge(sb,on='PATNO',how='inner')
pdf=df[df.group=='PD'].copy()

BASEF="age_at_visit + C(SEX) + EDUCYRS + C(batch)"
rows=[]
for oc,ocn,nets in [("updrs3_score","UPDRS-III",["vsl","vnt","slt","lng"]),
                    ("moca","MoCA",["fpn","vnt"])]:
    print(f"\n{'='*74}\n{ocn}: 方差分解 (同一子样本)\n{'='*74}")
    for r in nets:
        pad=f"{r}_z"
        d=pdf.dropna(subset=[oc,pad,'putamen','age_at_visit','SEX','EDUCYRS','batch']).copy()
        m_base=smf.ols(f"{oc} ~ {BASEF}",d).fit()
        m_pad =smf.ols(f"{oc} ~ {BASEF} + {pad}",d).fit()
        m_sbr =smf.ols(f"{oc} ~ {BASEF} + putamen",d).fit()
        m_both=smf.ols(f"{oc} ~ {BASEF} + {pad} + putamen",d).fit()
        R0,Rp,Rs,Rb=m_base.rsquared,m_pad.rsquared,m_sbr.rsquared,m_both.rsquared
        uniq_pad=Rb-Rs; uniq_sbr=Rb-Rp
        total=Rb-R0; shared=total-uniq_pad-uniq_sbr
        # 嵌套F检验: PAD的独有贡献 (both vs sbr); SBR的独有贡献 (both vs pad)
        f_pad=anova_lm(m_sbr,m_both); p_pad=f_pad['Pr(>F)'].iloc[1]
        f_sbr=anova_lm(m_pad,m_both); p_sbr=f_sbr['Pr(>F)'].iloc[1]
        print(f"\n  {REG[r]}  (n={len(d)})")
        print(f"    base R²={R0:.4f} | +PAD {Rp:.4f} | +SBR {Rs:.4f} | +both {Rb:.4f}")
        print(f"    ── 脑龄独有 ΔR² = {uniq_pad*100:5.2f}%   (F检验 p={p_pad:.4f})")
        print(f"    ── DAT 独有 ΔR² = {uniq_sbr*100:5.2f}%   (F检验 p={p_sbr:.4f})")
        print(f"    ── 两者共享 ΔR² = {shared*100:5.2f}%")
        print(f"    ── 合计新增   = {total*100:5.2f}%   |  脑龄:DAT 独有贡献比 = {uniq_pad/uniq_sbr:.2f}" if uniq_sbr>1e-9 else "")
        rows.append([ocn,REG[r],len(d),R0,Rp,Rs,Rb,uniq_pad,uniq_sbr,shared,total,p_pad,p_sbr])
T=pd.DataFrame(rows,columns=['outcome','network','n','R2_base','R2_PAD','R2_SBR','R2_both',
    'unique_PAD','unique_SBR','shared','total_added','p_PAD_unique','p_SBR_unique'])
T.to_csv(f'{BASE}/DAT_variance_partition.csv',index=False,encoding='utf-8-sig')
print(f"\n\n{'='*74}\n小结\n{'='*74}")
for ocn in T.outcome.unique():
    s=T[T.outcome==ocn]
    print(f"  {ocn}: 脑龄独有均值 {s.unique_PAD.mean()*100:.2f}% vs DAT独有均值 {s.unique_SBR.mean()*100:.2f}%")
print("\nDONE -> DAT_variance_partition.csv")
