# -*- coding: utf-8 -*-
"""全队列临床主表: 560 PD + 61 HC 的基线结局/协变量。"""
import pandas as pd, numpy as np

CUR = r"F:/Deep_learning_derived_ MRI regional brain age/PPMI_Curated_Data_Cut_Public_20260511.xlsx"
BASE = r"F:/PD_brainage_data"

def rd(p): return [l.strip() for l in open(p) if l.strip()]
pd_ids = [int(x) for x in rd(f"{BASE}/pd_done_idlist.txt")]
hc_ids = [int(x) for x in rd(f"{BASE}/hc_all_done.txt")]
group = {**{i: "PD" for i in pd_ids}, **{i: "HC" for i in hc_ids}}
allids = list(group)
print(f"目标: PD={len(pd_ids)}  HC={len(hc_ids)}  合计={len(allids)}")

cur = pd.read_excel(CUR, sheet_name="20260511")
cols = ['PATNO','COHORT','EVENT_ID','age_at_visit','SEX','EDUCYRS','duration_yrs',
        'DOMSIDE','NHY','updrs3_score','updrs3_score_on','updrs_totscore','moca','LEDD']
cols = [c for c in cols if c in cur.columns]
sub = cur[cur['PATNO'].isin(allids)]
bl = sub[sub['EVENT_ID'] == 'BL'][cols].drop_duplicates('PATNO').copy()
bl['group'] = bl['PATNO'].map(group)

print(f"匹配到基线记录: {len(bl)} / {len(allids)}")
miss = set(allids) - set(bl['PATNO'])
print(f"缺基线记录的影像ID数: {len(miss)}")
print("\n各组计数:", bl['group'].value_counts().to_dict())
print("COHORT分布:", bl['COHORT'].value_counts().to_dict())
print("\n=== 关键变量缺失率 ===")
for c in ['age_at_visit','SEX','EDUCYRS','updrs3_score','moca','LEDD','DOMSIDE','duration_yrs']:
    if c in bl.columns:
        print(f"  {c}: 缺失 {bl[c].isna().sum()}/{len(bl)}")
print("\n=== PD vs HC 描述(均值) ===")
for c in ['age_at_visit','EDUCYRS','updrs3_score','moca']:
    if c in bl.columns:
        g = bl.groupby('group')[c].mean()
        print(f"  {c}: PD={g.get('PD',np.nan):.2f}  HC={g.get('HC',np.nan):.2f}")
out = f"{BASE}/master_clinical.csv"
bl.to_csv(out, index=False, encoding='utf-8-sig')
print(f"\n已保存: {out}")
