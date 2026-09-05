# -*- coding: utf-8 -*-
"""扩展临床主表: 基线 + 纵向随访, 覆盖全部临床结局、域特异认知与生物标志。"""
import pandas as pd, numpy as np

CUR  = r"F:/Deep_learning_derived_ MRI regional brain age/PPMI_Curated_Data_Cut_Public_20260511.xlsx"
BASE = r"F:/PD_brainage_data"

def rd(p): return [l.strip() for l in open(p) if l.strip()]
pd_ids = [int(x) for x in rd(f"{BASE}/pd_done_idlist.txt")]
hc_ids = [int(x) for x in rd(f"{BASE}/hc_all_done.txt")]
group  = {**{i: "PD" for i in pd_ids}, **{i: "HC" for i in hc_ids}}

BASIC = ['PATNO','COHORT','EVENT_ID','age_at_visit','SEX','EDUCYRS','duration_yrs',
         'DOMSIDE','NHY','updrs3_score','updrs3_score_on','updrs_totscore','moca','LEDD']
# R1(4)(5) / R2(2): 功能与残疾结局；R1(1)(2) / R2(1): 亚型与起病年龄
CLIN  = ['ageonset','td_pigd','td_pigd_on','pigd','pigd_on','sym_tremor',
         'updrs1_score','updrs2_score','updrs4_score','MSEADLG','NHY_ON']
# R1(4) / R2(3): MoCA 之外的域特异认知与 MCI 转化
COG   = ['cogstate','MCI_testscores','COG_COMPOSITE_INT','COG_COMPOSITE_EXT',
         'hvlt_immediaterecall','hvlt_retention','hvlt_discrimination','DVT_RECOG_DISC_INDEX',
         'bjlot','lns','DVS_LNS','SDMTOTAL','DVT_SDM','DVT_SFTANIM','DVS_SFTANIM']
# R1(最后一点): 共病理与神经退变生物标志
BIO   = ['abeta','tau','ptau','asyn','CSFSAA','ptau217_plasma','bd_tau_plasma',
         'Lilly_ptau217p','IU_ABeta42_CSF','IU_ABeta40_CSF','IU_pTau181_CSF','NFL_CSF']
# 非运动症状(备用)
NMS   = ['upsit','rem','gds','stai','scopa']

cur = pd.read_excel(CUR, sheet_name="20260511")
want = BASIC + CLIN + COG + BIO + NMS
have = [c for c in want if c in cur.columns]
missing = [c for c in want if c not in cur.columns]
if missing: print("[!] curated 表中缺失:", missing)

sub = cur[cur['PATNO'].isin(group)].copy()
sub['group'] = sub['PATNO'].map(group)

# --- 基线表 ---
bl = sub[sub['EVENT_ID'] == 'BL'][have + ['group']].drop_duplicates('PATNO').copy()
bl.to_csv(f"{BASE}/master_clinical_extended.csv", index=False, encoding='utf-8-sig')

# --- 纵向表(全随访): 用于 H&Y 进展 / MCI 转化 / 认知下降 ---
LONG = ['PATNO','EVENT_ID','age_at_visit','NHY','updrs2_score','updrs3_score',
        'moca','MSEADLG','cogstate','COG_COMPOSITE_INT','group']
lg = sub[[c for c in LONG if c in sub.columns]].copy()
lg.to_csv(f"{BASE}/master_clinical_longitudinal.csv", index=False, encoding='utf-8-sig')

print(f"基线表: {len(bl)} 例  ({bl['group'].value_counts().to_dict()})")
print(f"纵向表: {len(lg)} 行, {lg['PATNO'].nunique()} 例, {lg['EVENT_ID'].nunique()} 个访视点\n")

print("=== 新增变量在基线的可用性 (PD / HC 非缺失例数) ===")
for tag, cols in [("临床功能", CLIN), ("认知", COG), ("生物标志", BIO)]:
    print(f"\n-- {tag} --")
    for c in cols:
        if c not in bl.columns: continue
        n_pd = bl.loc[bl.group == 'PD', c].notna().sum()
        n_hc = bl.loc[bl.group == 'HC', c].notna().sum()
        print(f"  {c:24s} PD={n_pd:4d}  HC={n_hc:4d}")
