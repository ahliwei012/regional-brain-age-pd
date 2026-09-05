# -*- coding: utf-8 -*-
"""决定性验证: 61 HC 反推特征→脑龄 vs 实际年龄, 算MAE/r + delta校正。"""
import pandas as pd, numpy as np

import sys
_csv = sys.argv[1] if len(sys.argv) > 1 else r"F:/PD_brainage_data/brain_ages/hc_brain_ages.csv"
ba = pd.read_csv(_csv)
ba["PATNO"] = ba["subject"].str.extract(r"(\d+)").astype(int)

cur = pd.read_excel(r"F:/Deep_learning_derived_ MRI regional brain age/PPMI_Curated_Data_Cut_Public_20260511.xlsx",
                    sheet_name="20260511")
hc = cur[cur["PATNO"].isin(ba["PATNO"])]
bl = hc[hc["EVENT_ID"] == "BL"][["PATNO", "age_at_visit", "COHORT"]].drop_duplicates("PATNO")
m = ba.merge(bl, on="PATNO").dropna(subset=["age_at_visit"])
print(f"匹配到基线年龄的HC: {len(m)} / {len(ba)}")
print(f"COHORT分布: {m['COHORT'].value_counts().to_dict()}  (PPMI: 1=PD,2=HC)")

x = m["age_at_visit"].values.astype(float)
y = m["global"].values.astype(float)
mae = np.abs(y - x).mean(); r = np.corrcoef(x, y)[0, 1]
print(f"\n=== 原始(未校正) 全脑脑龄 vs 实际年龄  (n={len(m)}) ===")
print(f"  MAE = {mae:.2f} 岁    Pearson r = {r:.3f}    PAD均值 = {(y-x).mean():+.2f} 岁")

slope, intercept = np.polyfit(x, y, 1)
resid = y - (slope * x + intercept)
mae_c = np.abs(resid).mean()
print(f"\n=== delta校正后 (Smith 2019线性去偏) ===")
print(f"  拟合 brain = {slope:.3f}*age + {intercept:.2f}")
print(f"  校正后 MAE(残差) = {mae_c:.2f} 岁    PAD均值 = {resid.mean():+.2f} 岁")

print(f"\n参照: 论文模型健康人 MAE 2.9-3.1 岁, r 0.88-0.90")
verdict = "通过✓ 复刻确证" if (mae_c < 4.0 and r > 0.75) else ("接近, 需微调" if r > 0.6 else "未通过, 需检查")
print(f"\n>>> 判定: MAE校正后={mae_c:.2f}, r={r:.3f} -> {verdict}")

# 各区域脑龄的 r (看哪些区域信号最好)
print("\n各区域脑龄 vs 年龄 r:")
for col in ["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb","global"]:
    rr = np.corrcoef(x, m[col].values.astype(float))[0,1]
    print(f"  {col:7s}: r={rr:.3f}")
m.to_csv(r"F:/PD_brainage_data/brain_ages/hc_eval_merged.csv", index=False)
