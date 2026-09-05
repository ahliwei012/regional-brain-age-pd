# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
ba = pd.read_csv(r"F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main/Regional_and_global_brain_ages.csv")
ba["PATNO"] = ba["subject"].str.extract(r"(\d+)").astype(int)
cl = pd.read_csv(r"F:/PD_brainage_data/clinical_test5.csv")[["PATNO","age_at_visit","updrs3_score","moca"]]
m = ba.merge(cl, on="PATNO")
m["PAD_global"] = m["global"] - m["age_at_visit"]
print("=== 反推特征 → 模型脑龄 vs 实际年龄 (5例PD, 未做delta校正) ===")
print(m[["PATNO","age_at_visit","global","PAD_global","updrs3_score","moca"]]
      .round(2).to_string(index=False))
r = np.corrcoef(m["age_at_visit"], m["global"])[0,1]
mae = np.abs(m["global"] - m["age_at_visit"]).mean()
print(f"\n全脑脑龄 vs 实际年龄:  Pearson r = {r:.3f}   MAE = {mae:.2f} 岁")
print(f"PAD 均值 = {m['PAD_global'].mean():+.2f} 岁 (系统偏移, delta校正可消)")
