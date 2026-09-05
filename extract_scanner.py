# -*- coding: utf-8 -*-
"""从JSON提取每个受试的扫描仪信息, 交叉检验 扫描仪×分组 是否混杂。"""
import json, glob, pandas as pd, numpy as np

BASE = r"F:/PD_brainage_data"
CAND = [
 "G:/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P1/T1Img/{id}",
 "G:/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P2/T1Img/{id}",
 "G:/PPMI_HC_XING/BL/3DT1_only/{id}",
 "F:/Analysis_HC_BL_ppmi10_13/T1Img/{id}",
]
def find_json(i):
    for c in CAND:
        fs = glob.glob(c.format(id=i)+"/*.json")
        if fs: return fs[0]
    return None

cl = pd.read_csv(f"{BASE}/master_clinical.csv")[["PATNO","group"]]
rows=[]
for i in cl["PATNO"]:
    p = find_json(str(i))
    if not p: rows.append([i,None,None,None,None]); continue
    try:
        d = json.load(open(p,encoding="utf-8"))
        rows.append([i, d.get("Manufacturer"), d.get("ManufacturersModelName"),
                     d.get("DeviceSerialNumber"), d.get("MagneticFieldStrength")])
    except: rows.append([i,None,None,None,None])
sc = pd.DataFrame(rows, columns=["PATNO","Manuf","Model","Serial","Field"])
sc = sc.merge(cl, on="PATNO")
sc["scanner"] = sc["Manuf"].astype(str)+"_"+sc["Model"].astype(str)
sc.to_csv(f"{BASE}/scanner_info.csv", index=False, encoding="utf-8-sig")

print(f"提取到JSON: {sc['Manuf'].notna().sum()} / {len(sc)}")
print(f"\n=== 厂商分布 ===\n{sc['Manuf'].value_counts(dropna=False).to_string()}")
print(f"\n=== 唯一扫描仪台数(Serial): {sc['Serial'].nunique()}   (厂商×型号: {sc['scanner'].nunique()}) ===")
print(f"\n=== 厂商 × 分组 交叉表 (查混杂!) ===")
print(pd.crosstab(sc["Manuf"], sc["group"], margins=True).to_string())
print(f"\n=== 每台扫描仪(Serial)受试数分布 ===")
vc = sc["Serial"].value_counts()
print(f"  扫描仪台数={len(vc)}, 中位数样本/台={vc.median():.0f}, 最小={vc.min()}, <5例的台数={ (vc<5).sum() }")
print(f"\n=== 场强分布 ===\n{sc['Field'].value_counts(dropna=False).to_string()}")
