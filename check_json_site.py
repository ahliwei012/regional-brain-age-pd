# -*- coding: utf-8 -*-
import json
paths = {
 "PD_100005": r"G:/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P1/T1Img/100005/100005.json",
 "HC_100890": r"G:/PPMI_HC_XING/BL/3DT1_only/100890/100890.json",
}
keys = ["Manufacturer","ManufacturersModelName","MagneticFieldStrength","DeviceSerialNumber",
        "StationName","InstitutionName","SoftwareVersions","SeriesDescription","ImagingFrequency",
        "ReceiveCoilName","InstitutionAddress"]
for name,p in paths.items():
    print("="*30, name, "="*30)
    try:
        d = json.load(open(p, encoding="utf-8"))
        for k in keys:
            if k in d: print(f"  {k}: {d[k]}")
        print("  全部键:", list(d.keys()))
    except Exception as e:
        print("  读取失败:", e)
