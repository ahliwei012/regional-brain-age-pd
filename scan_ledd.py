# -*- coding: utf-8 -*-
import pandas as pd, re
path = r"F:\Deep_learning_derived_ MRI regional brain age\PPMI_Curated_Data_Cut_Public_20260511.xlsx"
xl = pd.ExcelFile(path)
print("SHEETS:", xl.sheet_names)
df = pd.read_excel(path, sheet_name=0)
print("MAIN shape:", df.shape)

pat = re.compile(r'ledd|levodopa|dopa|medic|pd_med|on_levo|drug|treat|dose', re.I)
print("\n用药/LEDD相关列:")
hits = [c for c in df.columns if pat.search(str(c))]
for c in hits: print("   ", c)
if not hits: print("   (无)")

print("\n关键分析列是否存在:")
for c in ['PATNO','COHORT','EVENT_ID','age_at_visit','SEX','EDUCYRS',
          'duration_yrs','DOMSIDE','NHY','updrs3_score','updrs3_score_on','moca']:
    print(f"   {c}: {c in df.columns}")
