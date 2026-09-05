# -*- coding: utf-8 -*-
import pandas as pd, re
path = r"<path to>/PPMI_Curated_Data_Cut_Public.xlsx"

xl = pd.ExcelFile(path)
print("SHEETS:", xl.sheet_names)

df = pd.read_excel(path, sheet_name=0)
print("MAIN SHEET shape:", df.shape)

pat = re.compile(r'updrs|np3|moca|montreal|age|sex|gender|patno|event|educ|durat|hoehn|hy|domin|onset|side|cohort|diag|visit|enroll', re.I)
print("\nRELEVANT COLUMNS:")
for c in df.columns:
    if pat.search(str(c)):
        print("   ", c)
print("\nTOTAL COLUMNS:", len(df.columns))

ids = []  # insert your own PPMI subject IDs (PATNO) here
# find the subject-id column
idcol = None
for cand in ['PATNO', 'patno', 'SUBJECT', 'subject_id', 'ID']:
    if cand in df.columns:
        idcol = cand; break
print("\nID column:", idcol)
if idcol:
    print("id dtype:", df[idcol].dtype)
    sub = df[df[idcol].astype(str).str.strip().isin([str(i) for i in ids])]
    print("rows matching our 5 test IDs:", len(sub))
    evcol = 'EVENT_ID' if 'EVENT_ID' in df.columns else None
    if evcol:
        print("EVENT_ID values for matched:", sorted(sub[evcol].astype(str).unique().tolist()))
    print("unique matched IDs:", sorted(sub[idcol].astype(str).unique().tolist()))
