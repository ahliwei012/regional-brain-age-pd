# -*- coding: utf-8 -*-
"""为测试批次的5例PD抽取基线(BL)临床数据，生成分析表。"""
import pandas as pd

path = r"<path to>/PPMI_Curated_Data_Cut_Public.xlsx"
out = r"F:\PD_brainage_data\clinical_test5.csv"
ids = []  # insert your own PPMI subject IDs (PATNO) here

df = pd.read_excel(path, sheet_name='20260316')

# 想要的列（仅保留实际存在的）
wanted = ['PATNO', 'COHORT', 'EVENT_ID', 'age_at_visit', 'SEX', 'EDUCYRS',
          'agediag', 'ageonset', 'duration_yrs', 'DOMSIDE', 'NHY', 'hy',
          'updrs3_score', 'updrs3_score_on', 'updrs_totscore', 'moca']
keep = [c for c in wanted if c in df.columns]

sub = df[df['PATNO'].isin(ids)]
bl = sub[sub['EVENT_ID'] == 'BL'][keep].sort_values('PATNO').reset_index(drop=True)

# COHORT 含义：PPMI 中 1=PD, 2=HC, 等
bl.to_csv(out, index=False, encoding='utf-8-sig')
print("已保存:", out)
print("\n=== 5例基线临床分析表 ===")
print(bl.to_string(index=False))
print("\n字段说明: updrs3_score=运动(OFF药), updrs3_score_on=运动(ON药), moca=认知, "
      "DOMSIDE=优势侧, duration_yrs=病程(年), NHY=Hoehn-Yahr分期")
