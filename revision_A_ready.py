# -*- coding: utf-8 -*-
"""补充分析 (第一批, 基于已生成的数据)
B1 各网络MAE / B3 HC组内关联 / B4 扫描仪方差分解 / B6 亚型与起病年龄
+ 年龄匹配敏感性 + MoCA控制运动severity。
"""
import pandas as pd, numpy as np, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REG  = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
ALL  = REG + ["global", "AD"]
EN   = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
        "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
        "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic",
        "global": "Global"}


def star(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


# ---------------- 数据装配 ----------------
# 沿用主分析数据集(含 ComBat 谐和后的 _h 列), 保证与正文数字一致;
# 仅补入本批分析新需要的临床变量与场强/型号元数据。
df = pd.read_csv(BASE + "/analysis_merged_final.csv")
NEWCL = [c for c in pd.read_csv(BASE + "/master_clinical_extended.csv", nrows=1).columns
         if c not in df.columns or c == "PATNO"]
cl = pd.read_csv(BASE + "/master_clinical_extended.csv")[NEWCL]
sc = pd.read_csv(BASE + "/scanner_info.csv")[["PATNO", "Model", "Serial", "Field"]]
df = df.merge(cl, on="PATNO", how="left").merge(sc, on="PATNO", how="left")


def norm_model(m):
    """PPMI 型号名在不同站点写法不一(空格/下划线/MAGNETOM前缀)。用于说明 ComBat 批次的构建规则。"""
    return str(m).upper().replace("MAGNETOM", "").replace("-", "").replace("_", "").replace(" ", "").strip()


df["Model_norm"] = df["Model"].map(norm_model)
vc = df["scanner"].value_counts()

print("样本 %d (PD=%d, HC=%d)" % (len(df), sum(df.group == "PD"), sum(df.group == "HC")))
print("扫描仪: 型号名称 %d 个 -> 拼写规范化后 %d 个实际型号 -> ComBat 批次 %d 个"
      % (df["Model"].nunique(), df["Model_norm"].nunique(), df["batch"].nunique()))
print("  (同一型号在不同站点被记录为不同名称时仍作为独立批次, 因此批次实际带有站点粒度)")
print("  n<10 的型号: %s" % vc[vc < 10].to_dict())
print("场强分布: %s\n" % df["Field"].value_counts(dropna=False).to_dict())

hc = df[df.group == "HC"]
hc2 = hc
for c in REG + ["global"]:
    a, b = np.polyfit(hc2["age_at_visit"], hc2[c + "_h"], 1)
    df[c + "_PADh"] = df[c + "_h"] - (a * df["age_at_visit"] + b)
    df[c + "_zh"] = (df[c + "_PADh"] - df[c + "_PADh"].mean()) / df[c + "_PADh"].std()

hc = df[df.group == "HC"]
pdf = df[df.group == "PD"].copy()


def fam(data, outcome, zsuf="_z", extra="", extra_vars=()):
    rows = []
    for c in REG + ["global"]:
        d = data.dropna(subset=[c + zsuf, outcome, "age_at_visit", "SEX", "EDUCYRS",
                                "batch"] + list(extra_vars)).rename(columns={c + zsuf: "x"})
        m = smf.mixedlm("%s ~ x + age_at_visit + C(SEX) + EDUCYRS%s" % (outcome, extra),
                        d, groups=d["batch"]).fit()
        ci = m.conf_int().loc["x"]
        rows.append([EN[c], m.params["x"], ci[0], ci[1], m.pvalues["x"], len(d)])
    r = pd.DataFrame(rows, columns=["Network", "beta", "lo", "hi", "p", "n"])
    r["fdr"] = multipletests(r["p"], method="fdr_bh")[1]
    r["sig"] = r["fdr"].map(star)
    return r


# ================= B1: 各网络 MAE (R2) =================
print("=" * 72)
print("B1  各皮层网络的年龄预测精度 (%d 名 HC)" % len(hc))
print("=" * 72)
# 稿件中的 "4.47 年" 是年龄偏倚校正后的残差 MAE, 与原始预测误差不同;
# 两者并列报告, 并明确标注定义。
rng = np.random.default_rng(42)
rows = []
for c in REG + ["global"]:
    yhat = hc[c].values.astype(float)
    age = hc["age_at_visit"].values.astype(float)
    err = np.abs(yhat - age)
    a, b = np.polyfit(age, yhat, 1)
    err_c = np.abs(yhat - (a * age + b))          # 去偏后残差
    idx = rng.integers(0, len(err), (1000, len(err)))
    bs = err_c[idx].mean(axis=1)
    rows.append([EN[c], err.mean(), err_c.mean(), np.percentile(bs, 2.5),
                 np.percentile(bs, 97.5), np.corrcoef(age, yhat)[0, 1],
                 float(np.sqrt(((yhat - age) ** 2).mean())), (yhat - age).mean(), a])
b1 = pd.DataFrame(rows, columns=["Network", "MAE_raw", "MAE_biascorr", "MAEc_lo",
                                 "MAEc_hi", "r", "RMSE_raw", "mean_offset", "slope"])
b1.to_csv(BASE + "/rev_B1_network_MAE.csv", index=False, encoding="utf-8-sig")
print(b1.round(3).to_string(index=False))
print("\n  注: 稿件正文的 MAE=4.47 对应本表 Global 的 MAE_biascorr;")
print("      原始预测误差 (MAE_raw) 更大, 因模型存在系统性偏移, 本文需明确区分两者。")

# ================= B4: 扫描仪/场强方差分解 (R1-6, R2-3) =================
print("\n" + "=" * 72)
print("B4  扫描仪与场强解释的 brain-PAD 方差 (ComBat 前 vs 后)")
print("=" * 72)


def icc_batch(d, y):
    m = smf.mixedlm("%s ~ age_at_visit + C(SEX) + EDUCYRS" % y, d, groups=d["batch"]).fit()
    vb = float(np.asarray(m.cov_re)[0, 0])
    return vb / (vb + float(m.scale))


d4 = df.dropna(subset=["EDUCYRS", "Field"]).copy()
d4["F15"] = (d4["Field"] < 2).astype(int)
rows = []
for c in REG + ["global"]:
    pre, post = icc_batch(d4, c + "_PAD"), icc_batch(d4, c + "_PADh")
    f1 = smf.ols("%s_PAD ~ F15 + age_at_visit + C(SEX)" % c, d4).fit()
    f2 = smf.ols("%s_PADh ~ F15 + age_at_visit + C(SEX)" % c, d4).fit()
    r0 = smf.ols("%s_PAD ~ age_at_visit + C(SEX)" % c, d4).fit().rsquared
    r0h = smf.ols("%s_PADh ~ age_at_visit + C(SEX)" % c, d4).fit().rsquared
    rows.append([EN[c], pre * 100, post * 100, (f1.rsquared - r0) * 100,
                 (f2.rsquared - r0h) * 100, f1.pvalues["F15"], f2.pvalues["F15"]])
b4 = pd.DataFrame(rows, columns=["Network", "Scanner_var_pre_%", "Scanner_var_post_%",
                                 "Field_R2_pre_%", "Field_R2_post_%",
                                 "P_field_pre", "P_field_post"])
b4.to_csv(BASE + "/rev_B4_scanner_variance.csv", index=False, encoding="utf-8-sig")
print(b4.round(3).to_string(index=False))
print("\n  扫描仪方差中位数: ComBat 前 %.2f%% -> 后 %.2f%%"
      % (b4["Scanner_var_pre_%"].median(), b4["Scanner_var_post_%"].median()))

# ---- B4b: 剔除 n<10 的扫描仪后重跑主关联 (R2 明确建议) ----
print("\n" + "-" * 72)
print("B4b 敏感性: 仅保留 n>=10 的扫描仪型号")
print("-" * 72)
keep = vc[vc >= 10].index
sub = pdf[pdf["scanner"].isin(keep)].copy()
print("  保留 %d 个型号, PD 例数 %d / %d" % (len(keep), len(sub), len(pdf)))
for oc, tag in [("updrs3_score", "UPDRS-III"), ("moca", "MoCA")]:
    r = fam(sub, oc)
    r.to_csv(BASE + "/rev_B4b_%s_bigscanners.csv" % oc, index=False, encoding="utf-8-sig")
    print("\n  [%s] FDR<0.05: %s" % (tag, r[r.fdr < 0.05]["Network"].tolist()))
    print(r[["Network", "beta", "p", "fdr", "sig", "n"]].round(4).to_string(index=False))

# ================= B3: HC 组内关联 + 组别交互 (R2) =================
print("\n" + "=" * 72)
print("B3  关联是否为 PD 特异: HC 组内 MoCA 关联 + group x brain-PAD 交互")
print("=" * 72)
rows = []
for c in REG + ["global"]:
    dh = hc.dropna(subset=[c + "_z", "moca", "age_at_visit", "SEX", "EDUCYRS",
                           "batch"]).rename(columns={c + "_z": "x"})
    mh = smf.mixedlm("moca ~ x + age_at_visit + C(SEX) + EDUCYRS", dh, groups=dh["batch"]).fit()
    da = df.dropna(subset=[c + "_z", "moca", "age_at_visit", "SEX", "EDUCYRS",
                           "batch"]).rename(columns={c + "_z": "x"}).copy()
    da["PD"] = (da.group == "PD").astype(int)
    mi = smf.mixedlm("moca ~ x*PD + age_at_visit + C(SEX) + EDUCYRS", da, groups=da["batch"]).fit()
    rows.append([EN[c], mh.params["x"], mh.pvalues["x"], len(dh),
                 mi.params.get("x:PD", np.nan), mi.pvalues.get("x:PD", np.nan)])
b3 = pd.DataFrame(rows, columns=["Network", "beta_HC_MoCA", "P_HC", "n_HC",
                                 "beta_interaction", "P_interaction"])
b3["fdr_HC"] = multipletests(b3["P_HC"], method="fdr_bh")[1]
b3.to_csv(BASE + "/rev_B3_HC_association.csv", index=False, encoding="utf-8-sig")
print(b3.round(4).to_string(index=False))

# ================= 年龄匹配敏感性 (R2 关于 Table 1 年龄差异) =================
print("\n" + "=" * 72)
print("年龄匹配子样本: PD vs HC 组间比较")
print("=" * 72)
# HC (n=261) 是较小的组, 因此以每名 HC 为基准, 从 PD 池中 1:1 无放回匹配最近年龄者
pdg, hcg = df[df.group == "PD"], df[df.group == "HC"]
used, pairs = set(), []
for j, row in hcg.sample(frac=1, random_state=42).iterrows():
    cand = pdg[(~pdg.index.isin(used)) & ((pdg.age_at_visit - row.age_at_visit).abs() <= 1)]
    if len(cand):
        i = (cand.age_at_visit - row.age_at_visit).abs().idxmin()
        used.add(i)
        pairs.append((i, j))
mi_df = df.loc[[p[0] for p in pairs] + [p[1] for p in pairs]].copy()
tt = stats.ttest_ind(mi_df[mi_df.group == "PD"].age_at_visit, mi_df[mi_df.group == "HC"].age_at_visit)
print("  匹配对数 %d; 匹配后年龄 PD %.1f vs HC %.1f, P = %.3f"
      % (len(pairs), mi_df[mi_df.group == "PD"].age_at_visit.mean(),
         mi_df[mi_df.group == "HC"].age_at_visit.mean(), tt.pvalue))
rows = []
for c in REG + ["global"]:
    d = mi_df.dropna(subset=[c + "_PAD", "EDUCYRS"]).copy()
    d["PD"] = (d.group == "PD").astype(int)
    m = smf.mixedlm("%s_PAD ~ PD + age_at_visit + C(SEX) + EDUCYRS" % c, d, groups=d["batch"]).fit()
    ci = m.conf_int().loc["PD"]
    rows.append([EN[c], m.params["PD"], ci[0], ci[1], m.pvalues["PD"], len(d)])
am = pd.DataFrame(rows, columns=["Network", "beta", "lo", "hi", "p", "n"])
am["fdr"] = multipletests(am["p"], method="fdr_bh")[1]
am["sig"] = am["fdr"].map(star)
am.to_csv(BASE + "/rev_agematched_pdvshc.csv", index=False, encoding="utf-8-sig")
print(am.round(4).to_string(index=False))

# ================= MoCA 控制运动严重度 (R2) =================
print("\n" + "=" * 72)
print("MoCA 关联在加入 UPDRS-III 为协变量后是否存活")
print("=" * 72)
r = fam(pdf, "moca", extra=" + updrs3_score", extra_vars=("updrs3_score",))
r.to_csv(BASE + "/rev_moca_adj_updrs3.csv", index=False, encoding="utf-8-sig")
print(r[["Network", "beta", "p", "fdr", "sig", "n"]].round(4).to_string(index=False))

# ================= B6: 亚型与起病年龄 =================
print("\n" + "=" * 72)
print("B6  TD/PIGD 亚型与起病年龄")
print("=" * 72)
print("  td_pigd 分布: %s" % pdf["td_pigd"].value_counts(dropna=False).to_dict())
rows = []
for c in REG + ["global"]:
    d = pdf.dropna(subset=[c + "_PAD", "td_pigd", "EDUCYRS"]).copy()
    m = smf.mixedlm("%s_PAD ~ C(td_pigd) + age_at_visit + C(SEX) + EDUCYRS" % c,
                    d, groups=d["batch"]).fit()
    k = [p for p in m.params.index if p.startswith("C(td_pigd)")][0]
    d2 = pdf.dropna(subset=[c + "_z", "updrs3_score", "ageonset",
                            "EDUCYRS"]).rename(columns={c + "_z": "x"})
    m2 = smf.mixedlm("updrs3_score ~ x*ageonset + age_at_visit + C(SEX) + EDUCYRS",
                     d2, groups=d2["batch"]).fit()
    rows.append([EN[c], m.params[k], m.pvalues[k], len(d),
                 m2.params.get("x:ageonset", np.nan), m2.pvalues.get("x:ageonset", np.nan)])
b6 = pd.DataFrame(rows, columns=["Network", "beta_PIGDvsTD", "P_subtype", "n",
                                 "beta_x_ageonset", "P_x_ageonset"])
b6["fdr_subtype"] = multipletests(b6["P_subtype"], method="fdr_bh")[1]
b6.to_csv(BASE + "/rev_B6_subtype_ageonset.csv", index=False, encoding="utf-8-sig")
print(b6.round(4).to_string(index=False))

df.to_csv(BASE + "/analysis_revision_base.csv", index=False, encoding="utf-8-sig")
print("\n[已保存] analysis_revision_base.csv 及 rev_*.csv")
