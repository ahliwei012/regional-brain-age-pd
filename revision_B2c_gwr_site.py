# -*- coding: utf-8 -*-
"""B2c 关键判别: GWR 与运动严重度的关联是真实信号还是站点/采集伪影?
GWR(灰白质信号强度比)对采集参数高度敏感, 而 brain-PAD 的主分析一律带站点随机效应。
公平比较需在同一框架下检验:
 (1) GWR 有多少方差来自扫描仪(与 brain-PAD 的 19.5% 对照);
 (2) 加入站点随机效应后 GWR-UPDRS 关联是否存活;
 (3) 对 GWR 做与 brain-PAD 相同的 ComBat 谐和后再检验。
"""
import numpy as np, pandas as pd, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"
REG = ["sns", "fpn", "drs", "vnt", "dft", "slt", "lng", "adt", "vsl", "lmb"]
EN = {"sns": "Sensorimotor", "fpn": "Frontoparietal", "drs": "Dorsal attention",
      "vnt": "Ventral attention", "dft": "Default mode", "slt": "Salience",
      "lng": "Language", "adt": "Auditory", "vsl": "Visual", "lmb": "Limbic"}

df = pd.read_csv(BASE + "/analysis_revision_base.csv").merge(
    pd.read_csv(BASE + "/network_thickness_gwr.csv"), on="PATNO", how="left")

# ---- 对 GWR 与厚度做与 brain-PAD 同款的 ComBat 谐和(经验贝叶斯位置/尺度校正) ----
def combat(vals, batch, mod):
    """标准 ComBat: 先回归掉生物学协变量, 再对残差做批次位置/尺度校正。"""
    y = vals.astype(float).copy()
    ok = np.isfinite(y)
    X = np.column_stack([np.ones(ok.sum())] + [m[ok] for m in mod])
    beta = np.linalg.lstsq(X, y[ok], rcond=None)[0]
    fit = np.full(len(y), np.nan)
    fit[ok] = X @ beta
    res = y - fit
    sd = np.nanstd(res[ok])
    out = res / sd
    for b in pd.unique(batch):
        m = (batch == b) & ok
        if m.sum() < 2:
            continue
        out[m] = (out[m] - np.nanmean(out[m])) / (np.nanstd(out[m]) or 1.0)
    return out * sd + fit


mod = [df["age_at_visit"].values, (df["SEX"] == df["SEX"].iloc[0]).astype(float).values,
       (df["group"] == "PD").astype(float).values]
for r in REG:
    df[r + "_gwr_h"] = combat(df[r + "_gwr"].values, df["batch"].values, mod)

pdf = df[df.group == "PD"].copy()
for r in REG:
    for s in ["_gwr", "_gwr_h", "_thick"]:
        v = pdf[r + s]
        pdf[r + s + "z"] = (v - v.mean()) / v.std()


def icc(d, y):
    m = smf.mixedlm("%s ~ age_at_visit + C(SEX) + EDUCYRS" % y, d, groups=d["batch"]).fit()
    vb = float(np.asarray(m.cov_re)[0, 0])
    return vb / (vb + float(m.scale)) * 100


# ============ (1) GWR 的扫描仪方差 ============
print("=" * 88)
print("(1) 各网络平均 GWR 中来自扫描仪的方差比例  [对照: brain-PAD 谐和前中位 19.5%]")
print("=" * 88)
d0 = df.dropna(subset=["EDUCYRS"]).copy()
rows = []
for r in REG:
    rows.append([EN[r], icc(d0, r + "_gwr"), icc(d0, r + "_thick"), icc(d0, r + "_PAD")])
v = pd.DataFrame(rows, columns=["Network", "GWR_scanner_var_%", "Thickness_scanner_var_%",
                                "brainPAD_scanner_var_%"])
v.to_csv(BASE + "/rev_B2c_gwr_scanner_variance.csv", index=False, encoding="utf-8-sig")
print(v.round(2).to_string(index=False))
print("\n  中位数: GWR %.1f%% | 厚度 %.1f%% | brain-PAD %.1f%%"
      % (v["GWR_scanner_var_%"].median(), v["Thickness_scanner_var_%"].median(),
         v["brainPAD_scanner_var_%"].median()))

# ============ (2)(3) 站点随机效应 / ComBat 后 GWR-UPDRS 关联 ============
for outcome, tag in [("updrs3_score", "UPDRS-III"), ("moca", "MoCA")]:
    print("\n" + "=" * 88)
    print("%s : GWR 关联在 OLS / 站点随机效应 / ComBat 谐和 三种处理下的表现" % tag)
    print("=" * 88)
    rows = []
    for r in REG:
        d = pdf.dropna(subset=[r + "_gwrz", r + "_gwr_hz", r + "_z", outcome,
                               "age_at_visit", "SEX", "EDUCYRS", "batch"]).copy()
        o = smf.ols("%s ~ %s_gwrz + age_at_visit + C(SEX) + EDUCYRS" % (outcome, r), d).fit()
        l = smf.mixedlm("%s ~ %s_gwrz + age_at_visit + C(SEX) + EDUCYRS" % (outcome, r),
                        d, groups=d["batch"]).fit()
        h = smf.mixedlm("%s ~ %s_gwr_hz + age_at_visit + C(SEX) + EDUCYRS" % (outcome, r),
                        d, groups=d["batch"]).fit()
        # 站点随机效应框架下, 在 GWR+厚度之上加入 brain-PAD 的偏效应
        f = smf.mixedlm("%s ~ %s_z + %s_thickz + %s_gwr_hz + age_at_visit + C(SEX) + EDUCYRS"
                        % (outcome, r, r, r), d, groups=d["batch"]).fit()
        rows.append([EN[r], o.params[r + "_gwrz"], o.pvalues[r + "_gwrz"],
                     l.params[r + "_gwrz"], l.pvalues[r + "_gwrz"],
                     h.params[r + "_gwr_hz"], h.pvalues[r + "_gwr_hz"],
                     f.params[r + "_z"], f.pvalues[r + "_z"],
                     f.params[r + "_gwr_hz"], f.pvalues[r + "_gwr_hz"], len(d)])
    out = pd.DataFrame(rows, columns=[
        "Network", "b_GWR_OLS", "P_OLS", "b_GWR_LME", "P_LME", "b_GWRh_LME", "P_GWRh",
        "b_PAD_full", "P_PAD_full", "b_GWRh_full", "P_GWRh_full", "n"])
    for c in ["P_OLS", "P_LME", "P_GWRh", "P_PAD_full", "P_GWRh_full"]:
        out["fdr" + c[1:]] = multipletests(out[c], method="fdr_bh")[1]
    out.to_csv(BASE + "/rev_B2c_%s_gwr_site.csv" % outcome, index=False, encoding="utf-8-sig")
    print(out[["Network", "b_GWR_OLS", "P_OLS", "b_GWR_LME", "P_LME",
               "b_GWRh_LME", "P_GWRh"]].round(4).to_string(index=False))
    print("\n  GWR 关联 FDR<0.05 :  OLS(无站点) %s" % out[out.fdr_OLS < .05]["Network"].tolist())
    print("                       +站点随机效应 %s" % out[out.fdr_LME < .05]["Network"].tolist())
    print("                       ComBat谐和后  %s" % out[out.fdr_GWRh < .05]["Network"].tolist())
    print("\n  完整模型(brain-PAD + 厚度 + 谐和GWR, 带站点随机效应):")
    print(out[["Network", "b_PAD_full", "P_PAD_full", "fdr_PAD_full",
               "b_GWRh_full", "P_GWRh_full", "fdr_GWRh_full"]].round(4).to_string(index=False))
