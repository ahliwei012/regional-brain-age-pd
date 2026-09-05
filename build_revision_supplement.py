# -*- coding: utf-8 -*-
"""补充材料:
1) 新建 Supplementary Figure 6 (扫描仪/场强方差分解, 谐和前后);
2) 按重排后的编号重命名全部补充图;
3) 生成补充材料 Word 文档(21 张表 + 7 张图注 + 独立参考文献表);
4) 生成 Brain Communications 要求的"补充图与图注合一"单个 PDF。
"""
import os, json, shutil
os.environ["MPLCONFIGDIR"] = r"G:/Temp/mpl"
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
import docx
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
import warnings
warnings.filterwarnings("ignore")
Image.MAX_IMAGE_PIXELS = None

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})

BASE = r"F:/PD_brainage_data"
R1 = BASE + "/SUBMISSION_BrainComms_R1"
FIG = R1 + "/Figures"
OLDFIG = BASE + "/figures_previous"
BLUE, ORANGE, GREY = '#0F4D92', '#D55E00', '#9A9A9A'
mp = json.load(open(R1 + "/supp_renumber_map.json"))

# ============ 1. 新 Supplementary Figure 6: 方差分解 ============
V = pd.read_csv(BASE + "/rev_B4_scanner_variance.csv")
fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.4), gridspec_kw={'wspace': 0.42})
x = np.arange(len(V))
ax = axes[0]
ax.bar(x - 0.2, V["Scanner_var_pre_%"], 0.4, color=ORANGE, label='Before harmonization')
ax.bar(x + 0.2, V["Scanner_var_post_%"], 0.4, color=BLUE, label='After harmonization')
ax.set_xticks(x)
ax.set_xticklabels(V["Network"], rotation=60, ha='right', fontsize=6.5)
ax.set_ylabel('Scanner-attributable variance\nin brain-PAD (%)', fontsize=8)
ax.set_title('A  Scanner', fontsize=9.5, loc='left', fontweight='bold', color='black')
ax.legend(fontsize=7, frameon=False, labelcolor='black')
ax.tick_params(colors='black', labelsize=7.5)
ax.yaxis.label.set_color('black')

ax = axes[1]
ax.bar(x - 0.2, V["Field_R2_pre_%"], 0.4, color=ORANGE, label='Before harmonization')
ax.bar(x + 0.2, V["Field_R2_post_%"], 0.4, color=BLUE, label='After harmonization')
for i, (p1, p2) in enumerate(zip(V["P_field_pre"], V["P_field_post"])):
    if p1 < 0.05:
        ax.text(i - 0.2, V["Field_R2_pre_%"][i] + 0.08, '*', ha='center', fontsize=8, color='black')
    if p2 < 0.05:
        ax.text(i + 0.2, V["Field_R2_post_%"][i] + 0.08, '*', ha='center', fontsize=8, color='black')
ax.set_xticks(x)
ax.set_xticklabels(V["Network"], rotation=60, ha='right', fontsize=6.5)
ax.set_ylabel('Field-strength partial R\u00b2\nin brain-PAD (%)', fontsize=8)
ax.set_title('B  Field strength (3 T, n = 754 versus 1.5 T, n = 49)', fontsize=8.6,
             loc='left', fontweight='bold', color='black')
ax.legend(fontsize=7, frameon=False, labelcolor='black')
ax.tick_params(colors='black', labelsize=7.5)
ax.yaxis.label.set_color('black')
fig.savefig(FIG + "/Figure_S6.png", dpi=400, bbox_inches='tight')
fig.savefig(FIG + "/Figure_S6.pdf", bbox_inches='tight')
plt.close(fig)
print("新建 Supplementary Figure 6 (方差分解)")

# ============ 2. 按新编号重排补充图 ============
# 新编号 <- 旧文件
FIGSRC = {1: OLDFIG + "/Figure_S6", 2: OLDFIG + "/Figure_S1", 3: OLDFIG + "/Figure_S3",
          4: OLDFIG + "/Figure_S2", 5: OLDFIG + "/Figure_S4", 6: FIG + "/Figure_S6",
          7: OLDFIG + "/Figure_S5"}
tmp = {}
for new, src in FIGSRC.items():
    for ext in ["png", "pdf"]:
        if os.path.exists(src + "." + ext):
            tmp[(new, ext)] = src + "." + ext
staged = {}
for (new, ext), src in tmp.items():
    dst = FIG + "/_stage_S%d.%s" % (new, ext)
    shutil.copy(src, dst)
    staged[(new, ext)] = dst
for n in range(1, 8):
    for ext in ["png", "pdf"]:
        old = FIG + "/Figure_S%d.%s" % (n, ext)
        if os.path.exists(old):
            os.remove(old)
for (new, ext), src in staged.items():
    os.rename(src, FIG + "/Figure_S%d.%s" % (new, ext))
# 尺寸检查与必要缩放
MAXW, MAXH = 170.0, 210.0
for n in range(1, 8):
    p = FIG + "/Figure_S%d.png" % n
    im = Image.open(p)
    dpi = im.info.get('dpi', (400, 400))[0]
    w, h = im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4
    k = max(w / MAXW, h / MAXH, 1.0)
    if k > 1.0:
        im.save(p, dpi=(dpi * k, dpi * k))
    print("  Supplementary Figure %d: %.0f x %.0f mm" % (n, w / k, h / k))

# ============ 3. 补充表内容 ============
def rd(f):
    return pd.read_csv(BASE + "/" + f)


def fmt(df, cols=None, r=3):
    d = df[cols] if cols else df.copy()
    return d.round(r)


TABLES = []


def T(title, caption, frame):
    TABLES.append((title, caption, frame))


T("Supplementary Table 1. MRI acquisition parameters by manufacturer",
  "Acquisition parameters of the T1-weighted sequences contributing to the analysis, summarized by "
  "scanner manufacturer.",
  pd.read_csv(BASE + "/table_acquisition_params.csv"))

T("Supplementary Table 2. Patients versus controls in an age-matched subsample",
  "Each of the 261 healthy controls was matched 1:1 without replacement to the nearest-age patient "
  "within a one-year caliper, giving 246 pairs of identical mean age (63.7 versus 63.7 years, "
  "P = 0.997). Linear mixed-effects models with scanner as a random effect, adjusted for age, sex "
  "and education. Beta is the group difference in brain-PAD in years (patients minus controls).",
  fmt(rd("rev_agematched_pdvshc.csv"), ["Network", "beta", "lo", "hi", "p", "fdr", "n"], 4))

T("Supplementary Table 3. Age-prediction accuracy of the transferred model, by network",
  "Accuracy in the 261 healthy controls. MAE raw is the mean absolute error of predicted against "
  "chronological age; MAE bias-corrected is the mean absolute residual after the linear age-bias "
  "correction described in the Methods, and is the dispersion of the brain-PAD measure used in all "
  "analyses. Slope is the regression of predicted on chronological age; unity would indicate no "
  "compression.",
  fmt(rd("rev_B1_network_MAE.csv"), ["Network", "MAE_raw", "MAE_biascorr", "MAEc_lo", "MAEc_hi",
                                     "r", "RMSE_raw", "mean_offset", "slope"], 3))

_b4u = rd("rev_B4b_updrs3_score_bigscanners.csv").add_suffix("_UPDRS")
_b4m = rd("rev_B4b_moca_bigscanners.csv").add_suffix("_MoCA")
T("Supplementary Table 4. Sensitivity analysis restricted to scanner models with at least ten participants",
  "Analyses repeated in the 534 of 542 patients scanned on the 15 scanner models contributing at "
  "least ten participants each, excluding the 11 models with fewer (37 participants, including the "
  "single Toshiba model with four).",
  pd.concat([_b4u[["Network_UPDRS", "beta_UPDRS", "p_UPDRS", "fdr_UPDRS"]].rename(
      columns={"Network_UPDRS": "Network"}),
      _b4m[["beta_MoCA", "p_MoCA", "fdr_MoCA"]]], axis=1).round(4))

_a2u = rd("rev_A2_updrs3_score_rank.csv")
_a2m = rd("rev_A2_moca_rank.csv")
T("Supplementary Table 5. Robustness of the associations to the scale of brain-PAD",
  "Primary analyses repeated with brain-PAD replaced by its rank-based inverse-normal transform "
  "(scale-invariant), by Huber robust regression, and by partial Spearman correlation on covariate "
  "residuals. Motor associations were recovered identically; among the cognitive associations only "
  "the frontoparietal network was robust to all three.",
  pd.concat([_a2u[["Network", "beta_LME", "fdr_LME", "fdr_robust", "fdr_spearman"]].rename(
      columns={"beta_LME": "UPDRS_beta", "fdr_LME": "UPDRS_fdr_rank",
               "fdr_robust": "UPDRS_fdr_robust", "fdr_spearman": "UPDRS_fdr_spearman"}),
      _a2m[["beta_LME", "fdr_LME", "fdr_robust", "fdr_spearman"]].rename(
          columns={"beta_LME": "MoCA_beta", "fdr_LME": "MoCA_fdr_rank",
                   "fdr_robust": "MoCA_fdr_robust", "fdr_spearman": "MoCA_fdr_spearman"})],
      axis=1).round(4))

T("Supplementary Table 6. Regional brain-PAD and motor severity (MDS-UPDRS III, off medication)",
  "Linear mixed-effects models with scanner as a random effect, adjusted for age, sex and education. "
  "Beta is the change in MDS-UPDRS III per standard deviation of brain-PAD.",
  fmt(rd("results_final_updrs.csv"), None, 4))

T("Supplementary Table 7. Regional brain-PAD and global cognition (MoCA)",
  "Model specification as in Supplementary Table 6. Beta is the change in MoCA per standard "
  "deviation of brain-PAD.",
  fmt(rd("results_final_moca.csv"), None, 4))

T("Supplementary Table 8. Cognitive associations adjusted for motor severity",
  "MoCA models repeated with MDS-UPDRS III added as a covariate, testing whether residual motor "
  "impairment affecting the visuomotor items of the MoCA could account for the association.",
  fmt(rd("rev_moca_adj_updrs3.csv"), ["Network", "beta", "lo", "hi", "p", "fdr", "n"], 4))

T("Supplementary Table 9. Cognitive associations within the healthy control group",
  "The MoCA analysis repeated within the 261 healthy controls, and the group-by-brain-PAD "
  "interaction fitted across the full sample. No control-group association and no interaction "
  "survived correction.",
  fmt(rd("rev_B3_HC_association.csv"), None, 4))

T("Supplementary Table 10. Regional brain-PAD in patients versus controls",
  "Linear mixed-effects models with scanner as a random effect, adjusted for age, sex and "
  "education, in all 803 participants. Beta is the group difference in brain-PAD in years.",
  fmt(rd("results_final_pdvshc.csv"), None, 4))

T("Supplementary Table 11. Variance in brain-PAD attributable to scanner and to field strength",
  "Scanner-attributable variance is the intraclass correlation from a mixed model with scanner as a "
  "random intercept. Field-strength partial R squared contrasts the 754 participants scanned at 3 T "
  "with the 49 at 1.5 T. Values are given before and after ComBat harmonization.",
  fmt(rd("rev_B4_scanner_variance.csv"), None, 3))

T("Supplementary Table 12. Sensitivity analysis adjusted for disease duration and levodopa equivalent daily dose",
  "Primary associations repeated with disease duration and levodopa equivalent daily dose added as "
  "covariates.",
  pd.concat([rd("sens_updrs3_score_durLEDD.csv"), rd("sens_moca_durLEDD.csv")],
            axis=0, ignore_index=True).round(4))

T("Supplementary Table 13. Regional brain age carries information beyond cortical thickness",
  "For each significant association, the brain-PAD coefficient before and after adding the same "
  "network's raw mean cortical thickness to the model, with the P value of the thickness term "
  "itself.",
  rd("DAT_incremental_validity.csv").round(4) if os.path.exists(BASE + "/DAT_incremental_validity.csv")
  else rd("G1_incremental_validity.csv").round(4))

T("Supplementary Table 14. Grey-to-white matter intensity ratio before and after harmonization",
  "The association between network-mean grey-to-white matter intensity ratio and motor severity, "
  "estimated by ordinary least squares, by a mixed model with scanner as a random effect, and after "
  "ComBat harmonization of the intensity ratio itself. The final columns give the brain-PAD and "
  "harmonized intensity-ratio coefficients in the full model containing brain-PAD, raw cortical "
  "thickness and harmonized intensity ratio.",
  fmt(rd("rev_B2c_updrs3_score_gwr_site.csv"),
      ["Network", "b_GWR_OLS", "P_OLS", "b_GWR_LME", "P_LME", "b_GWRh_LME", "P_GWRh",
       "b_PAD_full", "P_PAD_full", "b_GWRh_full", "P_GWRh_full", "n"], 4))

_c = rd("rev_B2_correlations.csv")
_v = rd("rev_B2_VIF.csv")
T("Supplementary Table 15. Correlations among brain-PAD and the model input features, and collinearity statistics",
  "Pearson correlations in the 542 patients between each network's brain-PAD and that network's mean "
  "cortical thickness and mean grey-to-white matter intensity ratio, and between the two inputs. "
  "Variance inflation factors are from the model containing brain-PAD, thickness and intensity ratio "
  "together with the covariates; all are well below the conventional threshold of 5.",
  pd.concat([_c[["Network", "r_PAD_thickness", "r_PAD_GWR", "r_thickness_GWR"]],
             _v[["VIF_brainPAD", "VIF_thickness", "VIF_GWR"]]], axis=1).round(3))

T("Supplementary Table 16. Regional brain age versus dopaminergic denervation",
  "Association between each network's brain-PAD and putaminal striatal binding ratio in the 533 "
  "patients with dopamine-transporter imaging.",
  rd("DAT_coupling.csv").round(4))

_saa = rd("rev_B7_SAA_pos.csv").assign(Marker="alpha-synuclein seeding (odds ratio)")
_ab = rd("rev_B7_Abeta_ratio.csv").assign(Marker="CSF Abeta42/Abeta40")
_pt = rd("rev_B7_IU_pTau181_CSF.csv").assign(Marker="CSF p-tau181")
_pa = rd("rev_B7_pTau_Abeta42.csv").assign(Marker="CSF p-tau181/Abeta42")
_bio = pd.concat([_saa.rename(columns={"OR_per_SD": "estimate"}),
                  _ab.rename(columns={"beta": "estimate"}),
                  _pt.rename(columns={"beta": "estimate"}),
                  _pa.rename(columns={"beta": "estimate"})], ignore_index=True)
T("Supplementary Table 17. Regional brain age and co-pathology markers",
  "Association between each network's brain-PAD and cerebrospinal fluid alpha-synuclein seed "
  "amplification status (logistic models, odds ratio per standard deviation) and the amyloid and tau "
  "panel (linear mixed-effects models). Each participant's earliest available measurement was used. "
  "No association survived correction for any marker.",
  _bio[["Marker", "Network", "estimate", "p", "fdr", "n"]].round(4))

T("Supplementary Table 18. Regional brain age and extended clinical outcomes",
  "Linear mixed-effects models with scanner as a random effect, adjusted for age, sex and education, "
  "for outcomes spanning non-motor experiences, motor experiences of daily living, functional "
  "independence, disease stage and domain-specific cognition. Networks listed are those surviving "
  "false discovery rate correction within each outcome.",
  rd("rev_B5_summary.csv"))

_sub = []
for k, nm in [("tremor", "Tremor"), ("rigidity", "Rigidity"),
              ("bradykinesia", "Bradykinesia"), ("axial_pigd", "Axial / PIGD")]:
    d = rd("rev_B8_%s.csv" % k).assign(Subdomain=nm)
    _sub.append(d[["Subdomain", "Network", "beta", "p", "fdr", "n"]])
T("Supplementary Table 19. Regional brain age and MDS-UPDRS III motor subdomains",
  "Subdomain scores derived from item-level MDS-UPDRS Part III in the off or untreated state at "
  "baseline (n = 529), using standard item groupings. Only rigidity showed network-specific "
  "associations.",
  pd.concat(_sub, ignore_index=True).round(4))

T("Supplementary Table 20. Motor phenotype and age at onset",
  "Difference in regional brain-PAD between tremor-dominant (n = 383) and postural-instability "
  "gait-difficulty (n = 145) patients, and the interaction between brain-PAD and age at onset on "
  "motor severity. Neither was significant in any network.",
  fmt(rd("rev_B6_subtype_ageonset.csv"), None, 4))

_lu = rd("rev_long_updrs3_score.csv").assign(Outcome="MDS-UPDRS III")
_lm = rd("rev_long_moca.csv").assign(Outcome="MoCA")
_l2 = rd("rev_long_updrs2_score.csv").assign(Outcome="MDS-UPDRS II")
T("Supplementary Table 21. Baseline regional brain age and subsequent clinical change",
  "Time-by-brain-PAD interactions from linear mixed models with random intercepts and slopes for "
  "participant, adjusted for baseline age, sex and education, in patients with at least three "
  "assessments. Analyses were rebuilt for this revision on the full curated release, raising the "
  "sample from 193-195 to 349-464 depending on outcome. Progression to Hoehn and Yahr stage 3 "
  "(Cox models, n = 458, 24 events) and conversion to mild cognitive impairment or dementia "
  "(n = 464, 97 conversions) were likewise null and are reported in the Results.",
  pd.concat([_lu, _lm, _l2], ignore_index=True)[
      ["Outcome", "Network", "beta_time_x_PAD", "lo", "hi", "p", "fdr", "n_subj"]].round(4))

print("\n补充表共 %d 张" % len(TABLES))

# ============ 4. 补充图注 ============
LEGENDS = [
 ("Supplementary Figure 1. Participant flow.",
  "Number of PPMI participants with a baseline T1-weighted scan submitted for cortical processing, "
  "the number successfully processed with valid cortical thickness and grey-to-white matter "
  "features, and the number excluded at each stage."),
 ("Supplementary Figure 2. Validation of the reconstructed feature-extraction pipeline.",
  "Distributions of the vertex-wise features produced by the reconstructed pipeline against the "
  "input distributions expected by the transferred model."),
 ("Supplementary Figure 3. Random-forest importance of the ten networks.",
  "Importance of each network for classifying patients with poor versus preserved motor outcome, "
  "over 1,000 bootstrap resamples. Bars show the mean and whiskers the interquartile range of the "
  "bootstrap distribution; n = 542 patients."),
 ("Supplementary Figure 4. Standardized associations across networks and outcomes.",
  "Heat map of standardized regression coefficients dissociating the motor (MDS-UPDRS III, n = 542) "
  "from the cognitive (MoCA, n = 539) contribution of each network. Asterisks mark false discovery "
  "rate corrected P < 0.05. The colour bar encodes the standardized coefficient in clinical points "
  "per standard deviation of brain-PAD, with zero at the midpoint."),
 ("Supplementary Figure 5. Sensitivity of the main associations to scanner harmonization.",
  "Linear mixed-effects estimates with scanner as a random effect compared with ComBat-harmonized "
  "estimates, for each network and both outcomes. Points are standardized regression coefficients "
  "and error bars 95% confidence intervals; n = 542 (MDS-UPDRS III) and n = 539 (MoCA)."),
 ("Supplementary Figure 6. Variance in regional brain age attributable to acquisition.",
  "(A) Percentage of brain-PAD variance attributable to scanner, estimated as the intraclass "
  "correlation from a mixed model with scanner as a random intercept, before and after ComBat "
  "harmonization. (B) Partial R squared for field strength, contrasting the 754 participants "
  "scanned at 3 T with the 49 scanned at 1.5 T, before and after harmonization; asterisks mark "
  "P < 0.05. n = 803 participants."),
 ("Supplementary Figure 7. Baseline regional brain age and the rate of subsequent clinical change.",
  "Time-by-brain-PAD interaction for each network, testing whether baseline regional brain age "
  "predicts the rate of change in (A) MDS-UPDRS III (n = 349) and (B) MoCA (n = 440) over a median "
  "three-year follow-up. Points are interaction coefficients from linear mixed models with random "
  "intercepts and slopes for participant, and error bars 95% confidence intervals. No network "
  "survived false discovery rate correction."),
]

# ============ 5. Word 文档 ============
doc = docx.Document()
st = doc.styles['Normal']
st.font.name = 'Arial'
st.font.size = Pt(9)
h = doc.add_heading("Supplementary Material", 0)
p = doc.add_paragraph("Regional cortical brain age is independent of dopaminergic denervation in "
                      "Parkinson's disease")
p.runs[0].italic = True
doc.add_paragraph("Hu X, Zhu S, Qian X, Li W.   Manuscript BRAINCOM-2026-949, revision 1.")
doc.add_paragraph("Supplementary figures are supplied together with their legends as a single PDF, "
                  "as specified in the journal's formatting checklist. Their legends are reproduced "
                  "below for reference.")

doc.add_heading("Supplementary Figure legends", 1)
for t, c in LEGENDS:
    q = doc.add_paragraph()
    q.add_run(t).bold = True
    q.add_run(" " + c)

doc.add_heading("Supplementary Tables", 1)
for title, cap, frame in TABLES:
    doc.add_heading(title, 2)
    q = doc.add_paragraph()
    q.add_run(cap).italic = True
    frame = frame.fillna("")
    tb = doc.add_table(rows=1, cols=len(frame.columns))
    tb.style = 'Light Grid Accent 1'
    for j, col in enumerate(frame.columns):
        cell = tb.rows[0].cells[j]
        cell.text = ""
        r = cell.paragraphs[0].add_run(str(col))
        r.bold = True
        r.font.size = Pt(7.5)
    for _, row in frame.iterrows():
        cells = tb.add_row().cells
        for j, v in enumerate(row):
            cells[j].text = ""
            rr = cells[j].paragraphs[0].add_run(str(v))
            rr.font.size = Pt(7.5)

doc.add_heading("References cited in the Supplementary Material", 1)
doc.add_paragraph("No references are cited in the supplementary material beyond those in the main "
                  "reference list.")
doc.save(R1 + "/06_Supplementary_Material.docx")
print("补充材料 Word 已保存")

# ============ 6. 补充图单个 PDF ============
pdf_path = R1 + "/07_Supplementary_Figures.pdf"
with PdfPages(pdf_path) as pdf:
    for i, (title, cap) in enumerate(LEGENDS, 1):
        img = Image.open(FIG + "/Figure_S%d.png" % i)
        dpi = img.info.get('dpi', (400, 400))[0]
        w_mm, h_mm = img.size[0] / dpi * 25.4, img.size[1] / dpi * 25.4
        f = plt.figure(figsize=(8.27, 11.69))          # A4
        axi = f.add_axes([0.08, 0.34, 0.84, 0.60])
        axi.imshow(img)
        axi.axis('off')
        sc = min(0.84 * 210 / w_mm, 0.60 * 297 / h_mm, 1.0)
        axi.set_position([0.5 - 0.84 * sc * w_mm / 210 / 2, 0.94 - 0.60 * sc * h_mm / 297,
                          0.84 * sc * w_mm / 210, 0.60 * sc * h_mm / 297])
        f.text(0.08, 0.26, title, fontsize=10, fontweight='bold', va='top', wrap=True)
        f.text(0.08, 0.235, cap, fontsize=9, va='top', wrap=True)
        pdf.savefig(f)
        plt.close(f)
print("补充图 PDF 已保存 ->", pdf_path)
