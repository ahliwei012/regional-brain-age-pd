# -*- coding: utf-8 -*-
"""仿参考文献(Lancet Digit Health)风格出初步结果图。"""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")

BASE = r"F:/PD_brainage_data"; FIG = f"{BASE}/figures"; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                     "axes.spines.right":False,"figure.dpi":300,"savefig.dpi":300})
REG = ["sns","fpn","drs","vnt","dft","slt","lng","adt","vsl","lmb"]
EN  = {"sns":"Sensorimotor","fpn":"Frontoparietal","drs":"Dorsal attention","vnt":"Ventral attention",
       "dft":"Default mode","slt":"Salience","lng":"Language","adt":"Auditory","vsl":"Visual","lmb":"Limbic"}
def star(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""

df = pd.read_csv(f"{BASE}/analysis_merged.csv"); pdf = df[df.group=="PD"]
def assoc(outcome):
    rows=[]
    for c in REG:
        d=pdf.dropna(subset=[c+"_PAD",outcome,"age_at_visit","SEX","EDUCYRS"]).copy()
        d["z"]=(d[c+"_PAD"]-d[c+"_PAD"].mean())/d[c+"_PAD"].std()
        m=smf.ols(f"{outcome} ~ z + age_at_visit + C(SEX) + EDUCYRS",d).fit()
        ci=m.conf_int().loc["z"]
        rows.append([c,EN[c],m.params["z"],ci[0],ci[1],m.pvalues["z"]])
    r=pd.DataFrame(rows,columns=["r","name","beta","lo","hi","p"]); r["fdr"]=multipletests(r["p"],method="fdr_bh")[1]
    return r
U=assoc("updrs3_score"); M=assoc("moca")

# ---------- Fig1: 森林图 区域脑龄PAD -> UPDRS-III ----------
fig,ax=plt.subplots(figsize=(6,4.2))
s=U.sort_values("beta").reset_index(drop=True)
y=np.arange(len(s))
for i in range(len(s)):
    c="#c0392b" if s["fdr"].iloc[i]<0.05 else "#95a5a6"
    ax.plot([s["lo"].iloc[i],s["hi"].iloc[i]],[i,i],color=c,lw=2,zorder=1)
    ax.scatter(s["beta"].iloc[i],i,color=c,s=45,zorder=2)
    if star(s["fdr"].iloc[i]): ax.text(s["hi"].iloc[i]+0.02,i,star(s["fdr"].iloc[i]),va="center",fontsize=11,color="#c0392b")
ax.axvline(0,color="k",lw=0.8,ls="--")
ax.set_yticks(y); ax.set_yticklabels(s["name"]); ax.set_xlabel("Standardized β (regional brain-PAD → UPDRS-III)")
ax.set_title("Regional cortical brain age vs motor severity in PD\n(adjusted age/sex/education; red = FDR p<0.05)",fontsize=10)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig1_forest_updrs.png"); plt.close()

# ---------- Fig2: 热图 区域 × {UPDRS-III, MoCA} ----------
mat=np.column_stack([U["beta"].values,M["beta"].values])
fig,ax=plt.subplots(figsize=(4.2,5))
norm=TwoSlopeNorm(vmin=min(mat.min(),-abs(mat).max()),vcenter=0,vmax=abs(mat).max())
im=ax.imshow(mat,cmap="RdBu_r",norm=norm,aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["UPDRS-III\n(motor)","MoCA\n(cognition)"])
ax.set_yticks(range(len(REG))); ax.set_yticklabels([EN[c] for c in REG])
for i in range(len(REG)):
    ax.text(0,i,star(U["fdr"].iloc[i]),ha="center",va="center",fontsize=12,color="k")
    ax.text(1,i,star(M["fdr"].iloc[i]),ha="center",va="center",fontsize=12,color="k")
ax.set_title("Regional brain-PAD ↔ clinical outcomes\n(standardized β; * FDR<0.05)",fontsize=10)
plt.colorbar(im,ax=ax,fraction=0.046,pad=0.04,label="Standardized β")
plt.tight_layout(); plt.savefig(f"{FIG}/Fig2_heatmap.png"); plt.close()

# ---------- Fig3: 随机森林特征重要性 (仿Fig4) ----------
rf=pd.read_csv(f"{BASE}/results_rf_importance.csv")  # 中文名, 换英文
zh2en={"感觉运动":"Sensorimotor","额顶":"Frontoparietal","背侧注意":"Dorsal attention","腹侧注意":"Ventral attention",
       "默认模式":"Default mode","突显":"Salience","语言":"Language","听觉":"Auditory","视觉":"Visual","边缘":"Limbic"}
rf["en"]=rf["region"].map(zh2en); rf=rf.sort_values("importance",ascending=True)
fig,ax=plt.subplots(figsize=(6,4.2))
ax.barh(rf["en"],rf["importance"],xerr=rf["sd"],color="#aed6f1",edgecolor="#2980b9",
        error_kw=dict(ecolor="#34495e",capsize=2))
ax.set_xlabel("Mean feature importance (Gini, 1000 bootstraps)")
ax.set_title("Regional brain-PAD importance for predicting poor motor outcome",fontsize=10)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig3_rf_importance.png"); plt.close()

# ---------- Fig4: 散点 顶区(视觉/腹侧注意)PAD vs UPDRS-III ----------
fig,axes=plt.subplots(1,2,figsize=(9,4))
for ax,c in zip(axes,["vsl","vnt"]):
    d=pdf.dropna(subset=[c+"_PAD","updrs3_score"])
    ax.scatter(d[c+"_PAD"],d["updrs3_score"],s=12,alpha=0.4,color="#2980b9")
    z=np.polyfit(d[c+"_PAD"],d["updrs3_score"],1); xx=np.linspace(d[c+"_PAD"].min(),d[c+"_PAD"].max(),50)
    ax.plot(xx,np.polyval(z,xx),color="#c0392b",lw=2)
    rr=np.corrcoef(d[c+"_PAD"],d["updrs3_score"])[0,1]
    ax.set_xlabel(f"{EN[c]} brain-PAD (years)"); ax.set_ylabel("UPDRS-III"); ax.set_title(f"r={rr:.2f}, n={len(d)}")
fig.suptitle("Regional brain age vs motor severity (PD)",fontsize=11)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig4_scatter.png"); plt.close()

print("图已生成:")
for f in ["Fig1_forest_updrs","Fig2_heatmap","Fig3_rf_importance","Fig4_scatter"]:
    print("  ", f"{FIG}/{f}.png")
print("\nUPDRS-III关联(FDR<0.05):", U[U.fdr<0.05]["name"].tolist())
print("MoCA关联(FDR<0.05):", M[M.fdr<0.05]["name"].tolist())
