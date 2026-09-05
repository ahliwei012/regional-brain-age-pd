# -*- coding: utf-8 -*-
"""反推灰白质强度比特征: 沿白质表面法线内外1mm采样t1_final, 试多个公式, 对比样例文件第2列。"""
import numpy as np, nibabel as nib, glob, os
from scipy.ndimage import map_coordinates

CIVET = "F:/PD_brainage_data/civet_out"
TESTDIR = "F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main/test_dataset"
IDS = []  # insert your own PPMI subject IDs (PATNO) here

def read_obj_points(path):
    with open(path) as f:
        first = f.readline().split()
        n = int(first[-1])
        rest = f.read().split()
    pts = np.array(rest[:3*n], dtype=float).reshape(n, 3)
    return pts, n

# ---------- 目标: 样例文件第2列 ----------
print("=== 目标: 样例文件(civet_0X)第2列分布 ===")
tcols = []
for f in sorted(glob.glob(TESTDIR + "/civet_0*_features_20k.txt")):
    a = np.loadtxt(f)
    tcols.append(a[:, 1])
    print(f"  {os.path.basename(f)}: mean={a[:,1].mean():+.4f} std={a[:,1].std():.4f}")
tall = np.concatenate(tcols)
print(f"  >>> 合计 raw: mean={tall.mean():+.4f} std={tall.std():.4f}  abs_mean={np.abs(tall).mean():.4f}\n")

# ---------- 反推 ----------
formulas = {
    "GM/WM":           lambda gm, wm: gm/wm,
    "WM/GM":           lambda gm, wm: wm/gm,
    "(GM-WM)/WM":      lambda gm, wm: (gm-wm)/wm,
    "(GM-WM)/(GM+WM)": lambda gm, wm: (gm-wm)/(gm+wm),
    "log(GM/WM)":      lambda gm, wm: np.log(gm/wm),
}
agg = {k: [] for k in formulas}
wm_means, gm_means = [], []

for ID in IDS:
    d = f"{CIVET}/{ID}"
    img = nib.load(f"{d}/PPMI_{ID}_t1_final.mnc")
    data = img.get_fdata()
    inv = np.linalg.inv(img.affine)
    for hemi in ["left", "right"]:
        w, n = read_obj_points(f"{d}/PPMI_{ID}_white_surface_rsl_{hemi}_81920.obj")
        g, _ = read_obj_points(f"{d}/PPMI_{ID}_gray_surface_rsl_{hemi}_81920.obj")
        dirv = g - w
        L = np.linalg.norm(dirv, axis=1, keepdims=True); L[L == 0] = 1
        nrm = dirv / L
        wm_pts = w - 1.0*nrm   # 内1mm (白质侧)
        gm_pts = w + 1.0*nrm   # 外1mm (灰质侧)
        def samp(pts):
            vox = nib.affines.apply_affine(inv, pts)
            return map_coordinates(data, vox.T, order=1, mode="nearest")
        iwm = samp(wm_pts); igm = samp(gm_pts)
        ok = (iwm > 0) & (igm > 0) & np.isfinite(iwm) & np.isfinite(igm)
        wm_means.append(iwm[ok].mean()); gm_means.append(igm[ok].mean())
        for k, fn in formulas.items():
            agg[k].append(fn(igm[ok], iwm[ok]))

print(f"=== 对齐自检: 平均WM强度={np.mean(wm_means):.1f}  平均GM强度={np.mean(gm_means):.1f}  (WM应>GM) ===\n")
print(f"=== 各公式在5例(双侧40962)分布  [目标 raw mean={tall.mean():+.3f} std={tall.std():.3f}] ===")
best = None
for k in formulas:
    v = np.concatenate(agg[k]); v = v[np.isfinite(v)]
    diff = abs(v.mean() - tall.mean())
    print(f"  {k:18s}: mean={v.mean():+.4f} std={v.std():.4f}  |与目标均值差|={diff:.4f}")
    if best is None or diff < best[1]:
        best = (k, diff)
print(f"\n>>> 最接近目标均值的公式: {best[0]} (差 {best[1]:.4f})")
