# -*- coding: utf-8 -*-
"""构建模型输入特征(20484x2): col1=厚度(CIVET 20mm平滑,取前10242), col2=log(GM/WM)+表面平滑。
输出到 F:/PD_brainage_data/features/<ID>_features_20k.txt  (绝不写C盘)"""
import numpy as np, nibabel as nib, os
import scipy.sparse as sp
from scipy.ndimage import map_coordinates

CIVET = "F:/PD_brainage_data/civet_out"
SURF  = "F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main/surface_information"
OUTF  = "F:/PD_brainage_data/features"
os.makedirs(OUTF, exist_ok=True)
import sys
if len(sys.argv) > 1:
    IDS = [l.strip() for l in open(sys.argv[1]) if l.strip()]
else:
    IDS = []  # insert your own PPMI subject IDs (PATNO) here
NV  = 10242            # 每半球 ico5 顶点数
TARGET_STD = 0.064     # 样例文件第2列目标std

def read_obj_points(path):
    with open(path) as f:
        first = f.readline().split(); n = int(first[-1]); rest = f.read().split()
    return np.array(rest[:3*n], dtype=float).reshape(n, 3)

# ---- 20484顶点图邻接(edges20k, 0-based) ----
e = np.loadtxt(SURF + "/edges20k.txt", dtype=int) - 1   # 文件是1-based, 转0-based
N = 2 * NV
A = sp.coo_matrix((np.ones(len(e)*2),
                   (np.r_[e[:,0], e[:,1]], np.r_[e[:,1], e[:,0]])), shape=(N, N)).tocsr()
deg = np.asarray(A.sum(1)).ravel(); deg[deg == 0] = 1
def smooth(x, K):
    for _ in range(K):
        x = 0.5*x + 0.5*(A.dot(x)/deg)
    return x

def gmwm(ID, hemi, inv, data):
    d = f"{CIVET}/{ID}"
    w = read_obj_points(f"{d}/PPMI_{ID}_white_surface_rsl_{hemi}_81920.obj")
    g = read_obj_points(f"{d}/PPMI_{ID}_gray_surface_rsl_{hemi}_81920.obj")
    dirv = g - w; L = np.linalg.norm(dirv, axis=1, keepdims=True); L[L==0] = 1; nrm = dirv/L
    def samp(pts):
        return map_coordinates(data, nib.affines.apply_affine(inv, pts).T, order=1, mode="nearest")
    igm = np.clip(samp(w + nrm), 1, None)
    iwm = np.clip(samp(w - nrm), 1, None)
    return np.log(igm/iwm)[:NV]

# ---- 固定平滑迭代数 K=30 (与PD复刻一致; 之前在100005上校准 std≈0.062) ----
bestK = 30
print(f">>> 使用固定 K={bestK} (与PD一致), 共 {len(IDS)} 例")

# ---- 构建全部5例特征 ----
for ID in IDS:
    d = f"{CIVET}/{ID}"
    img = nib.load(f"{d}/PPMI_{ID}_t1_final.mnc"); data = img.get_fdata(); inv = np.linalg.inv(img.affine)
    thL = np.loadtxt(f"{d}/PPMI_{ID}_native_rms_rsl_tlink_20mm_left.txt")[:NV]
    thR = np.loadtxt(f"{d}/PPMI_{ID}_native_rms_rsl_tlink_20mm_right.txt")[:NV]
    th  = np.concatenate([thL, thR])
    gw  = np.concatenate([gmwm(ID,"left",inv,data), gmwm(ID,"right",inv,data)])
    gw  = smooth(gw, bestK)
    feat = np.column_stack([th, gw])
    out = f"{OUTF}/{ID}_features_20k.txt"
    np.savetxt(out, feat, fmt="%.5f", delimiter="\t")
    print(f"  {ID}: thick mean={th.mean():.3f}  gw mean={gw.mean():+.4f} std={gw.std():.4f} -> {out}")
print("\n完成. 特征文件在", OUTF)
