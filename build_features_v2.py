# -*- coding: utf-8 -*-
"""改进版: log(GM/WM)在40962细表面(edges80k)上平滑后再降采样到10242; 32核并行。
用法: python build_features_v2.py <idlist文件>"""
import numpy as np, nibabel as nib, os, sys
import scipy.sparse as sp
from scipy.ndimage import map_coordinates
from multiprocessing import Pool

CIVET = "F:/PD_brainage_data/civet_out"
SURF  = "F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main/surface_information"
OUTF  = "F:/PD_brainage_data/features_v2"
os.makedirs(OUTF, exist_ok=True)
NV, NF = 10242, 40962          # 每半球: 模型10242, CIVET 40962
TARGET_STD = 0.064

# ---- 40962x2 图邻接 (edges80k, 1-based->0) ----
_e = np.loadtxt(SURF + "/edges80k.txt", dtype=int) - 1
N80 = 2 * NF
A80 = sp.coo_matrix((np.ones(len(_e)*2), (np.r_[_e[:,0],_e[:,1]], np.r_[_e[:,1],_e[:,0]])),
                    shape=(N80, N80)).tocsr()
deg80 = np.asarray(A80.sum(1)).ravel(); deg80[deg80 == 0] = 1

def smooth80(x, K):
    for _ in range(K):
        x = 0.5*x + 0.5*(A80.dot(x)/deg80)
    return x

def read_obj_points(path):
    with open(path) as f:
        first = f.readline().split(); n = int(first[-1]); rest = f.read().split()
    return np.array(rest[:3*n], dtype=float).reshape(n, 3)

def gmwm_full(d, ID, hemi, inv, data):
    w = read_obj_points(f"{d}/PPMI_{ID}_white_surface_rsl_{hemi}_81920.obj")
    g = read_obj_points(f"{d}/PPMI_{ID}_gray_surface_rsl_{hemi}_81920.obj")
    dirv = g - w; L = np.linalg.norm(dirv, axis=1, keepdims=True); L[L==0] = 1; nrm = dirv/L
    def samp(pts):
        return map_coordinates(data, nib.affines.apply_affine(inv, pts).T, order=1, mode="nearest")
    igm = np.clip(samp(w + nrm), 1, None)
    iwm = np.clip(samp(w - nrm), 1, None)
    return np.log(igm/iwm)                     # 完整 40962

def raw_gw81924(ID):
    d = f"{CIVET}/{ID}"
    img = nib.load(f"{d}/PPMI_{ID}_t1_final.mnc"); data = img.get_fdata(); inv = np.linalg.inv(img.affine)
    return np.concatenate([gmwm_full(d, ID, "left", inv, data),
                           gmwm_full(d, ID, "right", inv, data)])   # 81924

def process(args):
    ID, K = args
    d = f"{CIVET}/{ID}"
    thL = np.loadtxt(f"{d}/PPMI_{ID}_native_rms_rsl_tlink_20mm_left.txt")[:NV]
    thR = np.loadtxt(f"{d}/PPMI_{ID}_native_rms_rsl_tlink_20mm_right.txt")[:NV]
    th = np.concatenate([thL, thR])
    gw_full = smooth80(raw_gw81924(ID), K)               # 40962上平滑
    gw = np.concatenate([gw_full[:NV], gw_full[NF:NF+NV]])  # 降采样: 左前10242 + 右前10242
    feat = np.column_stack([th, gw])
    np.savetxt(f"{OUTF}/{ID}_features_20k.txt", feat, fmt="%.5f", delimiter="\t")
    return f"{ID}: thick={th.mean():.3f} gw_mean={gw.mean():+.4f} gw_std={gw.std():.4f}"

if __name__ == "__main__":
    idlist = sys.argv[1] if len(sys.argv) > 1 else None
    IDS = [l.strip() for l in open(idlist)] if idlist else []  # insert your own PPMI subject IDs (PATNO) here
    IDS = [i for i in IDS if i]
    bestK = 20   # 固定平滑迭代数(=FWHM, 在HC上校准得到); PD/HC统一, 避免组间混杂
    print(f">>> 固定 K={bestK} (PD/HC一致, 共{len(IDS)}例)")
    nproc = min(32, len(IDS))
    print(f"并行 {nproc} 进程处理 {len(IDS)} 例...")
    with Pool(nproc) as p:
        for r in p.imap_unordered(process, [(i, bestK) for i in IDS]):
            print(" ", r)
    print("完成 ->", OUTF)
