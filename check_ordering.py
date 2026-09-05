# -*- coding: utf-8 -*-
"""验证'取CIVET前10242'是否匹配模型edges20k的顶点顺序: 用CIVET前10242坐标+edges20k算边长。"""
import numpy as np
SURF = "F:/Deep_learning_derived_ MRI regional brain age/regional_Brain_age-main/surface_information"
W = "F:/PD_brainage_data/civet_out/100890/PPMI_100890_white_surface_rsl_left_81920.obj"

def read_obj_points(path):
    with open(path) as f:
        first = f.readline().split(); n = int(first[-1]); rest = f.read().split()
    return np.array(rest[:3*n], dtype=float).reshape(n, 3)

w = read_obj_points(W)             # 40962
print("CIVET左白质表面顶点数:", len(w))
w10 = w[:10242]                    # 取前10242

edges = np.loadtxt(SURF + "/edges20k.txt", dtype=int) - 1   # 0-based, 0..20483
left_e = edges[(edges[:,0] < 10242) & (edges[:,1] < 10242)]
d = np.linalg.norm(w10[left_e[:,0]] - w10[left_e[:,1]], axis=1)
print(f"左半球edges20k边数: {len(left_e)}")
print(f"用CIVET前10242坐标算的边长: mean={d.mean():.2f}mm median={np.median(d):.2f} max={d.max():.2f} >10mm占比={100*(d>10).mean():.1f}%")
print("(正确10242皮层表面相邻边长应~2-4mm; 若均值大/长边多 => 顺序不匹配)")

# 对比: 全40962表面相邻顶点的典型边长(用前1000个顶点到其最近邻)
from scipy.spatial import cKDTree
t = cKDTree(w); dd,_ = t.query(w[:2000], k=2)
print(f"\n参照: CIVET 40962表面最近邻间距 ~{dd[:,1].mean():.2f}mm")
