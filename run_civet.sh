#!/bin/bash
# ============================================================================
# CIVET 处理脚本草稿 (待 SIF 就绪后核验)
# 用法: bash run_civet.sh <PATNO>          例: bash run_civet.sh 100005
# 前置: SIF = /mnt/g/apptainer/civet_2.1.1.sif
#       输入已按 <prefix>_<id>_t1.nii 暂存在 sourcedir
# ============================================================================
set -e
ID="$1"
SIF=/mnt/g/apptainer/civet_2.1.1.sif
SRC=/mnt/f/PD_brainage_data/civet/sourcedir     # 输入: PPMI_<ID>_t1.nii
OUT=/mnt/f/PD_brainage_data/civet/out           # 输出目录
PREFIX=PPMI
mkdir -p "$OUT"

# --- 选项说明 (★=需向作者/容器核验) ---
#  -resample-surfaces : ★必须, 才能得到顶点对应的标准表面(20484的前提)
#  -thickness tlink 30: ★方法=tlink(论文用内外表面欧氏距离), 平滑核待确认
#  -N3-distance       : ★3T 常用 200; 论文未明确, 待确认
#  -lsq12             : 12参数线性配准到 MNI
#  -surface-atlas     : 模型的 AAL 分区是后处理另做, CIVET 阶段可不加
# ----------------------------------------------------------------------------

echo "[CIVET] 处理 $PREFIX $ID  (开始: $(date))"
apptainer exec \
  -B "$SRC":/data_in \
  -B "$OUT":/data_out \
  "$SIF" \
  CIVET_Processing_Pipeline \
    -prefix "$PREFIX" \
    -sourcedir /data_in \
    -targetdir /data_out \
    -N3-distance 200 \
    -lsq12 \
    -resample-surfaces \
    -thickness tlink 30 \
    -no-surface-atlas \
    -run "$ID"
echo "[CIVET] 完成 $ID  (结束: $(date))"

# ============================================================================
# 待 SIF 就绪后要先做的核验 (不要直接跑全部):
#   1) apptainer exec $SIF CIVET_Processing_Pipeline -help   # 看真实选项名
#   2) apptainer exec $SIF cat /opt/CIVET*/VERSION           # 确认2.1.1
#   3) 单跑 1 例, 测时长 + 看 drvfs(/mnt) 是否报错/过慢
#      若 drvfs 太慢: 改用 WSL 原生 ext4 做 scratch, 完成后拷回 F/G
#   4) 确认输出里有: *_native_rms_rsl_tlink_*_left/right.txt (重采样厚度)
#      + 灰白质强度比仍需作者的特征提取脚本
# ============================================================================
