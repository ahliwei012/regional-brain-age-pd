#!/bin/bash
# 单例 CIVET 测试 (自挂载 ext4 + 转换 + 跑 CIVET, 全在一个调用内)
# 用法: bash run_civet_test.sh <ID> <sudo密码>
set -e
ID="$1"
PW="$2"
SIF=/mnt/g/apptainer/civet_2.1.1.sif
IMG=/mnt/g/civet_ext4.img
WORK=/mnt/civet
SRC="$WORK/sourcedir"
OUT="$WORK/out"
STAGE=/mnt/f/PD_brainage_data/civet/sourcedir
PREFIX=PPMI

echo "[0] 确保 ext4 scratch 已挂载"
echo "$PW" | sudo -S -p '' mkdir -p "$WORK"
if ! mountpoint -q "$WORK"; then
  echo "$PW" | sudo -S -p '' mount -o loop "$IMG" "$WORK"
  echo "$PW" | sudo -S -p '' chown -R "$USER":"$USER" "$WORK"
fi
mkdir -p "$SRC" "$OUT"
echo "    挂载: $(mountpoint "$WORK")"

echo "[1] 拷贝 nii 到 scratch"
cp -f "$STAGE/${PREFIX}_${ID}_t1.nii" "$SRC/"

echo "[2] nii -> mnc"
if [ ! -f "$SRC/${PREFIX}_${ID}_t1.mnc" ]; then
  apptainer exec -B "$SRC":"$SRC" "$SIF" \
    nii2mnc "$SRC/${PREFIX}_${ID}_t1.nii" "$SRC/${PREFIX}_${ID}_t1.mnc"
fi
rm -f "$SRC/${PREFIX}_${ID}_t1.nii"
ls -lh "$SRC/${PREFIX}_${ID}_t1.mnc"

echo "[3] CIVET 处理 $ID  开始: $(date)"
apptainer exec -B "$SRC":/data_in -B "$OUT":/data_out "$SIF" \
  CIVET_Processing_Pipeline \
    -sourcedir /data_in -targetdir /data_out -prefix "$PREFIX" \
    -N3-distance 200 -lsq12 -thickness tlink 20 -resample-surfaces \
    -spawn -run "$ID"
rc=$?
echo "[CIVET] 结束 $ID rc=$rc  $(date)"
echo "[输出树]"
ls -R "$OUT/$PREFIX/$ID" 2>/dev/null | head -50
