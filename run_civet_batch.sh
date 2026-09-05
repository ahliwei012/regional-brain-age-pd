#!/bin/bash
# 批量 CIVET: 自挂载ext4 + 并行 + 每例拷关键输出到F并清scratch + 断点续跑
# 用法: bash run_civet_batch.sh <密码> <并行数> <T1源根目录> <ID列表文件>
#   T1源: <RAWDIR>/<ID>/<ID>.nii
set -u
PW="$1"; NPAR="$2"; RAWDIR="$3"; IDLIST="$4"
SIF=/mnt/g/apptainer/civet_2.1.1.sif
IMG=/mnt/g/civet_ext4.img; WORK=/mnt/civet
SRC="$WORK/sourcedir"; OUT="$WORK/out"
DONE=/mnt/f/PD_brainage_data/civet_out
PREFIX=PPMI
export APPTAINER_TMPDIR=/mnt/g/apptainer/tmp

# 挂载 ext4
echo "$PW" | sudo -S -p '' mkdir -p "$WORK"
mountpoint -q "$WORK" || echo "$PW" | sudo -S -p '' mount -o loop "$IMG" "$WORK"
echo "$PW" | sudo -S -p '' chown -R "$USER":"$USER" "$WORK" 2>/dev/null
mkdir -p "$SRC" "$OUT" "$DONE"

export SIF SRC OUT DONE PREFIX RAWDIR
proc_one() {
  ID="$1"
  local key="$DONE/$ID/${PREFIX}_${ID}_native_rms_rsl_tlink_20mm_left.txt"
  [ -f "$key" ] && { echo "[skip] $ID"; return 0; }
  local nii tmpnii=""; nii=$(ls "$RAWDIR/$ID/"*.nii 2>/dev/null | head -1)
  if [ -z "$nii" ]; then
    local gz; gz=$(ls "$RAWDIR/$ID/"*.nii.gz 2>/dev/null | head -1)
    [ -z "$gz" ] && { echo "[err] $ID 无nii/nii.gz"; return 1; }
    tmpnii="$SRC/${PREFIX}_${ID}_t1.nii"; gunzip -c "$gz" > "$tmpnii"; nii="$tmpnii"
  fi
  local mnc="$SRC/${PREFIX}_${ID}_t1.mnc"
  [ -f "$mnc" ] || apptainer exec -B "$SRC":"$SRC" -B "$(dirname "$nii")":"$(dirname "$nii")" "$SIF" \
       nii2mnc "$nii" "$mnc" >/dev/null 2>&1
  [ -n "$tmpnii" ] && rm -f "$tmpnii"
  apptainer exec -B "$SRC":/data_in -B "$OUT":/data_out "$SIF" \
    CIVET_Processing_Pipeline -sourcedir /data_in -targetdir /data_out -prefix "$PREFIX" \
    -N3-distance 200 -lsq12 -thickness tlink 20 -resample-surfaces -spawn -run "$ID" \
    > "$OUT/${ID}.civetlog" 2>&1
  local d="$OUT/$ID"
  if [ -f "$d/thickness/${PREFIX}_${ID}_native_rms_rsl_tlink_20mm_left.txt" ]; then
    mkdir -p "$DONE/$ID"
    cp -f "$d"/thickness/*rsl*.txt "$DONE/$ID/" 2>/dev/null
    cp -f "$d"/surfaces/*_surface_rsl_*.obj "$DONE/$ID/" 2>/dev/null
    cp -f "$d"/final/*t1_final.mnc "$DONE/$ID/" 2>/dev/null
    cp -f "$d"/verify/*civet_qc.txt "$DONE/$ID/" 2>/dev/null
    cp -f "$OUT/${ID}.civetlog" "$DONE/$ID/" 2>/dev/null
    rm -rf "$d" "$mnc" "$OUT/${ID}.civetlog"
    echo "[ok] $ID"
  else
    echo "[FAIL] $ID 未产出厚度, 保留scratch待查"
  fi
}
export -f proc_one

echo "=== 批次开始 $(date)  并行=$NPAR  总数=$(wc -l < "$IDLIST") ==="
cat "$IDLIST" | xargs -P "$NPAR" -I {} bash -c 'proc_one "$@"' _ {}
echo "=== 批次结束 $(date)  已完成例数=$(ls "$DONE" | wc -l) ==="
