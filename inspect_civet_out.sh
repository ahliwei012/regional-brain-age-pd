#!/bin/bash
# 挂载ext4并核验某例CIVET输出, 把关键文件拷到F盘持久保存
# 用法: bash inspect_civet_out.sh <密码> <ID>
PW="$1"; ID="$2"
IMG=/mnt/g/civet_ext4.img; WORK=/mnt/civet
echo "$PW" | sudo -S -p '' mkdir -p "$WORK"
mountpoint -q "$WORK" || echo "$PW" | sudo -S -p '' mount -o loop "$IMG" "$WORK"
echo "$PW" | sudo -S -p '' chown -R "$USER":"$USER" "$WORK" 2>/dev/null
D="$WORK/out/$ID"

echo "=== 输出顶层目录 ==="
ls "$D" 2>/dev/null
echo
echo "=== thickness/ ==="
ls -lh "$D/thickness/" 2>/dev/null
echo
echo "=== 重采样厚度txt (模型厚度输入源, 行数应=40962/半球) ==="
for s in left right; do
  f=$(find "$D" -iname "*rsl*tlink*${s}*.txt" 2>/dev/null | head -1)
  [ -n "$f" ] && echo "$(wc -l < "$f") 行  $f"
done
echo
echo "=== 拷贝关键输出到 F 盘持久保存 ==="
OUT_F="/mnt/f/PD_brainage_data/civet_out/$ID"
mkdir -p "$OUT_F"
cp -f "$D"/thickness/*.txt "$OUT_F"/ 2>/dev/null
cp -f "$D"/surfaces/*rsl*.obj "$OUT_F"/ 2>/dev/null
cp -f "$D"/final/*t1_final.mnc "$OUT_F"/ 2>/dev/null
cp -rf "$D"/verify "$OUT_F"/ 2>/dev/null
echo "已拷到 $OUT_F :"; ls "$OUT_F" | head -20
du -sh "$D" 2>/dev/null | awk '{print "单例输出总大小: "$1}'
