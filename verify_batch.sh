#!/bin/bash
DONE=/mnt/f/PD_brainage_data/civet_out
echo "=== 批次日志(完整) ==="
cat /mnt/f/PD_brainage_data/logs/civet_batch_test6.log
echo
echo "=== F盘已完成例 + 厚度文件行数校验(应40962/半球) ==="
for d in "$DONE"/*/; do
  id=$(basename "$d")
  lf="$d/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt"
  rf="$d/PPMI_${id}_native_rms_rsl_tlink_20mm_right.txt"
  ln=$( [ -f "$lf" ] && wc -l < "$lf" || echo NA )
  rn=$( [ -f "$rf" ] && wc -l < "$rf" || echo NA )
  nfiles=$(ls "$d" | wc -l)
  echo "  $id : L=$ln R=$rn  文件数=$nfiles"
done
echo
echo "=== 每例拷出文件清单(以100006为例) ==="
ls -lh "$DONE/100006/" 2>/dev/null
echo
echo "=== CIVET QC(civet_qc.txt 关键行, 100006) ==="
cat "$DONE/100006/PPMI_100006_civet_qc.txt" 2>/dev/null | head -20
echo
echo "=== ext4 scratch 是否已清空(out应只剩极少) ==="
ls /mnt/civet/out/ 2>/dev/null | head
df -h /mnt/civet 2>/dev/null | tail -1
