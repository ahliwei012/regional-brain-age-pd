#!/bin/bash
LOG=/mnt/f/PD_brainage_data/logs/civet_hc_batch.log
echo "=== 批次首尾 ==="
head -1 "$LOG"; tail -2 "$LOG"
echo "=== 统计 ==="
echo "ok=$(grep -c '\[ok\]' "$LOG")  skip=$(grep -c '\[skip\]' "$LOG")  FAIL=$(grep -c '\[FAIL\]' "$LOG")  err=$(grep -c '\[err\]' "$LOG")"
echo "=== 失败/错误明细(若有) ==="
grep -E '\[FAIL\]|\[err\]' "$LOG" | head
echo "=== HC在civet_out中且有厚度文件的例数 ==="
n=0
while read id; do
  [ -f "/mnt/f/PD_brainage_data/civet_out/$id/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt" ] && n=$((n+1))
done < /mnt/f/PD_brainage_data/hc_idlist.txt
echo "HC完成(有厚度文件): $n / $(wc -l < /mnt/f/PD_brainage_data/hc_idlist.txt)"
echo "=== civet_out总目录数(PD+HC) ==="
ls -d /mnt/f/PD_brainage_data/civet_out/*/ 2>/dev/null | wc -l
