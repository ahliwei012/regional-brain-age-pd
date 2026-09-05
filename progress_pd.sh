#!/bin/bash
LOG=/mnt/f/PD_brainage_data/logs/civet_pd_all.log
DONE=/mnt/f/PD_brainage_data/civet_out
echo "=== 批次日志关键行 ==="
grep -E 'P1 \(362\)|P2 \(198\)|全量PD完成' "$LOG" 2>/dev/null
echo "ok=$(grep -c '\[ok\]' "$LOG" 2>/dev/null)  skip=$(grep -c '\[skip\]' "$LOG" 2>/dev/null)  FAIL=$(grep -c '\[FAIL\]' "$LOG" 2>/dev/null)  err=$(grep -c '\[err\]' "$LOG" 2>/dev/null)"
echo
echo "=== 已完成(有厚度文件)统计 ==="
pd=0; hc=0
while read id; do [ -f "$DONE/$id/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt" ] && pd=$((pd+1)); done < <(cat /mnt/f/PD_brainage_data/pd_p1_idlist.txt /mnt/f/PD_brainage_data/pd_p2_idlist.txt)
while read id; do [ -f "$DONE/$id/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt" ] && hc=$((hc+1)); done < /mnt/f/PD_brainage_data/hc_idlist.txt
echo "PD完成: $pd / 560     HC完成: $hc / 61"
echo
echo "=== 当前运行 ==="
echo "CIVET进程: $(pgrep -fc CIVET_Processing_Pipeline 2>/dev/null)   apptainer: $(pgrep -fc apptainer 2>/dev/null)"
echo "负载: $(cat /proc/loadavg)"
echo "ext4占用: $(df -h /mnt/civet 2>/dev/null | tail -1)"
echo "civet_out总目录: $(ls -d $DONE/*/ 2>/dev/null | wc -l)"
echo
echo "=== 最近完成的5例(时间戳) ==="
ls -dt $DONE/*/ 2>/dev/null | head -5 | while read d; do echo "  $(basename $d): $(stat -c %y "$d" 2>/dev/null | cut -d. -f1)"; done
