#!/bin/bash
LOG=/mnt/f/PD_brainage_data/logs/civet_newhc.log
DONE=/mnt/f/PD_brainage_data/civet_out
echo "ok=$(grep -c '\[ok\]' "$LOG" 2>/dev/null)  FAIL=$(grep -c '\[FAIL\]' "$LOG" 2>/dev/null)  err=$(grep -c '\[err\]' "$LOG" 2>/dev/null)"
n=0; while read id; do [ -n "$id" ] && [ -f "$DONE/$id/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt" ] && n=$((n+1)); done < /mnt/f/PD_brainage_data/hc_new_idlist.txt
echo "新HC完成: $n / 221"
echo "CIVET进程: $(pgrep -fc CIVET_Processing_Pipeline 2>/dev/null)  负载: $(cut -d' ' -f1-3 /proc/loadavg)"
echo "最近完成5例:"; ls -dt $DONE/*/ 2>/dev/null | head -3 | while read d; do echo "  $(basename $d): $(stat -c %y "$d" 2>/dev/null|cut -d. -f1)"; done
