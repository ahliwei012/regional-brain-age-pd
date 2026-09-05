#!/bin/bash
echo "civet_proc=$(pgrep -fc CIVET_Processing_Pipeline 2>/dev/null)"
echo "nii2mnc_proc=$(pgrep -fc nii2mnc 2>/dev/null)"
echo "apptainer_proc=$(pgrep -fc apptainer 2>/dev/null)"
echo "mnc_in_scratch=$(ls /mnt/civet/sourcedir/*.mnc 2>/dev/null | wc -l)"
echo "subjdir_in_scratch=$(ls -d /mnt/civet/out/*/ 2>/dev/null | wc -l)"
echo "mounted=$(mountpoint -q /mnt/civet && echo yes || echo NO)"
echo "load=$(cat /proc/loadavg)"
echo "--- 批次日志尾 ---"
tail -4 /mnt/f/PD_brainage_data/logs/civet_hc_batch.log
