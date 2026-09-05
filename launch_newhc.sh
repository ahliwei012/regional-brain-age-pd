#!/bin/bash
# 新增HC CIVET (221例, .nii.gz输入)
# 用法: bash launch_newhc.sh <并行数>
sed -i 's/\r$//' /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh
bash /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh \
  "0" "$1" \
  "/mnt/g/PPMI_HC_XING/BL/3DT1_only" \
  "/mnt/f/PD_brainage_data/hc_new_idlist.txt"
