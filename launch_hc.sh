#!/bin/bash
# HC CIVET 批次启动器
# 用法: bash launch_hc.sh <并行数> <ID列表文件>
sed -i 's/\r$//' /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh
bash /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh \
  "0" "$1" \
  "/mnt/f/Analysis_HC_BL_ppmi10_13/T1Img" \
  "$2"
