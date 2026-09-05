#!/bin/bash
# 启动器: 硬编码中文RAWDIR, 避免命令行传中文乱码
# 用法: bash launch_civet.sh <并行数> <ID列表文件>
sed -i 's/\r$//' /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh
bash /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh \
  "0" "$1" \
  "/mnt/g/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P1/T1Img" \
  "$2"
