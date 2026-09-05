#!/bin/bash
# 全量PD CIVET: P1(362)->P2(198), 高并行, 断点续跑
# 用法: bash launch_pd_all.sh <并行数>
NPAR="$1"
sed -i 's/\r$//' /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh
echo "######## P1 (362) start $(date) ########"
bash /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh "0" "$NPAR" \
  "/mnt/g/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P1/T1Img" \
  "/mnt/f/PD_brainage_data/pd_p1_idlist.txt"
echo "######## P2 (198) start $(date) ########"
bash /mnt/f/PD_brainage_data/scripts/run_civet_batch.sh "0" "$NPAR" \
  "/mnt/g/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P2/T1Img" \
  "/mnt/f/PD_brainage_data/pd_p2_idlist.txt"
echo "######## 全量PD完成 $(date) ########"
