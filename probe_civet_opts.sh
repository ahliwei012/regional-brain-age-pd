#!/bin/bash
S=/mnt/g/apptainer/civet_2.1.1.sif
echo "== nii2mnc / mnc2nii 是否可用 =="
apptainer exec "$S" bash -c 'which nii2mnc; which mnc2nii'
echo
echo "== 关键选项 (从 -help grep) =="
apptainer exec "$S" CIVET_Processing_Pipeline -help 2>&1 | \
  grep -iE 'thickness|N3-distance|resample-surface|surface-atlas|prefix|sourcedir|targetdir|-run|-spawn|input|stereotaxic|-q ' | head -50
