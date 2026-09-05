#!/bin/bash
RAWDIR="/mnt/g/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P1/T1Img"
OUTDIR=/mnt/f/PD_brainage_data
ls -1 "$RAWDIR" 2>/dev/null | sort > "$OUTDIR/pd_p1_idlist.txt"
head -6 "$OUTDIR/pd_p1_idlist.txt" > "$OUTDIR/pd_test6_idlist.txt"
echo "RAWDIR可读: $([ -d "$RAWDIR" ] && echo yes || echo NO)"
echo "P1 PD总数: $(wc -l < "$OUTDIR/pd_p1_idlist.txt")"
echo "前6 ID:"; cat "$OUTDIR/pd_test6_idlist.txt"
echo "抽查nii路径:"
for id in $(head -3 "$OUTDIR/pd_test6_idlist.txt"); do
  ls "$RAWDIR/$id/"*.nii 2>/dev/null | head -1 || echo "  $id 无nii"
done
