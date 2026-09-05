#!/bin/bash
RAWDIR="/mnt/g/小邢整理的原始PPMI影像PD数据/PD_BL_T1_BOLD_P2/T1Img"
OUTDIR=/mnt/f/PD_brainage_data
echo "P2 RAWDIR可读: $([ -d "$RAWDIR" ] && echo yes || echo NO)"
ls -1 "$RAWDIR" 2>/dev/null | sort > "$OUTDIR/pd_p2_idlist.txt"
echo "P2 PD总数: $(wc -l < "$OUTDIR/pd_p2_idlist.txt")"
echo "前3 ID + nii:"
for id in $(head -3 "$OUTDIR/pd_p2_idlist.txt"); do
  ls "$RAWDIR/$id/"*.nii 2>/dev/null | head -1 || echo "  $id 无nii"
done
echo "P1+P2 合计: $(( $(wc -l < "$OUTDIR/pd_p1_idlist.txt") + $(wc -l < "$OUTDIR/pd_p2_idlist.txt") ))"
# 检查P1/P2有无重复ID
dup=$(cat "$OUTDIR/pd_p1_idlist.txt" "$OUTDIR/pd_p2_idlist.txt" | sort | uniq -d | wc -l)
echo "P1/P2重复ID数: $dup"
