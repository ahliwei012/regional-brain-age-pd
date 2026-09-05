#!/bin/bash
RAWDIR="/mnt/f/Analysis_HC_BL_ppmi10_13/T1Img"
OUTDIR=/mnt/f/PD_brainage_data
echo "RAWDIR可读: $([ -d "$RAWDIR" ] && echo yes || echo NO)"
ls -1 "$RAWDIR" 2>/dev/null | sort > "$OUTDIR/hc_idlist.txt"
echo "HC总数: $(wc -l < "$OUTDIR/hc_idlist.txt")"
echo "前5 ID:"; head -5 "$OUTDIR/hc_idlist.txt"
echo "抽查nii:"
for id in $(head -3 "$OUTDIR/hc_idlist.txt"); do
  ls "$RAWDIR/$id/"*.nii 2>/dev/null | head -1 || echo "  $id 无nii"
done
