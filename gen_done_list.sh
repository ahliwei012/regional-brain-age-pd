#!/bin/bash
DONE=/mnt/f/PD_brainage_data/civet_out
OUT=/mnt/f/PD_brainage_data
> "$OUT/pd_done_idlist.txt"
while read id; do
  [ -n "$id" ] && [ -f "$DONE/$id/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt" ] && echo "$id" >> "$OUT/pd_done_idlist.txt"
done < <(cat "$OUT/pd_p1_idlist.txt" "$OUT/pd_p2_idlist.txt")
echo "PD完成数: $(wc -l < "$OUT/pd_done_idlist.txt")"
comm -23 <(cat "$OUT/pd_p1_idlist.txt" "$OUT/pd_p2_idlist.txt" | sort -u) <(sort -u "$OUT/pd_done_idlist.txt") > "$OUT/pd_fails_idlist.txt"
echo "未完成(待重试): $(wc -l < "$OUT/pd_fails_idlist.txt")"
echo -n "失败ID: "; cat "$OUT/pd_fails_idlist.txt" | tr '\n' ' '; echo
