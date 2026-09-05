#!/bin/bash
DONE=/mnt/f/PD_brainage_data/civet_out
OUT=/mnt/f/PD_brainage_data
tf(){ [ -f "$DONE/$1/PPMI_${1}_native_rms_rsl_tlink_20mm_left.txt" ]; }

# 新HC完成
> "$OUT/hc_new_done.txt"
while read id; do [ -n "$id" ] && tf "$id" && echo "$id" >> "$OUT/hc_new_done.txt"; done < "$OUT/hc_new_idlist.txt"
echo "新HC完成: $(wc -l < "$OUT/hc_new_done.txt") / 221"

# 全部HC完成(原60 + 新, 去重)
cat "$OUT/hc_idlist.txt" "$OUT/hc_new_done.txt" | grep -v '^[[:space:]]*$' | sort -u > /tmp/hc_union.txt
> "$OUT/hc_all_done.txt"
while read id; do tf "$id" && echo "$id" >> "$OUT/hc_all_done.txt"; done < /tmp/hc_union.txt
echo "全部HC完成(去重): $(wc -l < "$OUT/hc_all_done.txt")"

# PD完成(已有pd_done)
echo "PD完成: $(wc -l < "$OUT/pd_done_idlist.txt")"

# 全体完成 = PD + 全HC
cat "$OUT/pd_done_idlist.txt" "$OUT/hc_all_done.txt" | grep -v '^[[:space:]]*$' | sort -u > "$OUT/all_final_idlist.txt"
echo "全体完成: $(wc -l < "$OUT/all_final_idlist.txt")"
# 需新建特征的HC(不在features_v2里的)
> "$OUT/hc_need_feat.txt"
while read id; do [ -f "$OUT/features_v2/${id}_features_20k.txt" ] || echo "$id" >> "$OUT/hc_need_feat.txt"; done < "$OUT/hc_all_done.txt"
echo "需新建特征的HC: $(wc -l < "$OUT/hc_need_feat.txt")"
