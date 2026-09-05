#!/bin/bash
NEW=/mnt/g/PPMI_HC_XING/BL/3DT1_only
OUT=/mnt/f/PD_brainage_data
echo "=== 3DT1_only 子目录数 ==="
ls -1 "$NEW" 2>/dev/null | wc -l
echo "=== 抽查文件 ==="
for id in $(ls -1 "$NEW" | head -3); do ls -la "$NEW/$id/" 2>/dev/null | grep -iE 'nii'; done
echo
echo "=== 与现有61 HC的重叠/新增 ==="
ls -1 "$NEW" 2>/dev/null | sort > /tmp/newhc_all.txt
comm -12 /tmp/newhc_all.txt <(sort "$OUT/hc_idlist.txt") | wc -l | xargs echo "与现有HC重叠:"
comm -23 /tmp/newhc_all.txt <(sort "$OUT/hc_idlist.txt") > "$OUT/hc_new_idlist.txt"
echo "新增HC(不在现有61里): $(wc -l < "$OUT/hc_new_idlist.txt")"
echo "=== 也检查是否和PD重复(应0) ==="
comm -12 /tmp/newhc_all.txt <(cat "$OUT/pd_p1_idlist.txt" "$OUT/pd_p2_idlist.txt" | sort -u) | wc -l | xargs echo "与PD重叠:"
echo "=== 新增HC前5 + nii.gz确认 ==="
for id in $(head -5 "$OUT/hc_new_idlist.txt"); do ls "$NEW/$id/"*.nii.gz 2>/dev/null | head -1; done
echo "=== 另: BL/T1Img 是什么 ==="
ls -1 /mnt/g/PPMI_HC_XING/BL/T1Img 2>/dev/null | head -3
ls /mnt/g/PPMI_HC_XING/BL/T1Img/$(ls -1 /mnt/g/PPMI_HC_XING/BL/T1Img 2>/dev/null | head -1)/ 2>/dev/null | head -3
