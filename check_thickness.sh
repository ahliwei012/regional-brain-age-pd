#!/bin/bash
echo "=== CIVET厚度均值/范围 (模型scaler厚度均值=3.12mm 作参照) ==="
for id in 100006 100007 100018 100267 100268; do
  f="/mnt/f/PD_brainage_data/civet_out/$id/PPMI_${id}_native_rms_rsl_tlink_20mm_left.txt"
  [ -f "$f" ] || { echo "  $id 无文件"; continue; }
  awk -v id="$id" '{s+=$1; ss+=$1*$1; if($1>mx||NR==1)mx=$1; if($1<mn||NR==1)mn=$1; n++}
    END{m=s/n; sd=sqrt(ss/n-m*m); printf "  %s 左: mean=%.3f sd=%.3f min=%.3f max=%.3f n=%d\n", id, m, sd, mn, mx, n}' "$f"
done
