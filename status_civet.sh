#!/bin/bash
echo "=== 当前CIVET进程数 ==="
pgrep -fc CIVET_Processing_Pipeline 2>/dev/null || echo 0
echo
echo "=== 各例当前阶段 (从civetlog末尾) ==="
for f in /mnt/civet/out/*.civetlog; do
  [ -f "$f" ] || continue
  id=$(basename "$f" .civetlog)
  stage=$(grep -oE 'status of [a-z0-9_]+ in pipe' "$f" | tail -1 | sed 's/status of //; s/ in pipe//')
  done_n=$(grep -c 'to finished' "$f" 2>/dev/null)
  echo "  $id : 阶段=$stage  已完成阶段数=$done_n"
done
echo
echo "=== 运行时长/负载 ==="
uptime
echo "=== ext4 scratch 占用 ==="
df -h /mnt/civet 2>/dev/null | tail -1
