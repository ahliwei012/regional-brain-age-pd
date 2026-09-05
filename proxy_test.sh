#!/bin/bash
# 测试 WSL 能否够到 Windows 上的 Clash 代理 (7892)
GW=$(ip route show default | grep -oP 'via \K[0-9.]+' | head -1)
echo "gateway(Windows host) IP: $GW"
for P in "127.0.0.1:7892" "${GW}:7892"; do
  echo "--- proxy http://$P ---"
  code=$(curl -x "http://$P" -m 8 -s -o /dev/null -w "%{http_code}" https://registry-1.docker.io/v2/ 2>/dev/null)
  echo "  result: HTTP ${code:-FAIL}  (401=可达成功, 000/空=不通)"
done
echo "--- 直连(无代理) ---"
code=$(curl -m 8 -s -o /dev/null -w "%{http_code}" https://registry-1.docker.io/v2/ 2>/dev/null)
echo "  result: HTTP ${code:-FAIL}"
