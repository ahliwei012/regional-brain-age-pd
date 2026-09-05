#!/bin/bash
# 临时把 WSL DNS 改为真实公共DNS, 绕过 Clash fake-ip
PW="$1"
printf 'nameserver 223.5.5.5\nnameserver 119.29.29.29\nnameserver 8.8.8.8\n' > /tmp/resolv.new
echo "$PW" | sudo -S -p '' cp /etc/resolv.conf /etc/resolv.conf.bak 2>/dev/null
echo "$PW" | sudo -S -p '' cp /tmp/resolv.new /etc/resolv.conf
echo "=== new resolv.conf ==="
cat /etc/resolv.conf
echo "=== resolve test (应为真实IP, 不再是198.18.x) ==="
echo -n "daocloud: "; getent hosts docker.m.daocloud.io || echo FAIL
echo -n "registry-1.docker.io: "; getent hosts registry-1.docker.io || echo FAIL
echo -n "dockerproxy: "; getent hosts dockerproxy.com || echo FAIL
