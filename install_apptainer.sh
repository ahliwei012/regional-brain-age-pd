#!/bin/bash
# 自动安装 Apptainer v1.5.2 到 WSL Ubuntu 22.04
# 用法: bash install_apptainer.sh <sudo密码>
set -e
PW="$1"
WORK=/mnt/g/_apptainer_tmp
mkdir -p "$WORK"; cd "$WORK"

echo "[1/4] 下载 .deb 到 $WORK"
wget -q https://github.com/apptainer/apptainer/releases/download/v1.5.2/apptainer_1.5.2_amd64.deb -O apptainer_1.5.2_amd64.deb
ls -lh apptainer_1.5.2_amd64.deb

echo "[2/4] apt-get update"
echo "$PW" | sudo -S apt-get update -y

echo "[3/4] 安装 apptainer (含依赖)"
if ! ( echo "$PW" | sudo -S apt-get install -y "$WORK/apptainer_1.5.2_amd64.deb" ); then
  echo "apt 直装失败, 回退 dpkg + fix deps"
  echo "$PW" | sudo -S dpkg -i "$WORK/apptainer_1.5.2_amd64.deb" || true
  echo "$PW" | sudo -S apt-get -f install -y
fi

echo "[4/4] 验证"
apptainer --version
echo "DONE_OK"
