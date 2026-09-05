#!/bin/bash
# 在 G 盘建一个 ext4 镜像并 loop 挂载, 作为 CIVET 的高速 scratch
# (ext4 内核态批量写, 远快于 drvfs 小文件; 物理文件在 G 盘, 不占 C)
PW="$1"
IMG=/mnt/g/civet_ext4.img
MNT=/mnt/civet
SIZE=300G

if [ ! -f "$IMG" ]; then
  echo "[1] 创建 ${SIZE} 稀疏镜像 $IMG"
  truncate -s "$SIZE" "$IMG"
  echo "[2] 格式化 ext4"
  echo "$PW" | sudo -S -p '' mkfs.ext4 -q -F "$IMG"
fi

echo "[3] 挂载到 $MNT"
echo "$PW" | sudo -S -p '' mkdir -p "$MNT"
if ! mountpoint -q "$MNT"; then
  echo "$PW" | sudo -S -p '' mount -o loop "$IMG" "$MNT"
fi
echo "$PW" | sudo -S -p '' chown -R "$USER":"$USER" "$MNT"
mkdir -p "$MNT/sourcedir" "$MNT/out"

echo "[4] 结果:"
mountpoint "$MNT" && df -h "$MNT" | tail -1
echo "ext4-on-G 就绪: $MNT  (实际占用随写入增长, 稀疏)"
