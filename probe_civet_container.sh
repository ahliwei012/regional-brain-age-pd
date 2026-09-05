#!/bin/bash
S=/mnt/g/apptainer/civet_2.1.1.sif
echo "== CIVET 安装目录 =="
apptainer exec "$S" bash -c 'ls -d /opt/CIVET* /opt/civet* 2>/dev/null'
echo
echo "== VERSION 文件 =="
apptainer exec "$S" bash -c 'cat /opt/CIVET*/VERSION 2>/dev/null || find /opt -iname "VERSION" 2>/dev/null | head'
echo
echo "== CIVET_Processing_Pipeline 位置 =="
apptainer exec "$S" bash -c 'which CIVET_Processing_Pipeline 2>/dev/null; ls /opt/CIVET*/CIVET_Processing_Pipeline 2>/dev/null'
echo
echo "== 容器内 PATH 是否已含 CIVET (默认env) =="
apptainer exec "$S" bash -c 'echo PATH=$PATH' | tr ':' '\n' | grep -i civet || echo "PATH 未含 CIVET (可能需 source 初始化)"
echo
echo "== 找初始化脚本 =="
apptainer exec "$S" bash -c 'ls /opt/CIVET*/*init* /opt/minc*/minc-toolkit-config.sh /opt/quarantine* 2>/dev/null; find /opt -iname "*init-sh*" -o -iname "minc-toolkit-config.sh" 2>/dev/null | head'
echo
echo "== CIVET 帮助(前40行, 看真实选项名) =="
apptainer exec "$S" bash -c 'CIVET_Processing_Pipeline -help 2>&1 || /opt/CIVET*/CIVET_Processing_Pipeline -help 2>&1' | head -40
