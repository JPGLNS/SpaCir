#!/bin/bash
#
# decompress_fq.sh - 将指定目录下的 .fq.gz 文件解压为 fastq 文件
#
# 用法: bash decompress_fq.sh <工作目录>
# 示例: bash decompress_fq.sh /path/to/TRB
#
# 转换规则:
#   15-TRB_1.fq.gz  →  fastq/15_TRB_R1.fastq
#   23-IGH_2.fq.gz  →  fastq/23_IGH_R2.fastq
#   - → _,  _1 → _R1,  _2 → _R2
#

set -euo pipefail

# --- 参数检查 ---
if [ $# -ne 1 ]; then
    echo "用法: bash $(basename "$0") <工作目录>"
    exit 1
fi

WORKDIR="$1"

if [ ! -d "$WORKDIR" ]; then
    echo "错误: 目录不存在: $WORKDIR"
    exit 1
fi

cd "$WORKDIR"

# --- 查找 .fq.gz 文件 ---
shopt -s nullglob
files=( *.fq.gz )
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    echo "提示: 在 $WORKDIR 下未找到 .fq.gz 文件"
    exit 0
fi

echo "工作目录: $WORKDIR"
echo "找到 ${#files[@]} 个 .fq.gz 文件"
echo "---"

# --- 创建 fastq 输出目录 ---
mkdir -p fastq

# --- 逐个处理 ---
skipped=0
processed=0

for fq_gz in "${files[@]}"; do
    # 去掉 .fq.gz 扩展名
    base="${fq_gz%.fq.gz}"

    # 转换规则: - → _,  _1 → _R1,  _2 → _R2
    base="${base//-/_}"
    base="${base//_1/_R1}"
    base="${base//_2/_R2}"

    target="fastq/${base}.fastq"

    if [ -f "$target" ]; then
        echo "跳过: $target (已存在)"
        skipped=$((skipped + 1))
    else
        echo "解压: $fq_gz → $target"
        gzip -dc "$fq_gz" > "$target"
        processed=$((processed + 1))
    fi
done

echo "---"
echo "完成: 解压 $processed 个, 跳过 $skipped 个"
