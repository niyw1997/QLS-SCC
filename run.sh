#!/bin/sh
# 最直接的两幅 FITS 启动方式
#   ./run.sh /path/to/aia.fits /path/to/euvi.fits [输出目录]
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
if [ $# -lt 2 ]; then
    echo "用法: $0 /path/to/aia.fits /path/to/euvi.fits [输出目录]" >&2
    exit 2
fi
A=$1
B=$2
OUT=${3:-./out}
exec "$HERE/bin/sccm" show --obs1 "$A" --obs2 "$B" --out "$OUT"
