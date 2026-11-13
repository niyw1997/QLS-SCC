#!/bin/sh
# 建 .venv 并安装依赖。macOS 若缺 Tk，先 brew install python-tk@3.12
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
"$PY" -m venv "$HERE/.venv"
"$HERE/.venv/bin/python" -m pip install --upgrade pip
"$HERE/.venv/bin/python" -m pip install -r "$HERE/requirements.txt"
echo
echo "安装完成，试一下："
echo "  $HERE/bin/sccm info --obs1 $HERE/examples/20250416/data/aia193_20250416_cutout.fits --obs2 $HERE/examples/20250416/data/euvi195_20250416_cutout.fts"
echo "  $HERE/examples/20250416/run.sh"
