#!/bin/sh
# 打开 GUI，FoV 用 fov.json
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
exec "$ROOT/bin/sccm" show \
  --obs1 "$HERE/data/aia304_20250416_cutout.fits" \
  --obs2 "$HERE/data/euvi304_20250416_cutout.fts" \
  --out "$HERE/output" --fov-json "$HERE/fov.json"
