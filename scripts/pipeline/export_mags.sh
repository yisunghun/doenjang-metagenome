#!/usr/bin/env bash
# 대표 MAG 63개를 저장소 mags/ 폴더로 gzip 압축해 복사 (WSL에서 저장소 루트를 인자로 실행)
# 사용: bash scripts/pipeline/export_mags.sh /path/to/doenjang-metagenome
set -euo pipefail
REPO=${1:?저장소 경로를 지정하세요}
REPS="$HOME/meta_hong/mag_catalog/reps"
mkdir -p "$REPO/mags"
n=0
for f in "$REPS"/*.fa; do
  b=$(basename "$f")
  gzip -c "$(readlink -f "$f")" > "$REPO/mags/${b}.gz"; n=$((n+1))
done
cp "$HOME/meta_hong/mag_catalog/reps_gtdb.tsv" "$REPO/mags/representatives_gtdbtk.tsv"
echo "exported $n MAGs → $REPO/mags ($(du -sh "$REPO/mags" | cut -f1))"
