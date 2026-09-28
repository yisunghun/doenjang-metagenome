#!/usr/bin/env bash
# =============================================================================
# 06_mag_abundance_meta_hong.sh
#   32개 샘플 MAG → 종 수준(ANI 95%) 대표 유전체 카탈로그 → 샘플별 리드 매핑 존재비
#   목적: Kraken2 DB에 없는 미기재 종(예: Kroppenstedtia sp959600065)까지 포함한 정량
#
#   1) 수집 : 각 샘플 CheckM 결과에서 완성도 ≥50%, 오염 <10% bin (MIMAG medium 이상)
#   2) 대표 : galah로 ANI 95% 클러스터링 → 클러스터별 대표 MAG (품질 점수 최고)
#   3) 정량 : CoverM genome (minimap2-sr, identity ≥95%, aligned ≥75%)
#   4) 병합 : 샘플×MAG 상대존재비 표 + GTDB 분류 부착
#
# 사용법
#   bash 06_mag_abundance_meta_hong.sh --setup     # conda 환경 생성(최초 1회)
#   bash 06_mag_abundance_meta_hong.sh --catalog   # 1)+2) (수십 분)
#   nohup bash 06_mag_abundance_meta_hong.sh > mag_abundance.log 2>&1 &   # 3)+4)
# =============================================================================
set -uo pipefail

RESULTS_DIR="/mnt/f/meta_hong/results"
WORK="$HOME/meta_hong/mag_catalog"          # 카탈로그는 ext4(빠름)
ENV="mag_abundance"
THREADS=$(nproc)
MIN_COMP=50; MAX_CONT=10; ANI=95
OUT_SUB="13_mag_abundance_coverm"

SAMPLES=(Me_1_p Me_2_p Me_3_p Me_4_p Me_5_p Me_6_p
         0M_1_p 0M_2_p 0M_3_p 0M_4_p 0M_5_p 0M_6_p
         1M_1_p 1M_2_p 1M_3_p 1M_4_p 1M_5_p 1M_6_p
         3M_1_p 3M_2_p 3M_3_p 3M_4_p
         6M_1 6M_2 6M_3 6M_4 6M_5
         12M_1_p 12M_2_p 12M_3_p 12M_5_p 12M_6_p)

log(){ echo "[$(date '+%F %T')] $*"; }
die(){ log "ERROR: $*"; exit 1; }
act(){ source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "$ENV" || die "환경 $ENV 없음 (--setup 먼저)"; }

setup(){
  conda create -y -n "$ENV" -c conda-forge -c bioconda coverm galah minimap2 samtools || die "환경 생성 실패"
  log "완료: conda activate $ENV; coverm --version; galah --version"
}

catalog(){
  act
  mkdir -p "$WORK/all_bins"
  local info="$WORK/genomeInfo.csv"; echo "genome,completeness,contamination" > "$info"
  local n=0
  for s in "${SAMPLES[@]}"; do
    local qa="$RESULTS_DIR/$s/03_assembly/$s.checkm_qa.tsv"
    [[ -s "$qa" ]] || { log "WARNING: CheckM 결과 없음 $s"; continue; }
    # "(UIDxxx)" 뒤 4번째/5번째 필드가 Completeness/Contamination
    while read -r bin comp cont; do
      f="$RESULTS_DIR/$s/03_assembly/bins/$bin.fa"
      [[ -s "$f" ]] || continue
      cp -n "$f" "$WORK/all_bins/$bin.fa"
      echo "$bin.fa,$comp,$cont" >> "$info"; n=$((n+1))
    done < <(awk -v mc=$MIN_COMP -v mx=$MAX_CONT '
        $1 ~ /\.bin\./ { for(i=2;i<=NF;i++) if($i ~ /\(UID/){c=$(i+4); x=$(i+5); break}
                         if(c>=mc && x<mx) print $1, c, x }' "$qa")
  done
  log "수집된 MAG (완성도≥$MIN_COMP, 오염<$MAX_CONT): $n"

  # all_bins에 이전 실행의 잔여 파일이 섞이지 않도록 이번 목록에 있는 것만 남김
  local keep="$WORK/keep.txt"; tail -n +2 "$info" | cut -d, -f1 | sort -u > "$keep"
  for f in "$WORK"/all_bins/*.fa; do grep -qxF "$(basename "$f")" "$keep" || rm -f "$f"; done
  log "클러스터링 대상: $(ls "$WORK/all_bins" | wc -l)개 (품질표 $(wc -l < "$keep")행)"

  # galah 버전에 따라 genome-info의 이름 형식(확장자 없음/파일명/전체경로)이 달라 차례로 시도
  local ok=0
  for fmt in stem name path; do
    local gi="$WORK/genomeInfo.$fmt.csv"; echo "genome,completeness,contamination" > "$gi"
    tail -n +2 "$info" | while IFS=, read -r g c x; do
      case $fmt in
        stem) echo "${g%.fa},$c,$x" ;;
        name) echo "$g,$c,$x" ;;
        path) echo "$WORK/all_bins/$g,$c,$x" ;;
      esac
    done >> "$gi"
    log "galah 클러스터링 (ANI $ANI%, genome-info 이름 형식: $fmt)"
    rm -rf "$WORK/reps"
    if galah cluster --genome-fasta-directory "$WORK/all_bins" --genome-fasta-extension fa \
                     --genome-info "$gi" --ani "$ANI" --precluster-ani 90 \
                     --min-completeness "$MIN_COMP" --max-contamination "$MAX_CONT" \
                     --output-cluster-definition "$WORK/clusters.tsv" \
                     --output-representative-fasta-directory "$WORK/reps" \
                     --threads "$THREADS" 2> "$WORK/galah.$fmt.log"; then
      ok=1; log "galah 성공 (형식: $fmt)"; break
    else
      log "  형식 $fmt 실패: $(grep -m1 -o 'Failed to find CheckM statistics for [^ ]*' "$WORK/galah.$fmt.log")"
    fi
  done
  (( ok )) || die "galah 3가지 형식 모두 실패 → $WORK/galah.*.log 확인"
  log "대표 MAG(종 수준 클러스터): $(ls "$WORK/reps" | wc -l)개 → $WORK/reps"

  # 대표 MAG의 GTDB 분류 모으기
  : > "$WORK/reps_gtdb.tsv"
  for f in "$WORK"/reps/*.fa; do
    b=$(basename "$f" .fa); s=${b%%.bin.*}
    cls=$(cat "$RESULTS_DIR/$s"/09_mag_taxonomy_gtdbtk/gtdbtk_output_filtered/gtdbtk.*.summary.tsv 2>/dev/null \
          | awk -F'\t' -v b="$b" '$1==b{print $2; exit}')
    printf "%s\t%s\n" "$b" "${cls:-unclassified}" >> "$WORK/reps_gtdb.tsv"
  done
  log "GTDB 분류표: $WORK/reps_gtdb.tsv"
}

quantify(){
  act
  [[ -d "$WORK/reps" ]] || die "카탈로그 없음 (--catalog 먼저)"
  local fail=()
  for s in "${SAMPLES[@]}"; do
    local out="$RESULTS_DIR/$s/$OUT_SUB"; mkdir -p "$out"
    [[ -f "$out/.done" ]] && { log "SKIP $s"; continue; }
    log "START $s"
    coverm genome \
      -1 "$RESULTS_DIR/$s/02_host_removed/${s}_1.hostremoved.fastq.gz" \
      -2 "$RESULTS_DIR/$s/02_host_removed/${s}_2.hostremoved.fastq.gz" \
      --genome-fasta-directory "$WORK/reps" -x fa \
      --mapper minimap2-sr --min-read-percent-identity 95 --min-read-aligned-percent 75 \
      -m relative_abundance covered_fraction mean \
      --min-covered-fraction 10 -t "$THREADS" \
      -o "$out/$s.coverm.tsv" 2> "$out/$s.coverm.log" \
      && touch "$out/.done" && log "DONE  $s  (unmapped: $(awk -F'\t' '$1=="unmapped"{print $2"%"}' "$out/$s.coverm.tsv"))" \
      || { log "ERROR $s (로그: $out/$s.coverm.log)"; fail+=("$s"); }
  done

  # 병합: 행=MAG, 열=샘플, 값=relative_abundance(%) ; 'unmapped' 행 포함
  mkdir -p "$RESULTS_DIR/v3_merged"
  local merged="$RESULTS_DIR/v3_merged/mag_relative_abundance.tsv"
  awk -F'\t' '
    FNR==1{ split(FILENAME,p,"/"); s=p[length(p)]; sub(/\.coverm\.tsv$/,"",s); S[++ns]=s; next }
    { g[$1]=1; v[$1,s]=$2 }
    END{ printf "genome"; for(i=1;i<=ns;i++) printf "\t%s",S[i]; print "";
         for(k in g){ printf "%s",k; for(i=1;i<=ns;i++) printf "\t%s",((k,S[i]) in v? v[k,S[i]]:0); print "" } }' \
    "$RESULTS_DIR"/*/"$OUT_SUB"/*.coverm.tsv > "$merged"
  cp "$WORK/reps_gtdb.tsv" "$WORK/clusters.tsv" "$RESULTS_DIR/v3_merged/" 2>/dev/null
  log "병합표: $merged"
  (( ${#fail[@]} )) && { log "실패: ${fail[*]} (재실행 시 실패분만)"; exit 1; }
  log "===== 32개 완료 ====="
}

case "${1:-}" in
  --setup)   setup ;;
  --catalog) catalog ;;
  *)         quantify ;;
esac
