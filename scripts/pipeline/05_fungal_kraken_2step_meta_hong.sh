#!/usr/bin/env bash
# =============================================================================
# 05_fungal_kraken_2step_meta_hong.sh
#   Amazon S3(PlusPF) 차단 환경용 대안: NCBI RefSeq 진균 전용 Kraken2 DB를 직접 구축하고,
#   기존 Standard-16 결과에서 '미분류(U)'로 남은 리드만 2단계로 진균 분류
#
#   1단계(완료됨) : Standard-16 → 세균/고균/바이러스  (05_taxonomy_reads_kraken/)
#   2단계(이 스크립트): 1단계 미분류 리드 → RefSeq Fungi DB (05c_taxonomy_fungi_2step/)
#
# 사용법 (각 단계 독립 실행 가능, 완료된 단계는 자동 건너뜀)
#   bash 05_fungal_kraken_2step_meta_hong.sh --build-db   # 다운로드+DB빌드+Bracken빌드 (수 시간)
#   bash 05_fungal_kraken_2step_meta_hong.sh --check
#   nohup bash 05_fungal_kraken_2step_meta_hong.sh > fungi_2step.log 2>&1 &
# =============================================================================
set -uo pipefail

# ---------------------------- [설정] ----------------------------------------
RESULTS_DIR="/mnt/f/meta_hong/results"
CONDA_ENV="metagenomics_shotgun_kraken"
THREADS=$(nproc)
READ_LEN=100
CONFIDENCE=0.1          # 진균 전용 DB는 세균 리드의 우연한 k-mer 일치를 줄이기 위해 0.1 적용
MIN_HIT_GROUPS=3        # kraken2 --minimum-hit-groups
BRACKEN_THRESHOLD=10
DL_PARALLEL=8           # 유전체 동시 다운로드 수

# RefSeq 진균 유전체 선택 범위: "all"(최신판 전체, Mucor/Rhizopus 등 scaffold 수준 포함, 권장)
#                               "reference"(reference/representative genome만, 빠름)
ASM_SCOPE="all"

DB_ROOT="$HOME/db/kraken2"
DB_DIR="${DB_ROOT}/k2_refseq_fungi"
[[ -n "${FUNGI_DB:-}" ]] && DB_DIR="$FUNGI_DB"        # 이미 빌드했다면 FUNGI_DB=경로 로 지정 가능
OUT_SUB="05c_taxonomy_fungi_2step"
NCBI="https://ftp.ncbi.nlm.nih.gov"
# ---------------------------------------------------------------------------

log(){ echo "[$(date '+%F %T')] $*"; }
die(){ log "ERROR: $*"; exit 1; }

activate_env(){
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "$CONDA_ENV" \
    || die "conda 환경 '$CONDA_ENV' 활성화 실패"
}

# 방화벽 차단 페이지(HTML, 200 OK)를 진짜 파일로 오인하지 않도록 내용 검증
fetch(){  # fetch URL OUT
  wget -q -c --tries=5 --timeout=60 "$1" -O "$2" || return 1
  if head -c 300 "$2" | grep -qi '<html\|<iframe\|warning.html'; then
    rm -f "$2"; die "방화벽 차단 페이지 수신: $1"
  fi
}

build_db(){
  activate_env
  for t in kraken2-build bracken-build wget; do command -v $t >/dev/null || die "$t 없음"; done
  mkdir -p "$DB_DIR"/{taxonomy,library,_download/genomes}
  local dl="$DB_DIR/_download"

  # 1) 택소노미 (최신 taxdump; 2024년판 Standard-16의 names/nodes에는 최근 진균 taxid가 없을 수 있음)
  if [[ ! -s "$DB_DIR/taxonomy/nodes.dmp" ]]; then
    log "taxdump 다운로드"
    fetch "$NCBI/pub/taxonomy/taxdump.tar.gz" "$dl/taxdump.tar.gz" || die "taxdump 실패"
    tar -xzf "$dl/taxdump.tar.gz" -C "$DB_DIR/taxonomy" nodes.dmp names.dmp || die "taxdump 해제 실패"
  fi

  # 2) RefSeq 진균 assembly 목록
  fetch "$NCBI/genomes/refseq/fungi/assembly_summary.txt" "$dl/assembly_summary.txt" || die "assembly_summary 실패"
  grep -q '^#.*assembly_accession' "$dl/assembly_summary.txt" || die "assembly_summary 형식 이상(차단 여부 확인)"

  # 열: 1 accession, 5 refseq_category, 6 taxid, 11 version_status, 12 assembly_level, 20 ftp_path
  awk -F'\t' -v scope="$ASM_SCOPE" '
    /^#/ {next}
    $11=="latest" && $20!="na" {
      if (scope=="reference" && $5!="reference genome" && $5!="representative genome") next
      p=$20; sub(/\r$/,"",p); sub(/\/+$/,"",p)          # 끝의 "/" 제거 (NCBI 경로가 /로 끝나는 경우)
      n=split(p,a,"/"); print $6"\t"p"/"a[n]"_genomic.fna.gz"
    }' "$dl/assembly_summary.txt" > "$dl/selected.tsv"
  log "선택된 진균 assembly: $(wc -l < "$dl/selected.tsv")개 (scope=$ASM_SCOPE)"
  log "  장류 핵심 속 포함 여부:"
  for g in Aspergillus Mucor Rhizopus Saccharomyces Zygosaccharomyces Debaryomyces Candida Penicillium Wickerhamomyces Pichia; do
    printf "    %-18s %s\n" "$g" "$(awk -F'\t' -v g="$g" '!/^#/ && $11=="latest" && $8 ~ "^"g" " ' "$dl/assembly_summary.txt" | wc -l)"
  done

  # 3) 유전체 다운로드 (병렬, 이어받기)
  log "유전체 다운로드 시작 (병렬 $DL_PARALLEL)"
  # 첫 1개를 자세히 시험 → 실패 원인을 로그에 남기고 중단
  local first; first=$(head -1 "$dl/selected.tsv" | cut -f2)
  local ff="$dl/genomes/$(basename "$first")"
  if [[ ! -s "$ff" ]]; then
    log "시험 다운로드: $first"
    if ! wget --tries=3 --timeout=60 -O "$ff" "$first" 2> "$dl/test_wget.log"; then
      rm -f "$ff"; tail -15 "$dl/test_wget.log"; die "시험 다운로드 실패 (위 wget 메시지 확인)"
    fi
    if ! gzip -t "$ff" 2>/dev/null; then
      head -c 400 "$ff"; echo; rm -f "$ff"; die "받은 파일이 gzip이 아님 (방화벽 차단 페이지 가능성)"
    fi
    log "시험 다운로드 성공 ($(du -h "$ff" | cut -f1))"
  fi

  export dl
  cut -f2 "$dl/selected.tsv" | xargs -P "$DL_PARALLEL" -n1 bash -c '
      u="$0"; f="$dl/genomes/${u##*/}"
      [ -s "$f" ] && exit 0
      wget -q --tries=5 --timeout=120 -O "$f.part" "$u" && gzip -t "$f.part" 2>/dev/null && mv "$f.part" "$f" \
        || { rm -f "$f.part"; echo "FAIL $u" >> "$dl/failed_downloads.txt"; }'
  local got total; got=$(ls "$dl/genomes" | grep -c '\.fna\.gz$'); total=$(wc -l < "$dl/selected.tsv")
  log "다운로드 완료: $got / $total  (실패 목록: $dl/failed_downloads.txt)"
  (( got * 100 >= total * 90 )) || die "다운로드 90% 미만 → 라이브러리 구성 중단. 다시 실행하면 이어받기"

  # 4) taxid를 헤더에 넣어 라이브러리에 추가 (accession2taxid 대용량 파일 불필요)
  #    kraken2-build 내장 마스킹은 21GB 단일 파일을 코어 1개로 처리해 매우 느리고 중간에 끊기기 쉬움
  #    → 유전체별로 병렬 마스킹(dustmasker 또는 k2mask) 후 --no-masking 으로 추가
  if [[ ! -f "$DB_DIR/.library_done" ]]; then
    rm -rf "$DB_DIR/library"                              # 이전 미완료 라이브러리 정리
    local masker=""
    command -v dustmasker >/dev/null && masker="dustmasker"
    [[ -z "$masker" ]] && command -v k2mask >/dev/null && masker="k2mask"
    [[ -n "$masker" ]] || die "dustmasker/k2mask 없음 → conda install -n $CONDA_ENV -c bioconda blast"
    log "라이브러리 구성: 유전체별 taxid 태깅 + 저복잡도 마스킹($masker) 병렬 $THREADS"
    mkdir -p "$dl/masked"; export dl masker
    xargs -P "$THREADS" -n2 bash -c '
        taxid="$0"; u="$1"; b="${u##*/}"; g="$dl/genomes/$b"; o="$dl/masked/${b%.gz}"
        [ -s "$g" ] || exit 0
        [ -s "$o" ] && exit 0
        zcat "$g" | sed "s/^>/>kraken:taxid|${taxid}|/" \
          | $masker -in /dev/stdin -outfmt fasta 2>/dev/null \
          | sed -e "/^>/!s/[a-z]/x/g" > "$o.part" && [ -s "$o.part" ] && mv "$o.part" "$o" \
          || { rm -f "$o.part"; echo "MASK_FAIL $b" >&2; }' < "$dl/selected.tsv"
    local nm; nm=$(ls "$dl/masked" | grep -c '\.fna$')
    log "마스킹 완료: $nm / $(ls "$dl/genomes" | grep -c '\.fna\.gz$')"
    (( nm * 100 >= $(wc -l < "$dl/selected.tsv") * 90 )) || die "마스킹 90% 미만 (다시 실행하면 남은 것만 처리)"

    local tmp="$dl/tagged.masked.fna"
    cat "$dl"/masked/*.fna > "$tmp"
    kraken2-build --add-to-library "$tmp" --db "$DB_DIR" --no-masking --threads "$THREADS" \
      || die "add-to-library 실패"
    rm -f "$tmp" "$dl/tagged.fna"; touch "$DB_DIR/.library_done"
    log "라이브러리 완료 (확인 후 정리 가능: rm -rf $dl/masked)"
  fi

  # 5) DB 빌드
  if [[ ! -s "$DB_DIR/hash.k2d" ]]; then
    log "kraken2-build --build (수 시간 소요 가능)"
    kraken2-build --build --db "$DB_DIR" --threads "$THREADS" || die "build 실패"
  fi

  # 6) Bracken k-mer 분포
  if [[ ! -s "$DB_DIR/database${READ_LEN}mers.kmer_distrib" ]]; then
    log "bracken-build -l $READ_LEN"
    bracken-build -d "$DB_DIR" -t "$THREADS" -k 35 -l "$READ_LEN" || die "bracken-build 실패"
  fi

  kraken2-inspect --db "$DB_DIR" > "$DB_DIR/inspect.txt" 2>/dev/null
  du -sh "$DB_DIR"
  log "DB 준비 완료: $DB_DIR"
  log "정리(선택): rm -rf $DB_DIR/_download/genomes $DB_DIR/library  # 수십 GB 확보"
}

preflight(){
  local ok=1
  for t in kraken2 bracken seqkit; do
    command -v $t >/dev/null || { log "ERROR: $t 없음"; [[ $t == seqkit ]] && log "  설치: conda install -n $CONDA_ENV -c bioconda seqkit"; ok=0; }
  done
  for f in hash.k2d taxo.k2d opts.k2d "database${READ_LEN}mers.kmer_distrib"; do
    [[ -s "$DB_DIR/$f" ]] || { log "ERROR: DB 파일 없음 $DB_DIR/$f (--build-db 먼저)"; ok=0; }
  done
  local n=0
  for d in "$RESULTS_DIR"/*/; do
    s=$(basename "$d")
    [[ -f "$d/05_taxonomy_reads_kraken/${s}.kraken2.out" ]] || continue
    [[ -f "$d/02_host_removed/${s}_1.hostremoved.fastq.gz" ]] || { log "ERROR: 입력 없음 $s"; ok=0; }
    n=$((n+1))
  done
  log "1단계 결과(kraken2.out) 있는 샘플: $n"
  (( n == 32 )) || log "WARNING: 32개가 아님"
  (( ok == 1 ))
}

run_sample(){
  local s=$1 d="$RESULTS_DIR/$1" out="$RESULTS_DIR/$1/$OUT_SUB"
  mkdir -p "$out"
  [[ -f "$out/.done" ]] && { log "SKIP $s"; return 0; }
  log "START $s"

  # (a) 1단계 미분류 리드 ID
  awk -F'\t' '$1=="U"{print $2}' "$d/05_taxonomy_reads_kraken/${s}.kraken2.out" > "$out/${s}.U.ids"
  local nu; nu=$(wc -l < "$out/${s}.U.ids")

  # (b) 미분류 리드만 추출 (임시, 끝나면 삭제)
  for r in 1 2; do
    seqkit grep -j "$THREADS" -f "$out/${s}.U.ids" "$d/02_host_removed/${s}_${r}.hostremoved.fastq.gz" \
      -o "$out/${s}_${r}.U.fastq.gz" || { log "ERROR seqkit $s"; return 1; }
  done

  # (c) 진균 DB 분류
  kraken2 --db "$DB_DIR" --threads "$THREADS" --paired --gzip-compressed \
          --confidence "$CONFIDENCE" --minimum-hit-groups "$MIN_HIT_GROUPS" \
          --report "$out/${s}.fungi.kraken2.report" --output /dev/null \
          "$out/${s}_1.U.fastq.gz" "$out/${s}_2.U.fastq.gz" 2> "$out/${s}.fungi.kraken2.log" \
    || { log "ERROR kraken2 $s"; return 1; }

  for lvl in S G; do
    bracken -d "$DB_DIR" -i "$out/${s}.fungi.kraken2.report" -o "$out/${s}.fungi.bracken_${lvl}.txt" \
            -w "$out/${s}.fungi.bracken_${lvl}.report" -r "$READ_LEN" -l "$lvl" -t "$BRACKEN_THRESHOLD" \
            >> "$out/${s}.fungi.bracken.log" 2>&1 || log "WARNING bracken $lvl $s (진균 리드가 매우 적으면 발생 가능)"
  done

  # (d) 요약: 전체 리드 대비 진균 비율 (분모 = 1단계 전체 리드)
  local tot; tot=$(awk -F'\t' '$4=="U"||$4=="R"{t+=$2}END{print t}' "$d/05_taxonomy_reads_kraken/${s}.kraken2.report")
  local nf;  nf=$(awk -F'\t' '$5==4751{print $2}' "$out/${s}.fungi.kraken2.report"); nf=${nf:-0}
  awk -v s="$s" -v t="$tot" -v u="$nu" -v f="$nf" 'BEGIN{
      printf "%s\ttotal=%d\tstep1_unclassified=%d (%.2f%%)\tfungi=%d (%.3f%% of total, %.2f%% of step1-U)\n", s,t,u,100*u/t,f,100*f/t,(u?100*f/u:0)}' \
    | tee "$out/${s}.fungi.summary.txt"
  echo "Top 5 진균 속:"; awk -F'\t' '$4=="G"' "$out/${s}.fungi.kraken2.report" | sort -t$'\t' -k2,2nr | head -5 | awk -F'\t' '{gsub(/^ +/,"",$6); printf "    %-22s %d\n",$6,$2}'

  rm -f "$out/${s}_1.U.fastq.gz" "$out/${s}_2.U.fastq.gz" "$out/${s}.U.ids"
  touch "$out/.done"; log "DONE  $s"
}

# ------------------------------- main ---------------------------------------
case "${1:-}" in
  --build-db) build_db; exit $? ;;
  --check)    activate_env; preflight; exit $? ;;
esac

activate_env
log "===== 진균 2단계 분류 시작 (DB: $DB_DIR) ====="
preflight || die "점검 실패"

ORDER=(Me_1_p Me_2_p Me_3_p Me_4_p Me_5_p Me_6_p
       0M_1_p 0M_2_p 0M_3_p 0M_4_p 0M_5_p 0M_6_p
       1M_1_p 1M_2_p 1M_3_p 1M_4_p 1M_5_p 1M_6_p
       3M_1_p 3M_2_p 3M_3_p 3M_4_p
       6M_1 6M_2 6M_3 6M_4 6M_5
       12M_1_p 12M_2_p 12M_3_p 12M_5_p 12M_6_p)
fail=()
for s in "${ORDER[@]}"; do run_sample "$s" || fail+=("$s"); done

mkdir -p "$RESULTS_DIR/v3_merged"
cat "$RESULTS_DIR"/*/"$OUT_SUB"/*.fungi.summary.txt > "$RESULTS_DIR/v3_merged/fungi_2step_summary.tsv" 2>/dev/null
log "요약표: $RESULTS_DIR/v3_merged/fungi_2step_summary.tsv"
(( ${#fail[@]} )) && { log "실패: ${fail[*]} (재실행 시 실패분만 재시도)"; exit 1; }
log "===== 32개 전부 완료 ====="
