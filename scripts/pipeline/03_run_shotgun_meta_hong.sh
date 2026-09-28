#!/usr/bin/env bash
# ==============================================================================
#  Shotgun Metagenome Analysis Pipeline — meta_hong 전용판
#  Base    : run_shotgun_mnt.sh v3.0 (Production, sunghunyi, 2025-04) — 로직 그대로 재사용
#  Updated : 2026-09 (Claude) — meta_hong 작업공간용 경로 재구성 + 안전장치 추가
#
#  변경 사항 (원본 run_shotgun_mnt.sh 대비):
#   1. BASE_DIR(작업/스크립트 기준)는 ~/meta_hong(홈, 로컬 디스크) 고정,
#      RESULTS_DIR(실제 대용량 산출물)는 용량 문제로 F 드라이브의 신규 경로
#      (/mnt/f/meta_hong/results) 기본값으로 변경 — 단, 기존 예약 폴더인
#      ~/metagenomics_shotgun 및 /mnt/f/Metagenomics/Shotgun/results 는
#      "다른 분석"용으로 보존, 이 스크립트는 그 두 곳을 절대 읽지도 쓰지도 않음
#      (아래 safety_check_workspace()에서 강제 검증)
#      MEGAHIT 조립의 임시 파일만 여전히 BASE_DIR(로컬 디스크)에서 처리 후
#      F 드라이브로 이동 — WSL DrvFs FIFO 에러 방지 목적(원본 스크립트와 동일 기법)
#   2. FASTQ_DIR 기본값을 raw data 원본 경로로 직접 지정
#      (복사 없이 원본을 읽기 전용으로 참조 — 136GB 복사 시간/용량 절약)
#   3. Prokka(Perl XML::Simple)·eggNOG-mapper(diamond 경로) 두 가지 기존에
#      실패하던 이슈에 대한 사전 점검 + 자동 완화 로직 추가
#   4. 샘플 필터(SAMPLE_INCLUDE_REGEX) 추가 — 32개 전체 실행 전 1~2개로 파일럿 테스트 가능
#   5. 예상 디스크 사용량 사전 경고 추가
#   6. R1/R2 페어링 로직은 원본 그대로 유지 (2026-09-02 사용자 확인: "_1"/"_2" = R1/R2 맞음,
#      기존 basename 기반 페어링 로직에 버그 없음)
#
#  Pipeline Steps: 원본과 동일 (fastp → Bowtie2 host removal → MEGAHIT → QUAST →
#  MetaBAT2 → CheckM → Prokka → Kraken2/Bracken → Krona → DIAMOND(UniRef90) →
#  GTDB-Tk → eggNOG-mapper → HUMAnN3 → AMRFinderPlus → MultiQC)
#
#  실행 예시:
#   # 1) 파일럿 테스트 (권장, 샘플 1개만)
#   SAMPLE_INCLUDE_REGEX='^0M_1_p' bash 03_run_shotgun_meta_hong.sh 2>&1 | tee pilot_run.log
#
#   # 2) 문제 없으면 전체 32샘플 실행 (시간이 매우 오래 걸림 — nohup/tmux 권장)
#   nohup bash 03_run_shotgun_meta_hong.sh > full_run.log 2>&1 &
# ==============================================================================

set -euo pipefail
# set -x  # 디버깅 시에만 활성화

# ------------------------------------------------------------------------------
# 0. 색상 출력 헬퍼
# ------------------------------------------------------------------------------
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; MAGENTA='\033[0;35m'; RESET='\033[0m'

log_info()  { echo -e "${CYAN}[INFO]${RESET}  $*"; }
log_ok()    { echo -e "${GREEN}[OK]${RESET}    $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${RESET}  $*" >&2; }
log_error() { echo -e "${RED}[ERROR]${RESET} $*" >&2; }
log_step()  { echo -e "\n${MAGENTA}>>> $*${RESET}"; }

# ------------------------------------------------------------------------------
# 1. 사용자 설정 (환경에 맞게 수정 가능, 대부분 env var로 오버라이드 가능)
# ------------------------------------------------------------------------------

# --- 기본 디렉토리: 작업 기준은 meta_hong(로컬), 대용량 결과물은 F 드라이브 ---
BASE_DIR="${BASE_DIR:-$HOME/meta_hong}"
# raw FASTQ 원본을 그대로 참조 (복사하지 않음, 읽기 전용)
FASTQ_DIR="${FASTQ_DIR:-./raw_fastq}"
# 결과물(용량 큼)은 F 드라이브의 신규 전용 폴더에 기록
# (기존 예약 경로 /mnt/f/Metagenomics/Shotgun/results 와는 별개의 새 경로)
RESULTS_DIR="${RESULTS_DIR:-./results}"

mkdir -p "${BASE_DIR}"

# --- 안전장치: 기존 예약 폴더(홈 디렉토리 쪽 + F드라이브 쪽 둘 다)를 절대 건드리지 않도록 강제 검증 ---
FORBIDDEN_DIRS=(
    "$(realpath -m "$HOME/metagenomics_shotgun")"
    "$(realpath -m "/mnt/f/Metagenomics/Shotgun/results")"
)
safety_check_workspace() {
    local target_real forbidden
    for target in "${RESULTS_DIR}" "${BASE_DIR}"; do
        target_real="$(realpath -m "${target}")"
        for forbidden in "${FORBIDDEN_DIRS[@]}"; do
            if [[ "${target_real}" == "${forbidden}" || "${target_real}" == "${forbidden}/"* ]]; then
                log_error "안전장치 작동: '${target}' 경로가 기존 예약 폴더(${forbidden}) 하위입니다."
                log_error "이 폴더는 다른 분석용으로 보존하기로 되어 있어 이 스크립트는 절대 쓰지 않습니다."
                log_error "RESULTS_DIR / BASE_DIR을 /mnt/f/meta_hong/results 등 다른 신규 경로로 지정해 주세요."
                exit 1
            fi
        done
    done
    log_ok "작업공간 안전장치 통과 (예약 폴더 미사용 확인)"
}
safety_check_workspace

mkdir -p "$RESULTS_DIR"

# --- 스레드 ---
THREADS="${THREADS:-24}"

# --- Conda 설치 루트 (eggNOG diamond PATH 보정용) ---
CONDA_ROOT="${CONDA_ROOT:-$HOME/apps/anaconda}"

# --- Conda 환경 이름 ---
ENV_QC="metagenomics_shotgun_qc"
ENV_ASSEMBLY="metagenomics_shotgun_assembly"
ENV_REPORTING="metagenomics_shotgun_reporting"
ENV_QUAST="metagenomics_shotgun_quast"
ENV_PROKKA="metagenomics_shotgun_prokka"
ENV_KRAKEN="metagenomics_shotgun_kraken"
ENV_DIAMOND="metagenomics_shotgun_diamond"
ENV_METABAT="metagenomics_shotgun_metabat"
ENV_CHECKM="metagenomics_shotgun_checkm"
ENV_GTDBTK="metagenomics_shotgun_gtdbtk"
ENV_EGGNOG="metagenomics_shotgun_eggnog"
ENV_HUMANN="metagenomics_shotgun_humann"
ENV_AMRFINDER="metagenomics_shotgun_amrfinder"

# eggNOG-mapper 환경에 diamond가 없을 때 대비, DIAMOND 환경의 bin 디렉토리를
# PATH에 추가로 넘겨줄 수 있도록 미리 계산해 둔다 (Step 10에서 사용)
DIAMOND_BIN_DIR="${CONDA_ROOT}/envs/${ENV_DIAMOND}/bin"

# --- DB 경로 (환경변수로 오버라이드 가능) ---
DB_ROOT="${DB_ROOT:-./databases}"

KRAKEN2_DB_PATH="${KRAKEN2_DB_PATH:-${DB_ROOT}/Kraken_DB}"
DIAMOND_DB_UNIREF90="${DIAMOND_DB_UNIREF90:-${DB_ROOT}/Diamond_DB/uniref90.dmnd}"
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-${DB_ROOT}/Eggnog_DB}"
HUMANN_NUCLEOTIDE_DB_PATH="${HUMANN_NUCLEOTIDE_DB_PATH:-${DB_ROOT}/Humann_chocophlan_DB/chocophlan}"
HUMANN_PROTEIN_DB_PATH="${HUMANN_PROTEIN_DB_PATH:-${DB_ROOT}/Humann_uniref_DB/uniref}"
HUMANN_TMP_DIR="${HUMANN_TMP_DIR:-/tmp/humann}"

export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-${DB_ROOT}/Gtdbtk_DB}"
export AMRFINDER_DB="${AMRFINDER_DB:-${DB_ROOT}/Amrfinderplus_DB/data/2025-03-25.1}"

# --- 호스트 게놈 목록 ---
HOSTS=(
    "glycine_max"
    "homo_sapiens"
    "hordeum_vulgare"
    "oryza_sativa"
    "triticum_aestivum"
)
HOST_GENOME_ROOT="${HOST_GENOME_ROOT:-${DB_ROOT}/Host_genome}"

# --- 분석 단계 ON/OFF 플래그 ---
RUN_HOST_REMOVAL="${RUN_HOST_REMOVAL:-true}"
RUN_GTDBTK="${RUN_GTDBTK:-true}"
RUN_EGGNOG_MAPPER="${RUN_EGGNOG_MAPPER:-true}"
RUN_HUMANN="${RUN_HUMANN:-true}"
RUN_AMRFINDERPLUS="${RUN_AMRFINDERPLUS:-true}"
RUN_KRONA_KRAKEN="${RUN_KRONA_KRAKEN:-true}"

# --- Prokka 옵션 ---
PROKKA_KINGDOM="${PROKKA_KINGDOM:-Bacteria}"
PROKKA_METAGENOME_MODE="${PROKKA_METAGENOME_MODE:-false}"

# --- fastp 파라미터 ---
FASTP_QUAL="${FASTP_QUAL:-20}"
FASTP_UNQUALIFIED="${FASTP_UNQUALIFIED:-20}"
FASTP_MIN_LEN="${FASTP_MIN_LEN:-50}"

# --- Kraken2/Bracken 파라미터 ---
BRACKEN_READ_LEN="${BRACKEN_READ_LEN:-100}"

# --- CheckM 필터링 기준 ---
CHECKM_MIN_COMPLETENESS="${CHECKM_MIN_COMPLETENESS:-50}"
CHECKM_MAX_CONTAMINATION="${CHECKM_MAX_CONTAMINATION:-10}"
CHECKM_HQ_COMPLETENESS="${CHECKM_HQ_COMPLETENESS:-90}"
CHECKM_HQ_CONTAMINATION="${CHECKM_HQ_CONTAMINATION:-5}"

# --- MAG 필터링 기준 ---
MAG_MIN_CONTIGS="${MAG_MIN_CONTIGS:-2}"
MAG_MIN_SIZE_BP="${MAG_MIN_SIZE_BP:-200000}"

# --- 샘플 필터 (파일럿 테스트용, 정규식. 비어 있으면 전체 샘플 실행) ---
# 예: SAMPLE_INCLUDE_REGEX='^0M_1_p' 로 지정하면 0M_1_p 샘플 하나만 실행
SAMPLE_INCLUDE_REGEX="${SAMPLE_INCLUDE_REGEX:-}"

# --- Summary 파일 ---
SUMMARY_FILE="${RESULTS_DIR}/summary_stats.md"

# ==============================================================================
# 2. 사전 검증
# ==============================================================================

validate_environment() {
    log_step "사전 환경 검증 시작..."

    local errors=0

    if [ ! -d "${FASTQ_DIR}" ]; then
        log_error "FASTQ 디렉토리가 없습니다: ${FASTQ_DIR}"
        (( errors++ )) || true
    fi

    if [ ! -d "${GTDBTK_DATA_PATH}" ]; then
        log_error "GTDBTK_DATA_PATH 디렉토리 없음: ${GTDBTK_DATA_PATH}"
        (( errors++ )) || true
    fi

    if [ ! -d "${AMRFINDER_DB}" ]; then
        log_error "AMRFINDER_DB 디렉토리 없음: ${AMRFINDER_DB}"
        (( errors++ )) || true
    fi

    if [ ! -f "${DIAMOND_DB_UNIREF90}" ]; then
        log_error "DIAMOND UniRef90 DB 없음: ${DIAMOND_DB_UNIREF90}"
        (( errors++ )) || true
    fi

    if [ ! -d "${HUMANN_NUCLEOTIDE_DB_PATH}" ]; then
        log_error "HUMAnN3 nucleotide DB 없음: ${HUMANN_NUCLEOTIDE_DB_PATH}"
        (( errors++ )) || true
    fi

    if [ ! -d "${HUMANN_PROTEIN_DB_PATH}" ]; then
        log_error "HUMAnN3 protein DB 없음: ${HUMANN_PROTEIN_DB_PATH}"
        (( errors++ )) || true
    fi

    if [ ! -d "${KRAKEN2_DB_PATH}" ]; then
        log_warn "Kraken2 DB 없음: ${KRAKEN2_DB_PATH} (Step 5 건너뜀)"
    fi

    # --- 알려진 이슈 사전 점검 (기존 metagenomics_shotgun 파이프라인에서 실패했던 항목) ---
    log_info "알려진 이슈 사전 점검 중 (Prokka Perl 모듈 / eggNOG diamond 경로)..."

    if ! conda run -n "${ENV_PROKKA}" perl -e 'use XML::Simple;' >/dev/null 2>&1; then
        log_warn "Prokka 환경(${ENV_PROKKA})에 Perl 모듈 XML::Simple 이 없습니다."
        log_warn "  → Step 4(Prokka)가 실패할 수 있습니다. 해결 명령 예시:"
        log_warn "     conda run -n ${ENV_PROKKA} cpanm XML::Simple"
        log_warn "     (또는) conda install -n ${ENV_PROKKA} -c bioconda perl-xml-simple"
        log_warn "  Prokka가 실패해도 파이프라인은 계속 진행되며, 이후 Prokka 결과에 의존하는"
        log_warn "  DIAMOND/eggNOG/AMRFinderPlus 단계만 자동으로 건너뜁니다."
    else
        log_ok "Prokka Perl 모듈(XML::Simple) 확인됨"
    fi

    if conda run -n "${ENV_EGGNOG}" which diamond >/dev/null 2>&1; then
        log_ok "eggNOG-mapper 환경에 diamond 바이너리 존재 확인됨"
    elif [ -x "${DIAMOND_BIN_DIR}/diamond" ]; then
        log_warn "eggNOG-mapper 환경(${ENV_EGGNOG})에 diamond가 없어,"
        log_warn "  DIAMOND 환경(${ENV_DIAMOND})의 diamond 바이너리를 PATH로 임시 연결해 사용합니다."
        log_warn "  (${DIAMOND_BIN_DIR}/diamond) — Step 10에서 자동 적용됨"
    else
        log_warn "diamond 바이너리를 어디에서도 찾지 못했습니다 (${ENV_EGGNOG}, ${ENV_DIAMOND} 둘 다 없음)."
        log_warn "  → Step 10(eggNOG-mapper)이 실패할 수 있습니다."
    fi

    if [ "${errors}" -gt 0 ]; then
        log_error "검증 실패: ${errors}개 오류. 위 항목을 확인하고 다시 실행하세요."
        exit 1
    fi

    log_ok "모든 필수 DB 검증 완료"
}

# --- 예상 디스크 사용량 경고 (파일럿 실행 결과 기반 추정치) ---
warn_disk_estimate() {
    log_step "예상 디스크 사용량 확인..."
    local n_samples avail_kb avail_gb est_gb
    n_samples=$(find "${FASTQ_DIR}" -name "*_1.fastq.gz" 2>/dev/null | wc -l)
    # 기존 1~2개 샘플 파일럿 실행 결과(약 41GB/2샘플 ≈ 20GB/샘플)를 기준으로 대략 추정
    est_gb=$(( n_samples * 20 ))
    avail_kb=$(df -Pk "${RESULTS_DIR}" 2>/dev/null | awk 'NR==2{print $4}')
    avail_gb=$(( avail_kb / 1024 / 1024 ))
    log_info "발견된 샘플 수: ${n_samples}개, 대략적 예상 총 용량: 약 ${est_gb}GB (샘플당 ~20GB 가정, 실제는 다를 수 있음)"
    log_info "RESULTS_DIR(${RESULTS_DIR}) 소속 파티션 여유 공간: 약 ${avail_gb}GB"
    if [ "${avail_gb}" -lt "${est_gb}" ]; then
        log_warn "여유 공간이 예상 사용량보다 적을 수 있습니다. RESULTS_DIR을 더 큰 드라이브로 지정하는 것을 고려하세요."
        log_warn "  예: RESULTS_DIR=/mnt/f/meta_hong_results bash $0"
    fi
}

# ==============================================================================
# 3. 헬퍼 함수
# ==============================================================================

bowtie2_index_exists() {
    local prefix="$1"
    [ -f "${prefix}.1.bt2" ] || [ -f "${prefix}.1.bt2l" ]
}

bins_exist() {
    local dir="$1"
    [ -d "${dir}" ] && [ -n "$(find "${dir}" -maxdepth 1 -name '*.fa' -print -quit 2>/dev/null)" ]
}

init_summary() {
    mkdir -p "$(dirname "${SUMMARY_FILE}")"
    if [ ! -f "${SUMMARY_FILE}" ]; then
        cat > "${SUMMARY_FILE}" <<'EOF'
| Sample | Total Reads | Trimmed Reads | N50 (bp) | Contigs | HQ MAGs |
|--------|-------------|---------------|----------|---------|---------|
EOF
    fi
}

append_summary() {
    local sample="$1" total="$2" trimmed="$3" n50="$4" contigs="$5" hq_mags="$6"
    echo "| ${sample} | ${total:-N/A} | ${trimmed:-N/A} | ${n50:-N/A} | ${contigs:-N/A} | ${hq_mags:-N/A} |" >> "${SUMMARY_FILE}"
}

count_reads() {
    local fq="$1"
    if [[ "${fq}" == *.gz ]]; then
        conda run -n "${ENV_QC}" pigz -dc "${fq}" | awk 'NR%4==1{c++} END{print c}'
    else
        awk 'NR%4==1{c++} END{print c}' "${fq}"
    fi
}

# ==============================================================================
# 4. Host Bowtie2 인덱스 빌드 (샘플 루프 전에 1회 실행)
# ==============================================================================

build_host_indices() {
    log_step "Host Bowtie2 인덱스 확인/빌드..."

    for HOST in "${HOSTS[@]}"; do
        local HOST_GENOME_DIR="${HOST_GENOME_ROOT}/${HOST}"
        local HOST_GENOME_FASTA
        HOST_GENOME_FASTA=$(find "${HOST_GENOME_DIR}" -name "*.fna" 2>/dev/null | head -n1)
        local HOST_GENOME_BOWTIE2_IDX="${HOST_GENOME_DIR}/bowtie2_index/${HOST}"

        if [ -z "${HOST_GENOME_FASTA}" ]; then
            log_warn "${HOST}: FASTA 파일 없음 (${HOST_GENOME_DIR}), 건너뜀"
            continue
        fi

        if bowtie2_index_exists "${HOST_GENOME_BOWTIE2_IDX}"; then
            log_info "${HOST}: 인덱스 이미 존재, 건너뜀"
            continue
        fi

        log_info "${HOST}: Bowtie2 인덱스 빌드 중..."
        mkdir -p "$(dirname "${HOST_GENOME_BOWTIE2_IDX}")"

        if conda run -n "${ENV_ASSEMBLY}" bowtie2-build \
                --threads "${THREADS}" \
                "${HOST_GENOME_FASTA}" \
                "${HOST_GENOME_BOWTIE2_IDX}"; then
            log_ok "${HOST}: 인덱스 빌드 완료"
        else
            log_warn "${HOST}: 인덱스 빌드 실패 (계속 진행)"
        fi
    done
}

# ==============================================================================
# 5. 샘플별 분석 함수 (원본 run_shotgun_mnt.sh와 동일 로직)
# ==============================================================================

process_sample() {
    local R1_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${R1_FILE}" _1.fastq.gz)
    local R2_FILE
    R2_FILE=$(dirname "${R1_FILE}")/${SAMPLE_NAME}_2.fastq.gz

    if [ ! -f "${R2_FILE}" ]; then
        log_warn "R2 파일 없음: ${R2_FILE} → 샘플 ${SAMPLE_NAME} 건너뜀"
        return 0
    fi

    echo ""
    echo -e "${MAGENTA}============================================================${RESET}"
    echo -e "${MAGENTA}  샘플 처리: ${SAMPLE_NAME}${RESET}"
    echo -e "${MAGENTA}============================================================${RESET}"

    local SAMPLE_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
    local QC_DIR="${SAMPLE_DIR}/01_qc_and_trimmed"
    local HOST_DIR="${SAMPLE_DIR}/02_host_removed"
    local ASSEMBLY_DIR="${SAMPLE_DIR}/03_assembly"
    local ANNOTATION_DIR="${SAMPLE_DIR}/04_annotation_prokka"
    local TAXONOMY_DIR="${SAMPLE_DIR}/05_taxonomy_reads_kraken"
    local DIAMOND_DIR="${SAMPLE_DIR}/06_function_diamond_uniref90"
    local BINNING_DIR="${SAMPLE_DIR}/03_assembly"
    local MAG_QC_DIR="${SAMPLE_DIR}/03_assembly"
    local MAG_TAXONOMY_DIR="${SAMPLE_DIR}/09_mag_taxonomy_gtdbtk"
    local EGGNOG_DIR="${SAMPLE_DIR}/10_function_eggnog"
    local HUMANN_DIR="${SAMPLE_DIR}/11_function_community_humann"
    local AMRFINDER_DIR="${SAMPLE_DIR}/12_resistance_virulence_amrfinder"

    mkdir -p "${QC_DIR}" "${HOST_DIR}" "${ASSEMBLY_DIR}" \
             "${ANNOTATION_DIR}" "${TAXONOMY_DIR}" "${DIAMOND_DIR}" \
             "${BINNING_DIR}" "${MAG_QC_DIR}" "${MAG_TAXONOMY_DIR}" \
             "${EGGNOG_DIR}" "${HUMANN_DIR}" "${AMRFINDER_DIR}"

    local CURRENT_R1="${R1_FILE}"
    local CURRENT_R2="${R2_FILE}"
    local TOTAL_READS="" TRIMMED_READS="" N50="" CONTIG_COUNT="" HQ_MAG_COUNT=""

    # ------------------------------------------------------------------
    # Step 1. QC & Trimming (fastp)
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 1: fastp QC/Trimming (env: ${ENV_QC})"

    local FASTP_R1="${QC_DIR}/${SAMPLE_NAME}_1.trimmed.fastq.gz"
    local FASTP_R2="${QC_DIR}/${SAMPLE_NAME}_2.trimmed.fastq.gz"
    local FASTP_HTML="${QC_DIR}/${SAMPLE_NAME}.fastp.html"
    local FASTP_JSON="${QC_DIR}/${SAMPLE_NAME}.fastp.json"

    if [ ! -f "${FASTP_R1}" ] || [ ! -f "${FASTP_R2}" ]; then
        TOTAL_READS=$(count_reads "${CURRENT_R1}")
        conda run -n "${ENV_QC}" fastp \
            -i "${CURRENT_R1}" -I "${CURRENT_R2}" \
            -o "${FASTP_R1}" -O "${FASTP_R2}" \
            --html "${FASTP_HTML}" \
            --json "${FASTP_JSON}" \
            --thread 16 \
            -q "${FASTP_QUAL}" \
            -u "${FASTP_UNQUALIFIED}" \
            -l "${FASTP_MIN_LEN}"
        log_ok "[${SAMPLE_NAME}] Step 1: 완료"
    else
        log_info "[${SAMPLE_NAME}] Step 1: 결과 존재, 건너뜀"
    fi

    CURRENT_R1="${FASTP_R1}"
    CURRENT_R2="${FASTP_R2}"

    if [ -f "${FASTP_JSON}" ]; then
        TRIMMED_READS=$(python3 -c "
import json, sys
with open('${FASTP_JSON}') as f:
    d = json.load(f)
print(d.get('filtering_result', {}).get('passed_filter_reads', 'N/A'))
" 2>/dev/null || echo "N/A")
    fi

    # ------------------------------------------------------------------
    # Step 2. Host DNA Removal (Bowtie2 - 독립 매핑 후 합집합 방식)
    # ------------------------------------------------------------------
    if [ "${RUN_HOST_REMOVAL}" = true ]; then
        log_step "[${SAMPLE_NAME}] Step 2: Host DNA Removal - 독립 매핑 합집합 방식 (env: ${ENV_ASSEMBLY})"

        local FINAL_HR_R1="${HOST_DIR}/${SAMPLE_NAME}_1.hostremoved.fastq.gz"
        local FINAL_HR_R2="${HOST_DIR}/${SAMPLE_NAME}_2.hostremoved.fastq.gz"
        local COMBINED_HOST_READIDS="${HOST_DIR}/${SAMPLE_NAME}.all_host_readids.txt"

        if [ -f "${FINAL_HR_R1}" ] && [ -f "${FINAL_HR_R2}" ]; then
            log_info "[${SAMPLE_NAME}] Step 2: 최종 결과 존재, 건너뜀"
            CURRENT_R1="${FINAL_HR_R1}"
            CURRENT_R2="${FINAL_HR_R2}"
        else
            > "${COMBINED_HOST_READIDS}"

            for HOST in "${HOSTS[@]}"; do
                local HOST_IDX="${HOST_GENOME_ROOT}/${HOST}/bowtie2_index/${HOST}"
                if ! bowtie2_index_exists "${HOST_IDX}"; then
                    log_warn "[${SAMPLE_NAME}] ${HOST} 인덱스 없음, 건너뜀"
                    continue
                fi

                local BT2_LOG="${HOST_DIR}/${SAMPLE_NAME}.${HOST}.bowtie2.log"
                local BT2_METRICS="${HOST_DIR}/${SAMPLE_NAME}.${HOST}.bowtie2.metrics.txt"
                local HOST_BAM="${HOST_DIR}/${SAMPLE_NAME}.${HOST}.host_mapping.bam"

                log_info "[${SAMPLE_NAME}] ${HOST} 독립 매핑 중..."

                conda run -n "${ENV_ASSEMBLY}" bash -c \
                    "bowtie2 \
                        -x '${HOST_IDX}' \
                        -1 '${CURRENT_R1}' \
                        -2 '${CURRENT_R2}' \
                        --threads '${THREADS}' \
                        --very-sensitive-local \
                        --no-unal \
                        --met-file '${BT2_METRICS}' \
                        -S /dev/stdout \
                        2> '${BT2_LOG}' \
                     | samtools view -@ '${THREADS}' -bS -F 4 - \
                     > '${HOST_BAM}'"

                if [ -s "${HOST_BAM}" ]; then
                    conda run -n "${ENV_ASSEMBLY}" samtools view \
                        -@ "${THREADS}" "${HOST_BAM}" \
                        | awk '{gsub(/\/[12]$/, "", $1); print $1}' \
                        >> "${COMBINED_HOST_READIDS}"
                    log_info "[${SAMPLE_NAME}] ${HOST}: mapped reads 누적 완료"
                else
                    log_warn "[${SAMPLE_NAME}] ${HOST}: mapped BAM 비어 있음 (호스트 오염 없음)"
                fi

                python3 - <<PYEOF
import json, re
summary = {}
try:
    with open("${BT2_LOG}") as f:
        for line in f:
            line = line.strip()
            if re.search(r"reads; of these:", line):
                summary["total_reads"] = int(line.split()[0])
            elif "aligned concordantly 0 times" in line:
                summary["aligned_0"] = int(line.split()[0])
            elif "aligned concordantly exactly 1 time" in line:
                summary["aligned_1"] = int(line.split()[0])
            elif "aligned concordantly >1 times" in line:
                summary["aligned_gt1"] = int(line.split()[0])
            elif "overall alignment rate" in line:
                m = re.search(r"([\d.]+)%", line)
                if m:
                    summary["overall_alignment_rate"] = float(m.group(1))
    out_json = "${BT2_LOG%.log}.json"
    with open(out_json, 'w') as out:
        json.dump(summary, out, indent=4)
    print(f"[INFO] ${HOST} Bowtie2 summary → {out_json}")
except Exception as e:
    print(f"[WARN] ${HOST} JSON 생성 실패: {e}")
PYEOF

                rm -f "${HOST_BAM}"

            done

            sort -u "${COMBINED_HOST_READIDS}" -o "${COMBINED_HOST_READIDS}"
            local TOTAL_HOST_READS
            TOTAL_HOST_READS=$(wc -l < "${COMBINED_HOST_READIDS}")
            log_info "[${SAMPLE_NAME}] 전체 호스트 매핑 read ID 합계: ${TOTAL_HOST_READS}개 (중복 제거 후)"

            log_info "[${SAMPLE_NAME}] 호스트 read 최종 제거 중 (filterbyname)..."
            conda run -n "${ENV_ASSEMBLY}" filterbyname.sh \
                in="${CURRENT_R1}" \
                in2="${CURRENT_R2}" \
                out="${FINAL_HR_R1}" \
                out2="${FINAL_HR_R2}" \
                names="${COMBINED_HOST_READIDS}" \
                include=f \
                threads="${THREADS}" \
                -Xmx16g

            rm -f "${COMBINED_HOST_READIDS}"

            if [ -f "${FINAL_HR_R1}" ] && [ -f "${FINAL_HR_R2}" ]; then
                log_ok "[${SAMPLE_NAME}] Step 2: 호스트 제거 완료 → ${FINAL_HR_R1}"
            else
                log_error "[${SAMPLE_NAME}] Step 2: 최종 호스트 제거 파일 생성 실패"
                exit 1
            fi

            CURRENT_R1="${FINAL_HR_R1}"
            CURRENT_R2="${FINAL_HR_R2}"
        fi
        log_ok "[${SAMPLE_NAME}] Step 2: 완료"
    fi

    # ------------------------------------------------------------------
    # Step 3. Assembly (MEGAHIT)
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 3: MEGAHIT Assembly (env: ${ENV_ASSEMBLY})"

    local MEGAHIT_DIR="${ASSEMBLY_DIR}/megahit_output"
    local FINAL_CONTIGS="${MEGAHIT_DIR}/final.contigs.fa"
    local LOCAL_TMP_DIR="${BASE_DIR}/tmp_assembly_${SAMPLE_NAME}"

    if [ ! -f "${FINAL_CONTIGS}" ]; then
        rm -rf "${MEGAHIT_DIR}"
        rm -rf "${LOCAL_TMP_DIR}"

        conda run -n "${ENV_ASSEMBLY}" megahit \
            -1 "${CURRENT_R1}" \
            -2 "${CURRENT_R2}" \
            -o "${LOCAL_TMP_DIR}" \
            -t "${THREADS}" \
            --min-contig-len 1000

        if [ -d "${LOCAL_TMP_DIR}" ]; then
            mkdir -p "$(dirname "${MEGAHIT_DIR}")"
            mv "${LOCAL_TMP_DIR}" "${MEGAHIT_DIR}"
            log_ok "[${SAMPLE_NAME}] Step 3: Assembly 완료 및 결과 이동 성공"
        else
            log_warn "[${SAMPLE_NAME}] Step 3: MEGAHIT 실행 실패 (결과 디렉토리 없음)"
            exit 1
        fi
    else
        log_info "[${SAMPLE_NAME}] Step 3: 결과 존재, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 3.1. QUAST
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 3.1: QUAST (env: ${ENV_QUAST})"

    local QUAST_DIR="${ASSEMBLY_DIR}/quast"
    local QUAST_REPORT="${QUAST_DIR}/report.tsv"

    if [ ! -f "${QUAST_REPORT}" ]; then
        mkdir -p "${QUAST_DIR}"
        conda run -n "${ENV_QUAST}" quast.py \
            -o "${QUAST_DIR}" \
            "${FINAL_CONTIGS}" \
            --threads "${THREADS}" \
            --silent
        log_ok "[${SAMPLE_NAME}] Step 3.1: QUAST 완료"
    else
        log_info "[${SAMPLE_NAME}] Step 3.1: 결과 존재, 건너뜀"
    fi

    if [ -f "${QUAST_REPORT}" ]; then
        N50=$(awk -F'\t' '$1=="N50" {print $2}' "${QUAST_REPORT}" 2>/dev/null || echo "N/A")
        CONTIG_COUNT=$(awk -F'\t' '$1=="# contigs (>= 0 bp)" {print $2}' "${QUAST_REPORT}" 2>/dev/null || echo "N/A")
    fi

    # ------------------------------------------------------------------
    # Step 3.2. MAG Binning (MetaBAT2)
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 3.2: MetaBAT2 Binning (env: ${ENV_METABAT})"

    local BINS_DIR="${BINNING_DIR}/bins"
    local CONTIG_IDX="${MEGAHIT_DIR}/final.contigs"
    local BAM_PREFIX="${BINNING_DIR}/${SAMPLE_NAME}"
    local SORTED_BAM="${BAM_PREFIX}.sorted.bam"
    local DEPTH_FILE="${BINNING_DIR}/${SAMPLE_NAME}.depth.txt"
    local CLEANED_DEPTH="${BINNING_DIR}/${SAMPLE_NAME}.depth.cleaned.txt"

    if [ -s "${FINAL_CONTIGS}" ] && ! bins_exist "${BINS_DIR}"; then
        mkdir -p "${BINS_DIR}"

        if ! bowtie2_index_exists "${CONTIG_IDX}"; then
            conda run -n "${ENV_ASSEMBLY}" bowtie2-build \
                --threads "${THREADS}" \
                "${FINAL_CONTIGS}" \
                "${CONTIG_IDX}"
        fi

        conda run -n "${ENV_ASSEMBLY}" bash -c \
            "bowtie2 -x \"${CONTIG_IDX}\" \
                -1 \"${CURRENT_R1}\" -2 \"${CURRENT_R2}\" \
                --no-unal -p \"${THREADS}\" -S /dev/stdout \
             | samtools view -@ \"${THREADS}\" -bS - \
             > \"${BAM_PREFIX}.bam\""

        conda run -n "${ENV_ASSEMBLY}" samtools sort \
            -@ "${THREADS}" -o "${SORTED_BAM}" "${BAM_PREFIX}.bam"
        rm -f "${BAM_PREFIX}.bam"
        conda run -n "${ENV_ASSEMBLY}" samtools index \
            -@ "${THREADS}" "${SORTED_BAM}"

        conda run -n "${ENV_METABAT}" jgi_summarize_bam_contig_depths \
            --outputDepth "${DEPTH_FILE}" "${SORTED_BAM}"

        awk 'NR==1 || ($3 >= 0 && $3 < 1e6)' "${DEPTH_FILE}" > "${CLEANED_DEPTH}"

        conda run -n "${ENV_METABAT}" metabat2 \
            -i "${FINAL_CONTIGS}" \
            -a "${CLEANED_DEPTH}" \
            -o "${BINS_DIR}/${SAMPLE_NAME}.bin" \
            -t "${THREADS}" \
            -m 1500

        log_ok "[${SAMPLE_NAME}] Step 3.2: Binning 완료"
    else
        log_info "[${SAMPLE_NAME}] Step 3.2: 결과 존재 또는 contig 없음, 건너뜀"
    fi

    local FILTERED_BINS_DIR="${BINS_DIR}/filtered"
    if bins_exist "${BINS_DIR}" && [ ! -d "${FILTERED_BINS_DIR}" ]; then
        log_info "[${SAMPLE_NAME}] 저품질 bin 필터링 (contig < ${MAG_MIN_CONTIGS} or size < ${MAG_MIN_SIZE_BP} bp)..."
        mkdir -p "${FILTERED_BINS_DIR}"
        for BIN_FA in "${BINS_DIR}"/*.fa; do
            [ -f "${BIN_FA}" ] || continue
            local BN
            BN=$(basename "${BIN_FA}")
            local CC BZ
            CC=$(grep -c '^>' "${BIN_FA}")
            BZ=$(awk '/^>/{if(s)t+=s; s=0; next}{s+=length($0)} END{t+=s; print t}' "${BIN_FA}")
            if [ "${CC}" -ge "${MAG_MIN_CONTIGS}" ] && [ "${BZ}" -ge "${MAG_MIN_SIZE_BP}" ]; then
                cp "${BIN_FA}" "${FILTERED_BINS_DIR}/${BN}"
            else
                log_info "  [필터됨] ${BN}: ${CC} contigs, ${BZ} bp"
            fi
        done
        log_ok "[${SAMPLE_NAME}] Bin 필터링 완료 → ${FILTERED_BINS_DIR}"
    fi

    # ------------------------------------------------------------------
    # Step 3.3. CheckM
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 3.3: CheckM (env: ${ENV_CHECKM})"

    local CHECKM_OUT="${MAG_QC_DIR}/checkm_output"
    local CHECKM_SUMMARY="${MAG_QC_DIR}/${SAMPLE_NAME}.checkm_summary.tsv"
    local CHECKM_QA="${MAG_QC_DIR}/${SAMPLE_NAME}.checkm_qa.tsv"
    local CHECKM_QA_CSV="${MAG_QC_DIR}/${SAMPLE_NAME}.checkm_qa.csv"
    local CHECKM_QA_XLSX="${MAG_QC_DIR}/${SAMPLE_NAME}.checkm_qa.xlsx"

    local BINS_FOR_CHECKM="${BINS_DIR}"
    bins_exist "${FILTERED_BINS_DIR}" && BINS_FOR_CHECKM="${FILTERED_BINS_DIR}"

    if bins_exist "${BINS_FOR_CHECKM}" && \
       ( [ ! -d "${CHECKM_OUT}" ] || [ -z "$(ls -A "${CHECKM_OUT}" 2>/dev/null)" ] ); then

        mkdir -p "${CHECKM_OUT}"
        log_info "[${SAMPLE_NAME}] CheckM lineage_wf 실행 중..."
        conda run -n "${ENV_CHECKM}" checkm lineage_wf \
            -x fa "${BINS_FOR_CHECKM}" \
            "${CHECKM_OUT}" \
            --threads "${THREADS}" \
            --pplacer_threads "${THREADS}" \
            -f "${CHECKM_SUMMARY}"
        log_ok "[${SAMPLE_NAME}] Step 3.3: CheckM lineage_wf 완료"
    else
        log_info "[${SAMPLE_NAME}] Step 3.3: CheckM 결과 존재 또는 bin 없음, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 3.4. CheckM QA 리포트 + CSV/XLSX 변환
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 3.4: CheckM QA 리포트 (env: ${ENV_CHECKM})"

    if [ -f "${CHECKM_OUT}/lineage.ms" ] && \
       ( [ ! -f "${CHECKM_QA}" ] || [ ! -s "${CHECKM_QA}" ] ); then
        conda run -n "${ENV_CHECKM}" checkm qa \
            "${CHECKM_OUT}/lineage.ms" \
            "${CHECKM_OUT}" \
            --out_format 2 \
            -f "${CHECKM_QA}" \
            --threads "${THREADS}"
        log_ok "[${SAMPLE_NAME}] Step 3.4: QA 리포트 생성 완료"
    elif [ ! -f "${CHECKM_OUT}/lineage.ms" ]; then
        log_warn "[${SAMPLE_NAME}] lineage.ms 없음, QA 건너뜀"
    fi

    if [ -f "${CHECKM_QA}" ] && [ ! -f "${CHECKM_QA_XLSX}" ]; then
        conda run -n "${ENV_REPORTING}" python3 - <<PYEOF
import pandas as pd, sys, traceback
try:
    df = pd.read_csv("${CHECKM_QA}", sep='\t', comment='-', skip_blank_lines=True)
    df.to_csv("${CHECKM_QA_CSV}", index=False)
    df.to_excel("${CHECKM_QA_XLSX}", index=False)
    print("[INFO] CheckM QA CSV/XLSX 저장 완료")
except Exception:
    traceback.print_exc()
    sys.exit(1)
PYEOF
    fi

    if [ -f "${CHECKM_QA}" ]; then
        HQ_MAG_COUNT=$(grep -v '^[-[:space:]]*$' "${CHECKM_QA}" \
                       | grep -v '^-' \
                       | awk -v comp="${CHECKM_HQ_COMPLETENESS}" \
                             -v cont="${CHECKM_HQ_CONTAMINATION}" \
                       'NR>1 && $7+0 >= comp && $8+0 <= cont {c++} END{print c+0}')
    fi
    log_info "[${SAMPLE_NAME}] HQ MAGs: ${HQ_MAG_COUNT:-N/A} (Completeness≥${CHECKM_HQ_COMPLETENESS}%, Contamination≤${CHECKM_HQ_CONTAMINATION}%)"

    # ------------------------------------------------------------------
    # Step 4. Prokka
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 4: Prokka Annotation (env: ${ENV_PROKKA})"

    local PROKKA_FAA="${ANNOTATION_DIR}/${SAMPLE_NAME}.faa"
    local PROKKA_GFF="${ANNOTATION_DIR}/${SAMPLE_NAME}.gff"

    if [ -s "${FINAL_CONTIGS}" ] && [ ! -f "${PROKKA_FAA}" ]; then
        local PROKKA_OPTS="--outdir ${ANNOTATION_DIR} --prefix ${SAMPLE_NAME} --cpus ${THREADS} --force"
        if [ "${PROKKA_METAGENOME_MODE}" = true ]; then
            PROKKA_OPTS="${PROKKA_OPTS} --metagenome"
        else
            PROKKA_OPTS="${PROKKA_OPTS} --kingdom ${PROKKA_KINGDOM}"
        fi
        # shellcheck disable=SC2086
        if ! conda run -n "${ENV_PROKKA}" prokka ${PROKKA_OPTS} "${FINAL_CONTIGS}"; then
            log_warn "[${SAMPLE_NAME}] Step 4: Prokka 실패 — XML::Simple 등 Perl 모듈 문제일 수 있습니다."
            log_warn "  → DIAMOND/eggNOG/AMRFinderPlus 단계는 이 샘플에서 자동으로 건너뜁니다."
        else
            log_ok "[${SAMPLE_NAME}] Step 4: Prokka 완료"
        fi
    else
        log_info "[${SAMPLE_NAME}] Step 4: 결과 존재 또는 contig 없음, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 5. Kraken2 + Bracken
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 5: Kraken2/Bracken Taxonomy (env: ${ENV_KRAKEN})"

    local KRAKEN_REPORT="${TAXONOMY_DIR}/${SAMPLE_NAME}.kraken2.report"
    local KRAKEN_OUT="${TAXONOMY_DIR}/${SAMPLE_NAME}.kraken2.out"

    if [ -d "${KRAKEN2_DB_PATH}" ] && [ ! -f "${KRAKEN_REPORT}" ]; then
        conda run -n "${ENV_KRAKEN}" kraken2 \
            --db "${KRAKEN2_DB_PATH}" \
            --threads "${THREADS}" \
            --paired "${CURRENT_R1}" "${CURRENT_R2}" \
            --report "${KRAKEN_REPORT}" \
            --output "${KRAKEN_OUT}" \
            --gzip-compressed

        local BRACKEN_KMER_DISTRIB="${KRAKEN2_DB_PATH}/database${BRACKEN_READ_LEN}mers.kmer_distrib"
        if [ -f "${BRACKEN_KMER_DISTRIB}" ]; then
            conda run -n "${ENV_KRAKEN}" bracken \
                -d "${KRAKEN2_DB_PATH}" \
                -i "${KRAKEN_REPORT}" \
                -o "${TAXONOMY_DIR}/${SAMPLE_NAME}.bracken_S.txt" \
                -r "${BRACKEN_READ_LEN}" \
                -l S
        else
            log_warn "[${SAMPLE_NAME}] Bracken kmer_distrib 없음: ${BRACKEN_KMER_DISTRIB}"
        fi
        log_ok "[${SAMPLE_NAME}] Step 5: Kraken2/Bracken 완료"
    else
        log_info "[${SAMPLE_NAME}] Step 5: 결과 존재 또는 DB 없음, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 5.1/5.2. Krona
    # ------------------------------------------------------------------
    if [ "${RUN_KRONA_KRAKEN}" = true ] && [ -s "${KRAKEN_REPORT}" ]; then
        log_step "[${SAMPLE_NAME}] Step 5.1/5.2: Krona 시각화 (env: ${ENV_KRAKEN})"

        local KRONA_TEXT_HTML="${TAXONOMY_DIR}/${SAMPLE_NAME}.kraken2.krona.text.html"
        local KRONA_TAX_HTML="${TAXONOMY_DIR}/${SAMPLE_NAME}.kraken2.krona.tax.html"
        local KRONA_TXT="${TAXONOMY_DIR}/${SAMPLE_NAME}.kraken_for_krona.txt"

        if [ ! -f "${KRONA_TEXT_HTML}" ]; then
            awk -F'\t' '$4 == "U" || $4 == "R" || $4 == "D" || $4 == "P" || \
                        $4 == "C" || $4 == "O" || $4 == "F" || $4 == "G" || $4 == "S" \
                        {print $2 "\t" $6}' "${KRAKEN_REPORT}" \
                | sed 's/^ *//; s/ /_/g' > "${KRONA_TXT}"
            conda run -n "${ENV_KRAKEN}" ktImportText "${KRONA_TXT}" -o "${KRONA_TEXT_HTML}"
        fi

        if [ ! -f "${KRONA_TAX_HTML}" ]; then
            conda run -n "${ENV_KRAKEN}" ktImportTaxonomy \
                "${KRAKEN_REPORT}" -o "${KRONA_TAX_HTML}"
        fi
        log_ok "[${SAMPLE_NAME}] Step 5.1/5.2: Krona 완료"
    fi

    # ------------------------------------------------------------------
    # Step 6. DIAMOND vs UniRef90
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 6: DIAMOND vs UniRef90 (env: ${ENV_DIAMOND})"

    local DIAMOND_OUT="${DIAMOND_DIR}/${SAMPLE_NAME}.prokka.diamond_uniref90.paf"

    if [ -s "${PROKKA_FAA}" ] && [ -f "${DIAMOND_DB_UNIREF90}" ] && [ ! -f "${DIAMOND_OUT}" ]; then
        conda run -n "${ENV_DIAMOND}" diamond blastp \
            --db "${DIAMOND_DB_UNIREF90}" \
            --query "${PROKKA_FAA}" \
            --out "${DIAMOND_OUT}" \
            --outfmt 101 \
            --threads "${THREADS}" \
            --sensitive \
            --max-target-seqs 1 \
            --evalue 1e-5
        log_ok "[${SAMPLE_NAME}] Step 6: DIAMOND 완료"
    else
        log_info "[${SAMPLE_NAME}] Step 6: 결과 존재 또는 입력 없음, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 9. GTDB-Tk (CheckM 연동)
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 9: GTDB-Tk MAG Taxonomy (env: ${ENV_GTDBTK})"

    local GTDBTK_BINS_DIR="${BINS_DIR}/gtdbtk_input"
    local GTDBTK_OUT="${MAG_TAXONOMY_DIR}/gtdbtk_output_filtered"

    if [ "${RUN_GTDBTK}" = true ] && bins_exist "${BINS_FOR_CHECKM:-${BINS_DIR}}"; then
        mkdir -p "${GTDBTK_BINS_DIR}"

        if [ -f "${CHECKM_QA}" ]; then
            log_info "[${SAMPLE_NAME}] CheckM 필터 적용 (Completeness≥${CHECKM_MIN_COMPLETENESS}%, Contamination≤${CHECKM_MAX_CONTAMINATION}%)"
            while IFS= read -r binname; do
                local BIN_FA="${BINS_FOR_CHECKM:-${BINS_DIR}}/${binname}.fa"
                if [ -f "${BIN_FA}" ]; then
                    ln -sf "${BIN_FA}" "${GTDBTK_BINS_DIR}/"
                else
                    log_warn "[${SAMPLE_NAME}] bin 파일 없음: ${BIN_FA}"
                fi
            done < <(grep -v '^[-[:space:]]*$' "${CHECKM_QA}" \
                     | grep -v '^-' \
                     | awk -v comp="${CHECKM_MIN_COMPLETENESS}" \
                           -v cont="${CHECKM_MAX_CONTAMINATION}" \
                       'NR>1 && $7+0 >= comp && $8+0 <= cont {print $1}')
        else
            log_warn "[${SAMPLE_NAME}] CheckM QA 없음 → 전체 bin으로 GTDB-Tk 실행"
            GTDBTK_BINS_DIR="${BINS_FOR_CHECKM:-${BINS_DIR}}"
        fi

        if [ ! -d "${GTDBTK_OUT}/classify" ] || \
           [ -z "$(ls -A "${GTDBTK_OUT}/classify" 2>/dev/null)" ]; then
            conda run -n "${ENV_GTDBTK}" gtdbtk classify_wf \
                --genome_dir "${GTDBTK_BINS_DIR}" \
                --out_dir "${GTDBTK_OUT}" \
                --cpus "${THREADS}" \
                --extension fa \
                --skip_ani_screen \
                || log_warn "[${SAMPLE_NAME}] GTDB-Tk 실패, 다음 단계 계속..."
            log_ok "[${SAMPLE_NAME}] Step 9: GTDB-Tk 완료"
        else
            log_info "[${SAMPLE_NAME}] Step 9: 결과 존재, 건너뜀"
        fi
    else
        log_info "[${SAMPLE_NAME}] Step 9: bin 없음 또는 RUN_GTDBTK=false, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 10. eggNOG-mapper (diamond 경로 보정 적용)
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 10: eggNOG-mapper (env: ${ENV_EGGNOG})"

    local EGGNOG_PREFIX="${EGGNOG_DIR}/${SAMPLE_NAME}.eggnog"
    local EGGNOG_ANNOT="${EGGNOG_PREFIX}.emapper.annotations"
    local EGGNOG_TMP="${EGGNOG_DIR}/tmp_${SAMPLE_NAME}"

    if [ "${RUN_EGGNOG_MAPPER}" = true ] && [ -s "${PROKKA_FAA}" ] && [ ! -f "${EGGNOG_ANNOT}" ]; then
        mkdir -p "${EGGNOG_TMP}"
        # eggNOG 환경 자체에 diamond가 없을 수 있으므로, DIAMOND 환경의 bin을
        # PATH 앞쪽에 추가해 emapper.py가 diamond를 찾도록 보정한다.
        if ! conda run -n "${ENV_EGGNOG}" bash -c \
            "export PATH=\"${DIAMOND_BIN_DIR}:\$PATH\"; emapper.py \
                -i '${PROKKA_FAA}' \
                --output '${EGGNOG_PREFIX}' \
                -m diamond \
                --data_dir '${EGGNOG_DATA_DIR}' \
                --cpu '${THREADS}' \
                --temp_dir '${EGGNOG_TMP}' \
                --override \
                --target_orthologs all \
                --tax_scope auto \
                --annotate_hits_table yes \
                --go_evidence all"; then
            log_warn "[${SAMPLE_NAME}] Step 10: eggNOG-mapper 실패 (diamond 경로 문제일 수 있음)"
        else
            log_ok "[${SAMPLE_NAME}] Step 10: eggNOG-mapper 완료"
        fi
        rm -rf "${EGGNOG_TMP}"
    else
        log_info "[${SAMPLE_NAME}] Step 10: 결과 존재 또는 조건 미충족, 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 11. HUMAnN3
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 11: HUMAnN3 (env: ${ENV_HUMANN})"

    local HUMANN_OUT="${HUMANN_DIR}/humann_output"
    local HUMANN_GENEFAMILIES="${HUMANN_OUT}/${SAMPLE_NAME}_genefamilies.tsv"
    local COMBINED_FQ="${HUMANN_DIR}/${SAMPLE_NAME}.combined_for_humann.fastq.gz"

    if [ "${RUN_HUMANN}" = true ] && [ -s "${CURRENT_R1}" ]; then
        mkdir -p "${HUMANN_OUT}" "${HUMANN_TMP_DIR}"

        if [ ! -f "${COMBINED_FQ}" ]; then
            log_info "[${SAMPLE_NAME}] HUMAnN3 입력용 R1+R2 합치는 중..."
            local TEMP_R1="${HUMANN_DIR}/${SAMPLE_NAME}.temp.R1.fq"
            local TEMP_R2="${HUMANN_DIR}/${SAMPLE_NAME}.temp.R2.fq"
            conda run -n "${ENV_QC}" bash -c "pigz -dc \"${CURRENT_R1}\" > \"${TEMP_R1}\""
            conda run -n "${ENV_QC}" bash -c "pigz -dc \"${CURRENT_R2}\" > \"${TEMP_R2}\""
            conda run -n "${ENV_QC}" bash -c \
                "cat \"${TEMP_R1}\" \"${TEMP_R2}\" | pigz -p ${THREADS} > \"${COMBINED_FQ}\""
            rm -f "${TEMP_R1}" "${TEMP_R2}"
        fi

        if [ ! -f "${HUMANN_GENEFAMILIES}" ]; then
            conda run -n "${ENV_HUMANN}" bash -c "
                export TMPDIR='${HUMANN_TMP_DIR}' && \
                humann \
                    --input '${COMBINED_FQ}' \
                    --output '${HUMANN_OUT}' \
                    --threads '${THREADS}' \
                    --nucleotide-database '${HUMANN_NUCLEOTIDE_DB_PATH}' \
                    --protein-database '${HUMANN_PROTEIN_DB_PATH}' \
                    --output-basename '${SAMPLE_NAME}' \
                    --remove-temp-output"
            log_ok "[${SAMPLE_NAME}] Step 11: HUMAnN3 완료"
        else
            log_info "[${SAMPLE_NAME}] Step 11: 결과 존재, 건너뜀"
        fi
    else
        log_info "[${SAMPLE_NAME}] Step 11: 조건 미충족 (RUN_HUMANN=false 또는 read 없음), 건너뜀"
    fi

    # ------------------------------------------------------------------
    # Step 12. AMRFinderPlus
    # ------------------------------------------------------------------
    log_step "[${SAMPLE_NAME}] Step 12: AMRFinderPlus (env: ${ENV_AMRFINDER})"

    local AMRFINDER_OUT="${AMRFINDER_DIR}/${SAMPLE_NAME}.amrfinder.protein_only.tsv"

    if [ "${RUN_AMRFINDERPLUS}" = true ] && [ -s "${PROKKA_FAA}" ]; then
        if [ ! -f "${AMRFINDER_OUT}" ]; then
            conda run -n "${ENV_AMRFINDER}" amrfinder \
                -p "${PROKKA_FAA}" \
                -o "${AMRFINDER_OUT}" \
                --threads "${THREADS}" \
                --plus
            log_ok "[${SAMPLE_NAME}] Step 12: AMRFinderPlus 완료"
        else
            log_info "[${SAMPLE_NAME}] Step 12: 결과 존재, 건너뜀"
        fi
    else
        log_warn "[${SAMPLE_NAME}] Step 12: Prokka FAA 없음, AMRFinderPlus 건너뜀"
    fi

    append_summary "${SAMPLE_NAME}" \
        "${TOTAL_READS:-N/A}" \
        "${TRIMMED_READS:-N/A}" \
        "${N50:-N/A}" \
        "${CONTIG_COUNT:-N/A}" \
        "${HQ_MAG_COUNT:-N/A}"

    log_ok "[${SAMPLE_NAME}] 모든 단계 완료"
}

# ==============================================================================
# 6. 메인 실행
# ==============================================================================

main() {
    echo -e "${MAGENTA}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║   Shotgun Metagenome Pipeline — meta_hong 전용판             ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${RESET}"
    log_info "BASE_DIR    = ${BASE_DIR}"
    log_info "FASTQ_DIR   = ${FASTQ_DIR} (읽기 전용 참조, 복사하지 않음)"
    log_info "RESULTS_DIR = ${RESULTS_DIR}"

    mkdir -p "${RESULTS_DIR}"

    validate_environment
    warn_disk_estimate
    init_summary
    build_host_indices

    local FASTQ_LIST
    mapfile -t FASTQ_LIST < <(find "${FASTQ_DIR}" -name "*_1.fastq.gz" | sort)

    if [ "${#FASTQ_LIST[@]}" -eq 0 ]; then
        log_error "FASTQ 파일을 찾을 수 없습니다: ${FASTQ_DIR}/**/*_1.fastq.gz"
        exit 1
    fi

    if [ -n "${SAMPLE_INCLUDE_REGEX}" ]; then
        log_warn "SAMPLE_INCLUDE_REGEX='${SAMPLE_INCLUDE_REGEX}' 적용됨 — 일부 샘플만 실행합니다 (파일럿 모드)"
        local FILTERED=()
        for F in "${FASTQ_LIST[@]}"; do
            local NAME
            NAME=$(basename "${F}" _1.fastq.gz)
            if [[ "${NAME}" =~ ${SAMPLE_INCLUDE_REGEX} ]]; then
                FILTERED+=("${F}")
            fi
        done
        FASTQ_LIST=("${FILTERED[@]}")
    fi

    if [ "${#FASTQ_LIST[@]}" -eq 0 ]; then
        log_error "필터 적용 후 남은 샘플이 없습니다 (SAMPLE_INCLUDE_REGEX 확인)"
        exit 1
    fi

    log_info "처리할 샘플: ${#FASTQ_LIST[@]}개"
    for F in "${FASTQ_LIST[@]}"; do
        log_info "  - $(basename "${F}" _1.fastq.gz)"
    done

    for R1 in "${FASTQ_LIST[@]}"; do
        process_sample "${R1}"
    done

    # ------------------------------------------------------------------
    # Final. MultiQC
    # ------------------------------------------------------------------
    log_step "Final: MultiQC 전체 리포트 생성 (env: ${ENV_QC})"

    local MULTIQC_DIR="${RESULTS_DIR}/multiqc_report_all"
    mkdir -p "${MULTIQC_DIR}"

    log_info "MultiQC 스캔 범위 제한 적용 중 (대용량 데이터 파일 제외)..."
    conda run --no-capture-output -n "${ENV_QC}" multiqc \
        "${RESULTS_DIR}" \
        -o "${MULTIQC_DIR}" \
        -f \
        --profile-runtime \
        --ignore-symlinks \
        --ignore "*.fastq.gz" \
        --ignore "*.fa" \
        --ignore "*.fasta" \
        --ignore "*.faa" \
        --ignore "*.ffn" \
        --ignore "*.fna" \
        --ignore "*.bam" \
        --ignore "*.bai" \
        --ignore "*.sam" \
        --ignore "*.paf" \
        --ignore "*kraken2.out" \
        --ignore "*/checkm_output/bins/*" \
        --ignore "*/checkm_output/storage/*" \
        --ignore "*/gtdbtk_output_filtered/align/*" \
        --ignore "*/gtdbtk_output_filtered/identify/*" \
        --ignore "*humann_output*" \
        --ignore "*.combined_for_humann.fastq.gz" \
        || log_warn "MultiQC 실행 중 일부 경고 발생 (리포트는 생성되었을 수 있음)"
    log_ok "MultiQC 리포트 → ${MULTIQC_DIR}"

    echo ""
    echo -e "${GREEN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║           파이프라인 완료!                                    ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${RESET}"
    echo "  결과 디렉토리 : ${RESULTS_DIR}"
    echo "  Summary 파일  : ${SUMMARY_FILE}"
    echo "  MultiQC 리포트: ${MULTIQC_DIR}/multiqc_report.html"
}

main "$@"
