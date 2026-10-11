#!/usr/bin/env bash
# =============================================================================
#  run.sh — ONT Phage Assembly Pipeline Launcher
#  Edit the CONFIG section below, then:  bash run.sh
# =============================================================================

set -euo pipefail

# =============================================================================
#  ▼▼▼  EDIT THIS SECTION  ▼▼▼
# =============================================================================

# ── Input / Output ────────────────────────────────────────────────────────────
INPUT="/path/to/your/reads/*.fastq.gz"  # ← EDIT THIS: full path to your ONT reads
                                       # Single file : "/data/phage1.fastq.gz"
                                       # Multiple    : "/data/*.fastq.gz"
                                       # Do NOT point at a folder — use a glob
OUTDIR="results"                       # Where all outputs will be written

# ── Phage genome ──────────────────────────────────────────────────────────────
GENOME_SIZE="50k"               # Expected phage genome size
                                #   Small phages        → 5k – 50k
                                #   Sipho / Myo phage   → 40k – 150k
                                #   Jumbo phage         → 150k – 500k

# ── Read quality ──────────────────────────────────────────────────────────────
READ_TYPE="nano-raw"            # nano-raw   R9 / standard basecalling
                                # nano-hq    R10 / SUP-called (Q20+)
                                # nano-corr  error-corrected reads
                                # Note: no leading "--" — Nextflow CLI cannot
                                # parse values that start with "--"

# ── Resources ─────────────────────────────────────────────────────────────────
THREADS=16

# ── Conda environments ────────────────────────────────────────────────────────
# Nextflow needs FULL PATHS to activate existing envs (not just names).
# Set CONDA_BASE (find it with: conda info --base), then set the env NAME
# for each tool — paths are constructed automatically below.
CONDA_BASE="$(conda info --base 2>/dev/null || echo "")"   # auto-detected; override if needed

PORECHOP_ENV="base"             # porechop_abi v0.5.0 is in the base env
FILTLONG_ENV="filtlong"         # conda env containing filtlong
FLYE_ENV="flye"                 # conda env containing flye
MEDAKA_ENV="medaka"             # conda env containing medaka
PHABOX2_ENV="phabox2"           # conda env containing phabox2
PHAROKKA_ENV="pharokka"         # conda env containing pharokka
DNAAPLER_ENV="pharokka"         # dnaapler is bundled with pharokka env
STATS_ENV="assembly-stats"      # optional — falls back to pure awk if missing

PHABOX2_DB="/path/to/phabox_db_v2"            # ← EDIT THIS
PHAROKKA_DB="/path/to/pharokka_databases"     # ← EDIT THIS

MIN_LENGTH=1000                 # Filtlong: discard reads shorter than this (bp)
                                # Matches Flye's hard minimum overlap floor (≥1000)
TARGET_BASES=""                 # Filtlong: subsample to N bases, keeping longest reads
                                # "" = disabled (use all reads passing MIN_LENGTH)
                                # Example: 11000000 = ~100x of a 110 kb phage
                                # Set this if coverage is extremely high (e.g. host contamination)

# ── Host branch (optional) ────────────────────────────────────────────────────
ASSEMBLE_HOST=false             # true → assemble + type host from unmapped reads
POOL_HOST=false                 # true → also assemble pooled host reads (same host only)
HOST_GENOME_SIZE="5.5m"         # K. pneumoniae ≈ 5.5 Mb
MAP_ENV="snippy"                # conda env with minimap2 + samtools
KAPTIVE_ENV="kaptive"           # conda env with kaptive v3
MLST_ENV="mlst"                 # conda env with mlst

# Profile: conda | docker | standard
PROFILE="conda"

# ── Resume a previous run? ────────────────────────────────────────────────────
RESUME=false                    # true → restart from last checkpoint

# =============================================================================
#  ▲▲▲  END OF CONFIG  ▲▲▲
# =============================================================================

# ── Resolve env names to full paths ──────────────────────────────────────────
# Nextflow requires full paths to activate existing envs.
# 'base' is special: its path IS CONDA_BASE (not CONDA_BASE/envs/base).
resolve_conda_path() {
    local env_name="$1"
    if [[ "${env_name}" == "base" ]]; then
        echo "${CONDA_BASE}"
    else
        echo "${CONDA_BASE}/envs/${env_name}"
    fi
}

PORECHOP_CONDA=$(resolve_conda_path "${PORECHOP_ENV}")
FILTLONG_CONDA=$(resolve_conda_path "${FILTLONG_ENV}")
FLYE_CONDA=$(resolve_conda_path "${FLYE_ENV}")
MEDAKA_CONDA=$(resolve_conda_path "${MEDAKA_ENV}")
PHABOX2_CONDA=$(resolve_conda_path "${PHABOX2_ENV}")
PHAROKKA_CONDA=$(resolve_conda_path "${PHAROKKA_ENV}")
DNAAPLER_CONDA=$(resolve_conda_path "${DNAAPLER_ENV}")
MAP_CONDA=$(resolve_conda_path "${MAP_ENV}")
KAPTIVE_CONDA=$(resolve_conda_path "${KAPTIVE_ENV}")
MLST_CONDA=$(resolve_conda_path "${MLST_ENV}")
STATS_CONDA=$(resolve_conda_path "${STATS_ENV}")

# =============================================================================
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')]${RESET} $*"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔${RESET} $*"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠${RESET} $*"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘ ERROR:${RESET} $*" >&2; exit 1; }

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "\n${BOLD}════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}   ONT Phage Assembly Pipeline${RESET}"
echo -e "${BOLD}════════════════════════════════════════════════════${RESET}\n"

# =============================================================================
#  PRE-FLIGHT CHECKS
# =============================================================================

# ── 1. Nextflow ───────────────────────────────────────────────────────────────
log "Checking Nextflow..."
if ! command -v nextflow &>/dev/null; then
    warn "Nextflow not found. Installing to ${HOME}/.local/bin ..."
    mkdir -p "${HOME}/.local/bin"
    curl -fsSL https://get.nextflow.io | bash
    mv nextflow "${HOME}/.local/bin/nextflow"
    export PATH="${HOME}/.local/bin:${PATH}"
    ok "Nextflow installed: $(nextflow -version 2>&1 | head -1)"
else
    NF_VER=$(nextflow -version 2>&1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p' || echo "unknown")
    ok "Nextflow ${NF_VER} — $(command -v nextflow)"
    NF_MAJOR=$(echo "${NF_VER}" | cut -d. -f1)
    [[ "${NF_MAJOR}" -lt 23 ]] 2>/dev/null && \
        warn "Nextflow ≥ 23 is recommended. Run: nextflow self-update"
fi

# ── 2. Pipeline files ─────────────────────────────────────────────────────────
log "Checking pipeline files..."
for f in main.nf nextflow.config \
          modules/porechop_abi.nf \
          modules/filtlong.nf \
          modules/flye.nf \
          modules/medaka.nf \
          modules/dnaapler.nf \
          modules/phabox2.nf \
          modules/pharokka.nf \
          modules/assembly_stats.nf; do
    [[ -f "${PIPELINE_DIR}/${f}" ]] || err "Missing: ${PIPELINE_DIR}/${f}"
done
ok "All pipeline files present"

# ── 3. Input reads ────────────────────────────────────────────────────────────
log "Checking input reads..."

# Catch the default placeholder — user forgot to edit INPUT
if [[ "${INPUT}" == /path/to/your/* ]]; then
    err "INPUT is still the placeholder value. Edit run.sh and set INPUT to your actual reads file."
fi

INPUT_FILES=( ${INPUT} )

# No files matched
[[ ${#INPUT_FILES[@]} -eq 0 ]] || [[ ! -e "${INPUT_FILES[0]}" ]] && \
    err "No files found matching: ${INPUT}\n  Check the path is correct and the file exists."

# Reject directories — porechop_abi needs a file, not a folder
for f in "${INPUT_FILES[@]}"; do
    if [[ -d "${f}" ]]; then
        err "'${f}' is a directory, not a fastq file.\n" \
            "  If your reads are inside that folder, use a glob instead:\n" \
            "  INPUT=\"${f}/*.fastq.gz\""
    fi
    if [[ ! "${f}" =~ \.(fastq|fq|fastq\.gz|fq\.gz)$ ]]; then
        warn "'${f}' does not look like a fastq file (unexpected extension)."
    fi
done

ok "${#INPUT_FILES[@]} input file(s) found:"
for f in "${INPUT_FILES[@]}"; do
    SIZE=$(du -sh "${f}" 2>/dev/null | cut -f1 || echo "?")
    log "    ${f}  (${SIZE})"
done

# ── 4. Conda environments ─────────────────────────────────────────────────────
if [[ "${PROFILE}" == conda* ]]; then

    # Resolve conda command (prefer mamba)
    if   command -v mamba &>/dev/null; then CONDA_CMD="mamba"
    elif command -v conda &>/dev/null; then CONDA_CMD="conda"
    else err "Neither conda nor mamba found in PATH."
    fi
    ok "Using ${CONDA_CMD}"

    # ── Check a full-path env exists and contains a tool ─────────────────────
    # Direct binary check — no `conda run` needed, no env creation triggered.
    check_tool_in_env() {
        local env_path="$1"
        local tool="$2"

        if [[ ! -d "${env_path}" ]]; then
            err "Conda env path does not exist: ${env_path}" \
                $'\n  Check CONDA_BASE in run.sh (currently: '"${CONDA_BASE}"')' \
                $'\n  Run: conda info --base'
        fi

        if [[ ! -x "${env_path}/bin/${tool}" ]]; then
            local found_in
            found_in=$(conda env list 2>/dev/null | grep -v '^#' | awk '{print $NF}' | \
                while read -r p; do [[ -x "${p}/bin/${tool}" ]] && echo "${p}" && break; done || true)
            if [[ -n "${found_in}" ]]; then
                err "'${tool}' not in ${env_path}  (found in: ${found_in})" \
                    $'\n  Fix: update the matching ENV variable in run.sh'
            else
                err "'${tool}' not found in any conda env." \
                    $'\n  Install: conda install -c bioconda '"${tool}"
            fi
        fi
    }

    log "Checking conda environments..."

    log "  Checking porechop_abi in ${PORECHOP_CONDA} ..."
    check_tool_in_env "${PORECHOP_CONDA}" "porechop_abi"
    PABI_VER=$("${PORECHOP_CONDA}/bin/porechop_abi" --version 2>&1 | head -1 || echo "unknown")
    ok "  porechop_abi ${PABI_VER}  <- ${PORECHOP_CONDA}"

    log "  Checking filtlong in ${FILTLONG_CONDA} ..."
    check_tool_in_env "${FILTLONG_CONDA}" "filtlong"
    FL_VER=$("${FILTLONG_CONDA}/bin/filtlong" --version 2>&1 | head -1 || echo "unknown")
    ok "  filtlong ${FL_VER}  <- ${FILTLONG_CONDA}"

    log "  Checking flye in ${FLYE_CONDA} ..."
    check_tool_in_env "${FLYE_CONDA}" "flye"
    FLYE_VER=$("${FLYE_CONDA}/bin/flye" --version 2>&1 | head -1 || echo "unknown")
    ok "  flye ${FLYE_VER}  <- ${FLYE_CONDA}"

    log "  Checking medaka in ${MEDAKA_CONDA} ..."
    check_tool_in_env "${MEDAKA_CONDA}" "medaka_consensus"
    ok "  medaka  <- ${MEDAKA_CONDA}"

    log "  Checking phabox2 in ${PHABOX2_CONDA} ..."
    check_tool_in_env "${PHABOX2_CONDA}" "phabox2"
    ok "  phabox2  <- ${PHABOX2_CONDA}"
    [[ -d "${PHABOX2_DB}" ]] || err "PhaBOX2 database not found: ${PHABOX2_DB}"
    ok "  PhaBOX2 DB: ${PHABOX2_DB}"

    log "  Checking dnaapler in ${DNAAPLER_CONDA} ..."
    check_tool_in_env "${DNAAPLER_CONDA}" "dnaapler"
    ok "  dnaapler  <- ${DNAAPLER_CONDA}"

    log "  Checking pharokka in ${PHAROKKA_CONDA} ..."
    check_tool_in_env "${PHAROKKA_CONDA}" "pharokka.py"
    ok "  pharokka  <- ${PHAROKKA_CONDA}"
    [[ -d "${PHAROKKA_DB}" ]] || err "Pharokka database not found: ${PHAROKKA_DB}"
    ok "  Pharokka DB: ${PHAROKKA_DB}"

    if [[ "${ASSEMBLE_HOST}" == true ]]; then
        log "  Host branch: checking minimap2/samtools, kaptive, mlst ..."
        check_tool_in_env "${MAP_CONDA}" "minimap2"
        check_tool_in_env "${MAP_CONDA}" "samtools"
        check_tool_in_env "${KAPTIVE_CONDA}" "kaptive"
        check_tool_in_env "${MLST_CONDA}" "mlst"
        ok "  Host branch tools found"
    fi

    log "  Checking assembly-stats in ${STATS_CONDA} (optional)..."
    if [[ -d "${STATS_CONDA}" ]] && [[ -x "${STATS_CONDA}/bin/assembly-stats" ]]; then
        STATS_VER=$("${STATS_CONDA}/bin/assembly-stats" --version 2>&1 | head -1 || echo "unknown")
        ok "  assembly-stats ${STATS_VER}  <- ${STATS_CONDA}"
    else
        warn "  assembly-stats not found -- will use awk fallback"
        # Nextflow still needs a valid conda path for the process.
        # Fall back to base env — awk is always available there.
        STATS_CONDA="${PORECHOP_CONDA}"
    fi

elif [[ "${PROFILE}" == "standard" ]]; then
    log "Profile 'standard': checking tools in PATH..."
    command -v porechop_abi &>/dev/null || err "'porechop_abi' not found in PATH."
    command -v flye         &>/dev/null || err "'flye' not found in PATH."
    ok "porechop_abi and flye found in PATH"

elif [[ "${PROFILE}" == "docker" ]]; then
    command -v docker &>/dev/null || err "Docker not found."
    docker info &>/dev/null       || err "Docker daemon is not running."
    ok "Docker is running"

fi

# ── 5. Genome size format ─────────────────────────────────────────────────────
if [[ -n "${GENOME_SIZE}" ]]; then
    echo "${GENOME_SIZE}" | grep -qE '^[0-9]+(\.[0-9]+)?[kmgKMG]?$' || \
        err "Invalid GENOME_SIZE '${GENOME_SIZE}'. Examples: 50k, 150k, 2m"
    ok "Genome size: ${GENOME_SIZE}"
fi

# ── 6. Read type ──────────────────────────────────────────────────────────────
case "${READ_TYPE}" in
    nano-raw|nano-hq|nano-corr|pacbio-raw|pacbio-corr|pacbio-hifi)
        ok "Read type: ${READ_TYPE}" ;;
    *) err "Invalid READ_TYPE: ${READ_TYPE}" ;;
esac

# =============================================================================
#  BUILD NEXTFLOW COMMAND
# =============================================================================

CMD=(
    nextflow run "${PIPELINE_DIR}/main.nf"
    --input          "${INPUT}"
    --outdir         "${OUTDIR}"
    --read_type      "${READ_TYPE}"
    --threads        "${THREADS}"
    --no_split       true
    --porechop_conda "${PORECHOP_CONDA}"
    --filtlong_conda "${FILTLONG_CONDA}"
    --flye_conda     "${FLYE_CONDA}"
    --medaka_conda   "${MEDAKA_CONDA}"
    --phabox2_conda  "${PHABOX2_CONDA}"
    --phabox2_db     "${PHABOX2_DB}"
    --pharokka_conda "${PHAROKKA_CONDA}"
    --dnaapler_conda "${DNAAPLER_CONDA}"
    --pharokka_db    "${PHAROKKA_DB}"
    --stats_conda    "${STATS_CONDA}"
    --min_length     "${MIN_LENGTH}"
    -profile         "${PROFILE}"
    -with-report     "${OUTDIR}/pipeline_report.html"
    -with-timeline   "${OUTDIR}/pipeline_timeline.html"
    -with-dag        "${OUTDIR}/pipeline_dag.html"
    -with-trace      "${OUTDIR}/pipeline_trace.txt"
)

[[ -n "${GENOME_SIZE}" ]]          && CMD+=(--genome_size "${GENOME_SIZE}")
[[ -n "${TARGET_BASES}" ]]         && CMD+=(--target_bases "${TARGET_BASES}")
if [[ "${ASSEMBLE_HOST}" == true ]]; then
    CMD+=(--assemble_host --host_genome_size "${HOST_GENOME_SIZE}"
          --host_map_conda "${MAP_CONDA}" --kaptive_conda "${KAPTIVE_CONDA}" --mlst_conda "${MLST_CONDA}")
    [[ "${POOL_HOST}" == true ]] && CMD+=(--pool_host)
fi
[[ "${RESUME}" == true ]]          && CMD+=(-resume)
# CONDA_BASE already exported above

# =============================================================================
#  RUN SUMMARY
# =============================================================================

echo ""
echo -e "${BOLD}Run configuration:${RESET}"
echo -e "  Input           : ${CYAN}${INPUT}${RESET}"
echo -e "  Output dir      : ${CYAN}${OUTDIR}${RESET}"
echo -e "  Genome size     : ${CYAN}${GENOME_SIZE:-auto}${RESET}"
echo -e "  Read type       : ${CYAN}${READ_TYPE}${RESET}"
echo -e "  Threads         : ${CYAN}${THREADS}${RESET}"
echo -e "  Profile         : ${CYAN}${PROFILE}${RESET}"
echo -e "  Porechop conda  : ${CYAN}${PORECHOP_CONDA}${RESET}"
echo -e "  Filtlong conda  : ${CYAN}${FILTLONG_CONDA}${RESET}"
echo -e "  Min read length : ${CYAN}${MIN_LENGTH} bp${RESET}"
echo -e "  Flye conda      : ${CYAN}${FLYE_CONDA}${RESET}"
echo -e "  Medaka conda    : ${CYAN}${MEDAKA_CONDA}${RESET}"
echo -e "  PhaBOX2 conda   : ${CYAN}${PHABOX2_CONDA}${RESET}"
echo -e "  PhaBOX2 DB      : ${CYAN}${PHABOX2_DB}${RESET}"
echo -e "  Pharokka conda  : ${CYAN}${PHAROKKA_CONDA}${RESET}"
echo -e "  Pharokka DB     : ${CYAN}${PHAROKKA_DB}${RESET}"
echo -e "  Stats conda     : ${CYAN}${STATS_CONDA}${RESET}"
echo -e "  No split        : ${CYAN}true${RESET} (phage default)"
echo -e "  Resume          : ${CYAN}${RESUME}${RESET}"
echo ""
echo -e "${BOLD}Command:${RESET}"
echo -e "  ${CYAN}${CMD[*]}${RESET}\n"

# =============================================================================
#  LAUNCH
# =============================================================================

mkdir -p "${OUTDIR}"
START=$(date +%s)
log "Starting pipeline..."
echo -e "${BOLD}════════════════════════════════════════════════════${RESET}\n"

set +e
"${CMD[@]}"
EXIT_CODE=$?
set -e

echo -e "\n${BOLD}════════════════════════════════════════════════════${RESET}"
ELAPSED=$(( $(date +%s) - START ))
ELAPSED_FMT=$(printf '%dh %dm %ds' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))

if [[ ${EXIT_CODE} -eq 0 ]]; then
    ok "Pipeline finished in ${ELAPSED_FMT}"
    echo ""
    echo -e "${BOLD}Key outputs:${RESET}"
    echo -e "  Trimmed reads   : ${OUTDIR}/01_porechop_abi/"
    echo -e "  Assembly FASTA  : ${OUTDIR}/03_flye_assembly/<sample>/assembly.fasta"
    echo -e "  Assembly graph  : ${OUTDIR}/03_flye_assembly/<sample>/assembly_graph.gfa"
    echo -e "  Assembly info   : ${OUTDIR}/03_flye_assembly/<sample>/assembly_info.txt"
    echo -e "  Stats summary   : ${OUTDIR}/summary/assembly_summary.tsv"
    echo -e "  Run report      : ${OUTDIR}/pipeline_report.html"
    echo ""
    echo -e "${BOLD}Assembly topology (phage QC):${RESET}"
    ASSEMBLY_INFO=$(find "${OUTDIR}/03_flye_assembly" -name "assembly_info.txt" 2>/dev/null | head -1)
    if [[ -f "${ASSEMBLY_INFO}" ]]; then
        echo ""
        cat "${ASSEMBLY_INFO}"
        echo ""
        grep -q $'\tY\t' "${ASSEMBLY_INFO}" 2>/dev/null \
            && ok "Circular contig reported by Flye" \
            || warn "No circular contig reported. Many phages are linear — judge completeness by contig number, length and coverage vs. a reference."
    fi
else
    echo -e "  Set ${BOLD}RESUME=true${RESET} and re-run to restart from the last checkpoint."
    err "Pipeline failed (exit ${EXIT_CODE}). Check logs in ${OUTDIR}/."
fi
