#!/usr/bin/env bash
# Public portfolio version: paths, scheduler names, sample IDs, and project-specific defaults are parameterized.
# ============================================================================
# HiFi 流程 Step3a: Small variant 多样本联合分型 (GLnexus)
#
# 功能：
#   把 step2 输出的多个 GVCF 联合分型为 cohort VCF。
#
# 输入：
#   ${BASE_OUT}/dv_gvcf/*.dv.g.vcf.gz   (step2 产物)
#
# 输出：
#   ${BASE_OUT}/cohort/${COHORT_NAME}.small.vcf.gz  + .tbi
#
# 关键差异（对比二代 GLnexus）：
#   - short-read DeepVariant 用 --config DeepVariantWGS
#   - HiFi   用 --config DeepVariant_unfiltered
#     (二代 WGS 配置内置的过滤参数针对 Illumina 错误模式，套到 HiFi 上过度过滤)
#
# 已知坑（practical note）：
#   GLnexus 输出 FILTER 列是 "." 而非 "PASS"，用 bcftools 时需 -f PASS,.
#
# 投递方式（GPU 节点，因为 CPU 节点没 docker）：
#   ssh <gpu-node>
#   cd /path/to/hifi-project/sh
#   tmux new -s hifi_step3a
#   BASE_OUT=/path/to/hifi-project \
#   COHORT_NAME=cohort1 \
#   nohup bash 04_joint_small_variants.sh > step3a.main.log 2>&1 &
#
# env 覆盖：
#   BASE_OUT       项目根目录（必须）
#   COHORT_NAME    cohort 名字（默认 cohort）
#   GLNEXUS_IMAGE  GLnexus 镜像（默认 ghcr.io/dnanexus-rnd/glnexus:v1.4.1）
#   GLNEXUS_CONFIG GLnexus 配置（默认 DeepVariant_unfiltered，HiFi 专用）
#   GLNEXUS_MEM    GLnexus 内存上限 GB（默认 64）
#   GLNEXUS_THREADS 线程数（默认 16）
# ============================================================================
set -uo pipefail

###########################
# 配置
###########################
BASE_OUT="${BASE_OUT:-${PWD}/work}"
COHORT_NAME="${COHORT_NAME:-cohort}"

DV_GVCF_DIR="${DV_GVCF_DIR:-${BASE_OUT}/dv_gvcf}"
COHORT_DIR="${COHORT_DIR:-${BASE_OUT}/cohort}"
LOG_DIR="${LOG_DIR:-${BASE_OUT}/log/step3}"

# GLnexus 配置（HiFi 专用，不用 DeepVariantWGS）
GLNEXUS_IMAGE="${GLNEXUS_IMAGE:-ghcr.io/dnanexus-rnd/glnexus:v1.4.1}"
GLNEXUS_CONFIG="${GLNEXUS_CONFIG:-DeepVariant_unfiltered}"
GLNEXUS_MEM="${GLNEXUS_MEM:-64}"
GLNEXUS_THREADS="${GLNEXUS_THREADS:-16}"

# bcftools + tabix 用绝对路径（不同环境）
BCFTOOLS="${BCFTOOLS:-bcftools}"
TABIX="${TABIX:-tabix}"

DOCKER_USER="${DOCKER_USER:-$(id -u):$(id -g)}"

###########################
# 输出目录
###########################
mkdir -p "$COHORT_DIR" "$LOG_DIR"

FAIL_TSV="${LOG_DIR}/step3_failed.tsv"
[[ -f "$FAIL_TSV" ]] || echo -e "cohort\tstep\treason\ttime" > "$FAIL_TSV"

TIMING="${LOG_DIR}/${COHORT_NAME}.timing.log"
LOG_FILE="${LOG_DIR}/${COHORT_NAME}.glnexus.log"

###########################
# 工具函数
###########################
fmt_time() { local s=$1; printf '%dh%02dm%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }

STEP_START_TS=0
step_start() {
  STEP_START_TS=$(date +%s)
  echo "[$(date '+%F %T')] 开始 $1..."
  printf '[%s] START  %s\n' "$(date '+%F %T')" "$1" >> "${TIMING}"
}
step_end() {
  local end_ts elapsed
  end_ts=$(date +%s)
  elapsed=$((end_ts - STEP_START_TS))
  echo "[$(date '+%F %T')] 完成 $1，耗时 $(fmt_time ${elapsed})"
  printf '[%s] END    %-20s elapsed=%-7s (%s)\n' \
    "$(date '+%F %T')" "$1" "${elapsed}s" "$(fmt_time ${elapsed})" >> "${TIMING}"
}
log_fail() {
  echo -e "${COHORT_NAME}\t$1\t$2\t$(date '+%F %T')" >> "$FAIL_TSV"
  printf '[%s] FAILED %s reason=%s\n' "$(date '+%F %T')" "$1" "$2" >> "${TIMING}"
}

###########################
# 检查
###########################
[[ -d "${DV_GVCF_DIR}" ]] || { echo "[FATAL] GVCF 目录不存在：${DV_GVCF_DIR}" >&2; exit 1; }

# 收集所有 GVCF
mapfile -t gvcf_files < <(ls "${DV_GVCF_DIR}"/*.dv.g.vcf.gz 2>/dev/null | sort)
n_samples=${#gvcf_files[@]}
[[ ${n_samples} -ge 1 ]] || { echo "[FATAL] 没找到 GVCF 文件：${DV_GVCF_DIR}/*.dv.g.vcf.gz" >&2; exit 1; }

if [[ ${n_samples} -lt 2 ]]; then
  echo "[WARN] 只有 ${n_samples} 个样本，联合分型意义不大，但仍继续。" >&2
fi

# 检查 bcftools/tabix
command -v "${BCFTOOLS}" >/dev/null 2>&1 || { echo "[FATAL] bcftools not found: ${BCFTOOLS}" >&2; exit 1; }
command -v "${TABIX}" >/dev/null 2>&1 || { echo "[FATAL] tabix not found: ${TABIX}" >&2; exit 1; }

# 检查 docker
command -v docker >/dev/null 2>&1 || {
  echo "[FATAL] docker 不可用。GLnexus 需要 docker 镜像，请在 GPU 节点上运行。" >&2
  exit 1
}

# 输出文件
OUT_BCF="${COHORT_DIR}/${COHORT_NAME}.small.bcf"
OUT_VCF="${COHORT_DIR}/${COHORT_NAME}.small.vcf.gz"

###########################
# 断点续跑检查
###########################
if [[ -f "${OUT_VCF}" && -f "${OUT_VCF}.tbi" ]]; then
  if gunzip -c "${OUT_VCF}" | head -10 | grep -q "##fileformat=VCF"; then
    echo ">>> ${OUT_VCF} 已存在且完整，跳过"
    printf '[%s] SKIP   glnexus (cohort VCF exists)\n' "$(date '+%F %T')" >> "${TIMING}"
    printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${COHORT_NAME}" >> "${TIMING}"
    exit 0
  fi
fi

# 本次运行 reset timing
: > "${TIMING}"

###########################
# 头信息
###########################
echo "=========================================="
echo "HiFi Step3a: Small variant 联合分型 (GLnexus)"
echo "COHORT_NAME       : ${COHORT_NAME}"
echo "DV_GVCF_DIR       : ${DV_GVCF_DIR}"
echo "样本数            : ${n_samples}"
echo "COHORT_DIR        : ${COHORT_DIR}"
echo "GLNEXUS_IMAGE     : ${GLNEXUS_IMAGE}"
echo "GLNEXUS_CONFIG    : ${GLNEXUS_CONFIG}"
echo "GLNEXUS_MEM       : ${GLNEXUS_MEM} GB"
echo "GLNEXUS_THREADS   : ${GLNEXUS_THREADS}"
echo "OUT_VCF           : ${OUT_VCF}"
echo "=========================================="
echo ""
echo "样本列表:"
for f in "${gvcf_files[@]}"; do
  echo "  $(basename "$f")"
done
echo ""

###########################
# 清理 GLnexus 工作目录（每次跑必须是空的，GLnexus 强制要求）
###########################
GLNEXUS_WORKDIR="${COHORT_DIR}/.glnexus_workdir_${COHORT_NAME}"
rm -rf "${GLNEXUS_WORKDIR}"

###########################
# Step 1: GLnexus 联合分型
###########################
step_start glnexus

# GVCF 在容器内的路径列表
container_gvcfs=()
for f in "${gvcf_files[@]}"; do
  container_gvcfs+=("/workdir/dv_gvcf/$(basename "$f")")
done

# 输出 BCF（GLnexus 原生输出格式，后面转 VCF）
container_bcf="/workdir/cohort/${COHORT_NAME}.small.bcf"
container_workdir="/workdir/cohort/.glnexus_workdir_${COHORT_NAME}"

if ! docker run --rm --user "${DOCKER_USER}" \
      --workdir /workdir \
      -v "${BASE_OUT}":/workdir \
      "${GLNEXUS_IMAGE}" \
      /usr/local/bin/glnexus_cli \
        --config "${GLNEXUS_CONFIG}" \
        --dir "${container_workdir}" \
        --mem-gbytes "${GLNEXUS_MEM}" \
        --threads "${GLNEXUS_THREADS}" \
        "${container_gvcfs[@]}" \
      > "${OUT_BCF}" 2>> "${LOG_FILE}"; then
  log_fail "glnexus" "glnexus_exit_nonzero"
  echo "[FATAL] GLnexus 失败，查看日志：${LOG_FILE}" >&2
  exit 1
fi

step_end glnexus

# GLnexus 工作目录 GLnexus 自己会保留（内部数据库），跑完可以删除
rm -rf "${GLNEXUS_WORKDIR}"

###########################
# Step 2: BCF → VCF.gz + tabix
###########################
step_start bcf_to_vcf

if ! ${BCFTOOLS} view "${OUT_BCF}" -Oz -o "${OUT_VCF}" 2>> "${LOG_FILE}"; then
  log_fail "bcf_to_vcf" "bcftools_view_failed"
  exit 1
fi

if ! ${TABIX} -p vcf -f "${OUT_VCF}" 2>> "${LOG_FILE}"; then
  log_fail "bcf_to_vcf" "tabix_failed"
  exit 1
fi

step_end bcf_to_vcf

# 清理中间 BCF
rm -f "${OUT_BCF}"

###########################
# Step 3: 基础统计（不做过滤，只输出信息）
###########################
step_start stats

echo "" >> "${LOG_FILE}"
echo "===== Cohort VCF 统计 =====" >> "${LOG_FILE}"
echo "样本数: $(${BCFTOOLS} query -l "${OUT_VCF}" | wc -l)" >> "${LOG_FILE}"
echo "总变异数: $(${BCFTOOLS} view -H "${OUT_VCF}" | wc -l)" >> "${LOG_FILE}"
echo "SNP 数: $(${BCFTOOLS} view -H -v snps "${OUT_VCF}" | wc -l)" >> "${LOG_FILE}"
echo "INDEL 数: $(${BCFTOOLS} view -H -v indels "${OUT_VCF}" | wc -l)" >> "${LOG_FILE}"
echo "" >> "${LOG_FILE}"
echo "样本列表:" >> "${LOG_FILE}"
${BCFTOOLS} query -l "${OUT_VCF}" >> "${LOG_FILE}"

step_end stats

###########################
# 完成
###########################
printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${COHORT_NAME}" >> "${TIMING}"

echo
echo "=========================================="
echo "Step3a 完成"
echo "  Cohort VCF     : ${OUT_VCF}"
echo "  日志           : ${LOG_FILE}"
echo "  timing         : ${TIMING}"
echo ""
echo "  下游处理提醒（避坑）:"
echo "    GLnexus 输出 FILTER 列是 '.' 不是 'PASS'"
echo "    过滤时用：bcftools view -f PASS,. ..."
echo "    (直接 -f PASS 会过滤掉所有变异)"
echo ""
echo "  样本数: $(${BCFTOOLS} query -l "${OUT_VCF}" | wc -l)"
echo "  变异数: $(${BCFTOOLS} view -H "${OUT_VCF}" | wc -l)"
echo "=========================================="
