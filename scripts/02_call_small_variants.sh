#!/usr/bin/env bash
# Public portfolio version: paths, scheduler names, sample IDs, and project-specific defaults are parameterized.
# ============================================================================
# HiFi 流程 Step2a: Small variant 检测 (DeepVariant, GPU)
#
# 功能：
#   1) DeepVariant (GPU, Parabricks): 每个样本从 BAM 跑到 GVCF
#   2) 从 GVCF 派生单样本 VCF (bcftools, 供下游注释/分析)
#
# 说明：
#   - 本脚本只做 DeepVariant，SV (Sniffles) 拆到 step2b 走 SLURM 队列
#   - 输出 GVCF (兼容 step3a GLnexus 多样本联合分型) + 派生 VCF (单样本下游用)
#
# 输入：
#   ${BASE_OUT}/aligned_bam/${sample}.minimap2.bam
#
# 输出：
#   ${BASE_OUT}/dv_gvcf/${sample}.dv.g.vcf.gz  + .tbi   (GVCF)
#   ${BASE_OUT}/dv_vcf/${sample}.dv.vcf.gz     + .tbi   (派生 VCF)
#   ${BASE_OUT}/log/step2a/${sample}.log + .timing.log
#
# 投递方式（GPU 节点，无 SLURM，nohup 串行）：
#   ssh <gpu-node>
#   cd /path/to/hifi-project/sh
#   tmux new -s hifi_step2a
#   BASE_OUT=/path/to/hifi-project \
#   SAMPLE_LIST=/path/to/hifi-project/samples.list \
#   nohup bash 02_call_small_variants.sh > step2a.main.log 2>&1 &
#
# 进度排查：
#   跑完的样本 : grep -l 'ALL DONE' ${BASE_OUT}/log/step2a/*.timing.log | wc -l
#   未跑完     : grep -L 'ALL DONE' ${BASE_OUT}/log/step2a/*.timing.log
#   各步耗时   : grep 'END' ${BASE_OUT}/log/step2a/*.timing.log
#   失败清单   : ${BASE_OUT}/log/step2a/step2a_failed.tsv
# ============================================================================
set -uo pipefail

###########################
# 环境
###########################
ENV_FILE="${ENV_FILE:-}"
[[ -n "$ENV_FILE" && -s "$ENV_FILE" ]] && source "$ENV_FILE"

###########################
# 配置（允许 env 覆盖）
###########################
BASE_OUT="${BASE_OUT:-${PWD}/work}"
ALIGNED_DIR="${ALIGNED_DIR:-${BASE_OUT}/aligned_bam}"
SAMPLE_LIST="${SAMPLE_LIST:-${BASE_OUT}/samples.list}"

REF_FASTA="${REF_FASTA:-${BASE_OUT}/reference/GRCh38.fasta}"
REF_ROOT="${REF_ROOT:-$(dirname "${REF_FASTA}")}"

DV_MODEL_TYPE="${DV_MODEL_TYPE:-pacbio}"

IMAGE="${IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1}"
DOCKER_GPUS="${DOCKER_GPUS:-all}"
DOCKER_USER="${DOCKER_USER:-$(id -u):$(id -g)}"

# bcftools + tabix 绝对路径（用于从 GVCF 派生 VCF）
BCFTOOLS="${BCFTOOLS:-bcftools}"
TABIX="${TABIX:-tabix}"

# 是否派生单样本 VCF (1=派生, 0=只留 GVCF)
DERIVE_VCF="${DERIVE_VCF:-1}"

###########################
# 输出目录（平铺）
###########################
DV_GVCF_DIR="${BASE_OUT}/dv_gvcf"
DV_VCF_DIR="${BASE_OUT}/dv_vcf"
LOG_DIR="${BASE_OUT}/log/step2a"

mkdir -p "$DV_GVCF_DIR" "$DV_VCF_DIR" "$LOG_DIR"

FAIL_TSV="${LOG_DIR}/step2a_failed.tsv"
[[ -f "$FAIL_TSV" ]] || echo -e "sample\tstep\treason\ttime" > "$FAIL_TSV"

###########################
# 工具函数
###########################
fmt_time() { local s=$1; printf '%dh%02dm%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }

SAMPLE_CUR=""
TIMING=""
STEP_START_TS=0

step_start() {
  STEP_START_TS=$(date +%s)
  echo "[$(date '+%F %T')] 开始 ${SAMPLE_CUR} $1..."
  printf '[%s] START  %s\n' "$(date '+%F %T')" "$1" >> "${TIMING}"
}
step_end() {
  local end_ts elapsed
  end_ts=$(date +%s)
  elapsed=$((end_ts - STEP_START_TS))
  echo "[$(date '+%F %T')] 完成 ${SAMPLE_CUR} $1，耗时 $(fmt_time ${elapsed})"
  printf '[%s] END    %-16s elapsed=%-7s (%s)\n' \
    "$(date '+%F %T')" "$1" "${elapsed}s" "$(fmt_time ${elapsed})" >> "${TIMING}"
}
log_fail() {
  echo -e "$1\t$2\t$3\t$(date '+%F %T')" >> "$FAIL_TSV"
  printf '[%s] FAILED %s reason=%s\n' "$(date '+%F %T')" "$2" "$3" >> "${TIMING}"
}

docker_gpus_arg() {
  local g="$1"
  if [[ -z "$g" || "$g" == "all" ]]; then echo "all"
  elif [[ "$g" =~ ^[0-9]+(,[0-9]+)*$ ]]; then echo "device=${g}"
  elif [[ "$g" =~ ^device=.*$ ]]; then echo "$g"
  else echo "$g"; fi
}
DOCKER_GPUS_FMT="$(docker_gpus_arg "$DOCKER_GPUS")"

###########################
# 检查
###########################
[[ -s "${SAMPLE_LIST}" ]] || { echo "[FATAL] 样本列表不存在：${SAMPLE_LIST}" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]]   || { echo "[FATAL] 参考缺失：${REF_FASTA}" >&2; exit 1; }
[[ -d "${ALIGNED_DIR}" ]] || { echo "[FATAL] BAM 目录不存在：${ALIGNED_DIR}" >&2; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "[FATAL] docker 不可用，请在 GPU 节点运行" >&2; exit 1; }

if [[ "${DERIVE_VCF}" == "1" ]]; then
  command -v "${BCFTOOLS}" >/dev/null 2>&1 || { echo "[FATAL] bcftools not found: ${BCFTOOLS}（或设 DERIVE_VCF=0）" >&2; exit 1; }
  command -v "${TABIX}" >/dev/null 2>&1 || { echo "[FATAL] tabix not found: ${TABIX}（或设 DERIVE_VCF=0）" >&2; exit 1; }
fi

REF_IN="/refdir/$(basename "$REF_FASTA")"

# 一次性读取样本表，避免 NFS 上长期打开文件导致 stale file handle
mapfile -t SAMPLES < <(sed 's/\r$//' "${SAMPLE_LIST}" | awk 'NF')
[[ ${#SAMPLES[@]} -gt 0 ]] || {
  echo "[FATAL] 样本列表为空：${SAMPLE_LIST}" >&2
  exit 1
}

###########################
# 头信息
###########################
echo "=========================================="
echo "HiFi Step2a: DeepVariant (GPU)"
echo "BASE_OUT       : ${BASE_OUT}"
echo "SAMPLE_LIST    : ${SAMPLE_LIST}"
echo "REF_FASTA      : ${REF_FASTA}"
echo "DV_GVCF_DIR    : ${DV_GVCF_DIR}"
echo "DV_VCF_DIR     : ${DV_VCF_DIR}"
echo "IMAGE          : ${IMAGE}"
echo "DOCKER_GPUS    : ${DOCKER_GPUS} (--gpus ${DOCKER_GPUS_FMT})"
echo "DV_MODEL_TYPE  : ${DV_MODEL_TYPE} (输出 GVCF)"
echo "DERIVE_VCF     : ${DERIVE_VCF}"
echo "=========================================="

###########################
# 主循环
###########################
for sample in "${SAMPLES[@]}"; do
  [[ -z "${sample}" ]] && continue

  SAMPLE_CUR="${sample}"
  TIMING="${LOG_DIR}/${sample}.timing.log"
  log_file="${LOG_DIR}/${sample}.log"

  echo
  echo ">>> 处理样本: ${sample}    开始: $(date)"

  in_bam="${ALIGNED_DIR}/${sample}.minimap2.bam"
  out_gvcf="${DV_GVCF_DIR}/${sample}.dv.g.vcf.gz"
  out_vcf="${DV_VCF_DIR}/${sample}.dv.vcf.gz"

  # ---- 断点续跑检查 ----
  skip_dv=0
  gvcf_ok=0
  vcf_ok=0

  # GVCF、索引存在，并且 tabix 可以正常读取
  if [[ -s "${out_gvcf}" && -s "${out_gvcf}.tbi" ]] \
     && "${TABIX}" -l "${out_gvcf}" >/dev/null 2>&1; then
    gvcf_ok=1
  fi

  # 派生 VCF、索引存在，并且 tabix 可以正常读取
  if [[ -s "${out_vcf}" && -s "${out_vcf}.tbi" ]] \
     && "${TABIX}" -l "${out_vcf}" >/dev/null 2>&1; then
    vcf_ok=1
  fi

  if [[ ${gvcf_ok} -eq 1 ]]; then
    if [[ "${DERIVE_VCF}" == "1" && ${vcf_ok} -ne 1 ]]; then
      echo ">>> GVCF 已存在且有效，跳过 DeepVariant，仅补派生 VCF"
      skip_dv=1
    else
      echo ">>> DeepVariant GVCF + VCF 已存在且有效，跳过整个样本"
      printf '[%s] SKIP   deepvariant (all outputs valid)\n' \
        "$(date '+%F %T')" >> "${TIMING}"
      printf '[%s] ALL DONE  %s\n' \
        "$(date '+%F %T')" "${sample}" >> "${TIMING}"
      continue
    fi
  fi

  : > "${TIMING}"

  # ---- 输入 BAM 检查 ----
  if [[ ! -f "${in_bam}" ]]; then
    echo "[WARN] 输入 BAM 不存在：${in_bam}" >&2
    log_fail "${sample}" "input" "bam_not_found"
    continue
  fi

  in_bam_container="/workdir/aligned_bam/${sample}.minimap2.bam"

  # ---- DeepVariant (GPU, 出 GVCF) ----
  if [[ ${skip_dv} -eq 0 ]]; then
    step_start deepvariant
    if ! docker run --gpus "${DOCKER_GPUS_FMT}" --user "${DOCKER_USER}" --rm \
          --workdir /workdir \
          -v "${REF_ROOT}":/refdir \
          -v "${BASE_OUT}":/workdir \
          "${IMAGE}" pbrun deepvariant \
            --ref "${REF_IN}" \
            --in-bam "${in_bam_container}" \
            --mode "${DV_MODEL_TYPE}" \
            --gvcf \
            --out-variants "/workdir/dv_gvcf/${sample}.dv.g.vcf.gz" \
          >> "${log_file}" 2>&1; then
      log_fail "${sample}" "deepvariant" "deepvariant_exit_nonzero"
      continue
    fi
    step_end deepvariant

    # 索引 GVCF
    if [[ -f "${out_gvcf}" && ! -f "${out_gvcf}.tbi" ]]; then
      step_start gvcf_index
      if ${TABIX} -p vcf -f "${out_gvcf}" >> "${log_file}" 2>&1; then
        step_end gvcf_index
      else
        log_fail "${sample}" "gvcf_index" "tabix_failed"
      fi
    fi
  fi

  # ---- 从 GVCF 派生单样本 VCF ----
  if [[ "${DERIVE_VCF}" == "1" && ! -f "${out_vcf}" ]]; then
    step_start derive_vcf
    # 过滤掉 non-variant reference block (<NON_REF> 和 .)，只留真变异
    if ${BCFTOOLS} view -e 'ALT="<NON_REF>" || ALT="."' \
         -Oz -o "${out_vcf}" "${out_gvcf}" >> "${log_file}" 2>&1 \
       && ${TABIX} -p vcf -f "${out_vcf}" >> "${log_file}" 2>&1; then
      step_end derive_vcf
    else
      log_fail "${sample}" "derive_vcf" "bcftools_or_tabix_failed"
    fi
  fi

  echo ">>> 完成: ${sample}"
  printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${sample}" >> "${TIMING}"
  echo ">>> 结束: $(date)"
  echo "---"
done

echo
echo "=========================================="
echo "Step2a 完成"
echo "  GVCF          : ${DV_GVCF_DIR}"
echo "  派生 VCF      : ${DV_VCF_DIR}"
echo "  每样本 log    : ${LOG_DIR}/<sample>.log"
echo "  失败清单      : ${FAIL_TSV}"
echo ""
echo "  下一步: 跑 step2b (Sniffles SV, SLURM 队列)"
echo "         多样本联合分型: step3a (GLnexus, 用 dv_gvcf/*.dv.g.vcf.gz)"
echo "=========================================="
