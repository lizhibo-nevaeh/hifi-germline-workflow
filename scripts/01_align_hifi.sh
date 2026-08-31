#!/usr/bin/env bash
# Public portfolio version: paths, scheduler names, sample IDs, and project-specific defaults are parameterized.
# ============================================================================
# HiFi GPU 流程 Step1: Parabricks minimap2 (HiFi 比对)
#
# 蓝本：somatic_step1_gpu.sh (体细胞 GPU step1)
# 改动：
#   1) 替换比对工具：BWA fq2bam → minimap2 (HiFi 长读)
#   2) 输入是 unaligned HiFi BAM（不是 paired FASTQ）
#   3) 一个样本可以由多个 cell BAM 组成，先 samtools cat 合并再比对
#   4) 使用 --pbmm2 选项匹配 PacBio 官方 pbmm2 输出（HiFi 推荐）
#   5) 不做 BQSR、不做 markdup（HiFi 数据不需要：无 PCR dup、QV 已校准）
#   6) 输出结构对齐：01_minimap2/log/{sample}.log + {sample}.timing.log
#                  aligned_bam/{sample}.minimap2.bam (+ .bai)
#
# 投递方式（GPU 节点无 SLURM，nohup 后台串行）：
#   ssh <gpu-node>                 # 或对应 GPU 节点
#   cd /path/to/hifi-project
#   tmux new -s hifi_step1
#   nohup bash 01_align_hifi.sh > hifi_step1.main.log 2>&1 &
#
# env 覆盖示例：
#   BASE_OUT=/path/to/hifi-project \
#   SAMPLE_LIST=/path/to/hifi-project/samples.list \
#   DOCKER_GPUS=0 \
#   nohup bash 01_align_hifi.sh > step1.main.log 2>&1 &
#
# 输入格式：
#   ${SAMPLE_LIST} 每行一个样本名，例如:
#     SAMPLE01
#     SAMPLE02
#     SAMPLE03
#
#   每个样本的 cell BAM 支持两种放置（按优先级）：
#     A) 子目录方式（推荐，多样本必须用这个）:
#          ${DATA_DIR}/${sample}/*.hifi_reads.bam
#          例: data/SAMPLE02/SAMPLE02_PacBio-Revio_xxx.hifi_reads.bam
#     B) 平铺+样本名前缀方式（单样本/历史布局兼容）:
#          ${DATA_DIR}/${sample}*.hifi_reads.bam
#          例: data/SAMPLE01_PacBio-Revio_xxx.hifi_reads.bam
#
#   【安全保证】只匹配属于该 ${sample} 的 BAM：
#     - 子目录方式天然隔离
#     - 平铺方式要求文件名以 ${sample} 开头，不会误合并别的样本
#     - 绝不会无脑扫描整个 data/ 目录（旧版本的危险回退已移除）
#
#   多样本务必用子目录方式（A），避免平铺时样本名前缀撞车
#   （例如样本 HG00 会误匹配 SAMPLE01/SAMPLE02…，子目录方式无此风险）
#
# 进度排查：
#   跑完的样本 : grep -l 'ALL DONE' ${BASE_OUT}/01_minimap2/log/*.timing.log | wc -l
#   未跑完     : grep -L 'ALL DONE' ${BASE_OUT}/01_minimap2/log/*.timing.log
#   各步耗时   : grep 'END' ${BASE_OUT}/01_minimap2/log/*.timing.log
#   失败清单   : ${BASE_OUT}/01_minimap2/log/step1_failed.tsv
# ============================================================================
set -uo pipefail   # 不用 -e，单样本失败不中断

###########################
# 环境
###########################
ENV_FILE="${ENV_FILE:-}"
[[ -n "$ENV_FILE" && -s "$ENV_FILE" ]] && source "$ENV_FILE"

# 把 hifi env 的 PATH 也加进来（用到 samtools/sniffles 等 CPU 工具时备用）
HIFI_ENV_BIN="${HIFI_ENV_BIN:-}"
[[ -n "$HIFI_ENV_BIN" && -d "$HIFI_ENV_BIN" ]] && export PATH="${HIFI_ENV_BIN}:$PATH"

###########################
# 配置（允许 env 覆盖）
###########################
BASE_OUT="${BASE_OUT:-${PWD}/work}"
DATA_DIR="${DATA_DIR:-${BASE_OUT}/data}"
SAMPLE_LIST="${SAMPLE_LIST:-${BASE_OUT}/samples.list}"

# 每个 sample 下匹配 cell BAM 的 glob 模式（支持多个 cell）
# SAMPLE01 Revio: SAMPLE01_PacBio-Revio_m84039_*.hifi_reads.bam
BAM_GLOB="${BAM_GLOB:-*.hifi_reads.bam}"

REF_FASTA="${REF_FASTA:-${BASE_OUT}/reference/GRCh38.fasta}"
REF_ROOT="${REF_ROOT:-$(dirname "${REF_FASTA}")}"

IMAGE="${IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1}"

DOCKER_GPUS="${DOCKER_GPUS:-all}"
DOCKER_USER="${DOCKER_USER:-$(id -u):$(id -g)}"

# Read group（PacBio HiFi）
RG_PL="${RG_PL:-PACBIO}"
RG_LB="${RG_LB:-HiFi}"

# 中间合并 BAM 是否保留（默认删除，节省空间）
KEEP_MERGED_UBAM="${KEEP_MERGED_UBAM:-0}"

# 是否启用 --pbmm2（PacBio HiFi 推荐处理）
USE_PBMM2_OPT="${USE_PBMM2_OPT:-1}"

###########################
# 输出目录
###########################
ALIGNED_DIR="${BASE_OUT}/aligned_bam"
TMP_DIR="${BASE_OUT}/tmp"
LOG_DIR="${BASE_OUT}/01_minimap2/log"

mkdir -p "$ALIGNED_DIR" "$TMP_DIR" "$LOG_DIR"

FAIL_TSV="${LOG_DIR}/step1_failed.tsv"
[[ -f "$FAIL_TSV" ]] || echo -e "sample\tstep\treason\ttime" > "$FAIL_TSV"

###########################
# 工具函数（沿用体细胞 step1）
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
  printf '[%s] END    %-14s elapsed=%-7s (%s)\n' \
    "$(date '+%F %T')" "$1" "${elapsed}s" "$(fmt_time ${elapsed})" >> "${TIMING}"
}

log_fail() {
  local sample="$1" step="$2" reason="$3"
  echo -e "${sample}\t${step}\t${reason}\t$(date '+%F %T')" >> "$FAIL_TSV"
  printf '[%s] FAILED %s reason=%s\n' "$(date '+%F %T')" "$step" "$reason" >> "${TIMING}"
}

docker_gpus_arg() {
  local g="$1"
  if [[ -z "$g" || "$g" == "all" ]]; then echo "all"
  elif [[ "$g" =~ ^[0-9]+(,[0-9]+)*$ ]]; then echo "device=${g}"
  elif [[ "$g" =~ ^device=.*$ ]]; then echo "$g"
  else echo "$g"; fi
}
DOCKER_GPUS_FMT="$(docker_gpus_arg "$DOCKER_GPUS")"

HAS_SAMTOOLS=0
command -v samtools >/dev/null 2>&1 && HAS_SAMTOOLS=1

###########################
# 检查
###########################
[[ -s "${SAMPLE_LIST}" ]] || { echo "[FATAL] 样本列表不存在：${SAMPLE_LIST}" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]] || { echo "[FATAL] 参考缺失：${REF_FASTA}" >&2; exit 1; }
[[ "${HAS_SAMTOOLS}" -eq 1 ]] || { echo "[FATAL] samtools 不可用（合并 cell BAM 需要）" >&2; exit 1; }

# 容器内参考路径
REF_IN="/refdir/$(basename "$REF_FASTA")"

sed -i 's/\r$//' "${SAMPLE_LIST}" 2>/dev/null || true

###########################
# 头信息
###########################
echo "=========================================="
echo "HiFi GPU Step1: Parabricks minimap2"
echo "BASE_OUT         : ${BASE_OUT}"
echo "DATA_DIR         : ${DATA_DIR}"
echo "SAMPLE_LIST      : ${SAMPLE_LIST}"
echo "BAM_GLOB         : ${BAM_GLOB}"
echo "REF_FASTA        : ${REF_FASTA}"
echo "IMAGE            : ${IMAGE}"
echo "DOCKER_GPUS      : ${DOCKER_GPUS} (docker --gpus ${DOCKER_GPUS_FMT})"
echo "USE_PBMM2_OPT    : ${USE_PBMM2_OPT}  (1=加 --pbmm2 匹配 pbmm2 输出)"
echo "RG               : PL=${RG_PL} LB=${RG_LB}"
echo "LOG_DIR          : ${LOG_DIR}"
echo "=========================================="

###########################
# 主循环：逐样本
###########################
while read -r sample; do
  [[ -z "${sample}" ]] && continue

  SAMPLE_CUR="${sample}"
  TIMING="${LOG_DIR}/${sample}.timing.log"
  log_file="${LOG_DIR}/${sample}.log"

  echo
  echo ">>> 处理样本: ${sample}    开始: $(date)"

  out_bam="${ALIGNED_DIR}/${sample}.minimap2.bam"
  out_bai="${ALIGNED_DIR}/${sample}.minimap2.bam.bai"
  merged_ubam="${TMP_DIR}/${sample}.unaligned.merged.bam"

  # 断点续跑：已有 BAM+BAI 且 quickcheck 通过则跳过
  if [[ -f "${out_bam}" && -f "${out_bai}" ]]; then
    if samtools quickcheck -q "${out_bam}" >/dev/null 2>&1; then
      echo ">>> 跳过（已存在且 quickcheck 通过）: ${sample}"
      : > "${TIMING}"
      printf '[%s] SKIP   BAM already exists and quickcheck OK\n' "$(date '+%F %T')" >> "${TIMING}"
      printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${sample}" >> "${TIMING}"
      continue
    else
      echo ">>> 警告：BAM 存在但 quickcheck 失败，删除重跑: ${sample}" >&2
      rm -f "${out_bam}" "${out_bai}" 2>/dev/null || true
    fi
  fi

  # 本次运行 reset timing
  : > "${TIMING}"

  # ---- 收集该样本的所有 cell unaligned BAM ----
  # 安全策略：只匹配属于该 ${sample} 的 BAM，绝不扫描整个 data/ 目录
  #   优先级 A: ${DATA_DIR}/${sample}/${BAM_GLOB}          (子目录方式，推荐)
  #   兜底  B: ${DATA_DIR}/${sample}${BAM_GLOB}            (平铺+样本名前缀)
  #   兜底  B': ${DATA_DIR}/${sample}_${BAM_GLOB}          (平铺+样本名_前缀)
  # 三种都以 ${sample} 限定，不会误合并其他样本的数据
  shopt -s nullglob
  cell_bams=()

  # A) 子目录方式（多样本必须用这个，天然隔离）
  cell_bams=( "${DATA_DIR}/${sample}/"${BAM_GLOB} )

  # B) 平铺但文件名以样本名开头（兼容 SAMPLE01 历史布局）
  if [[ ${#cell_bams[@]} -eq 0 ]]; then
    cell_bams=( "${DATA_DIR}/${sample}"${BAM_GLOB} )
    # 若上面因为 glob 拼接（sample + *xxx）匹配范围过宽，再补一个更严格的 sample_ 前缀
    if [[ ${#cell_bams[@]} -eq 0 ]]; then
      cell_bams=( "${DATA_DIR}/${sample}_"${BAM_GLOB} )
    fi
  fi
  shopt -u nullglob

  if [[ ${#cell_bams[@]} -eq 0 ]]; then
    echo "[WARN] 找不到 ${sample} 的 unaligned BAM" >&2
    echo "       尝试过的位置:" >&2
    echo "         A) ${DATA_DIR}/${sample}/${BAM_GLOB}" >&2
    echo "         B) ${DATA_DIR}/${sample}${BAM_GLOB}" >&2
    echo "         B') ${DATA_DIR}/${sample}_${BAM_GLOB}" >&2
    log_fail "${sample}" "input" "no_cell_bam_found"
    continue
  fi

  # ---- 防呆检查：确认匹配到的 BAM 文件名都真的属于该样本 ----
  # 防止平铺场景下样本名前缀撞车（如 HG00 误匹配 SAMPLE01）
  bad_match=0
  for b in "${cell_bams[@]}"; do
    bname="$(basename "$b")"
    parent="$(basename "$(dirname "$b")")"
    # 合法情形：父目录 == sample（子目录方式），或 文件名以 sample 开头（平铺方式）
    if [[ "$parent" == "$sample" ]]; then
      continue
    elif [[ "$bname" == "${sample}"* ]]; then
      continue
    else
      echo "[FATAL] ${sample}: 匹配到不属于该样本的 BAM: $b" >&2
      echo "        这通常是样本名前缀撞车。请用子目录方式 data/${sample}/ 隔离。" >&2
      bad_match=1
    fi
  done
  if [[ ${bad_match} -eq 1 ]]; then
    log_fail "${sample}" "input" "sample_prefix_collision"
    continue
  fi

  echo ">>> 发现 ${#cell_bams[@]} 个 cell BAM (样本 ${sample}):"
  for b in "${cell_bams[@]}"; do echo "    - $b"; done

  # ---- 步骤 1: 合并多个 cell BAM (samtools cat) ----
  # samtools cat 适合 unaligned BAM 串接，比 samtools merge 快很多（不排序、不重写）
  if [[ ${#cell_bams[@]} -eq 1 ]]; then
    # 只有一个 cell，直接用，不需要合并
    in_ubam="${cell_bams[0]}"
    echo ">>> 单 cell BAM，跳过合并: $in_ubam"
    printf '[%s] SKIP   merge_cells (single cell BAM)\n' "$(date '+%F %T')" >> "${TIMING}"
  else
    # 多个 cell，合并
    in_ubam="${merged_ubam}"
    if [[ -s "${merged_ubam}" ]] && samtools quickcheck -q "${merged_ubam}" >/dev/null 2>&1; then
      echo ">>> 已合并的 unaligned BAM 存在，跳过 samtools cat"
      printf '[%s] SKIP   merge_cells (merged BAM exists)\n' "$(date '+%F %T')" >> "${TIMING}"
    else
      step_start merge_cells
      if ! samtools cat -o "${merged_ubam}" "${cell_bams[@]}" >> "${log_file}" 2>&1; then
        log_fail "${sample}" "merge_cells" "samtools_cat_failed"
        rm -f "${merged_ubam}" 2>/dev/null || true
        continue
      fi
      step_end merge_cells
    fi
  fi

  # ---- 步骤 2: Parabricks minimap2 ----
  step_start minimap2

  # 把 in_ubam 所在目录挂载进容器（如果是 merged ubam 在 tmp 下，已经在 BASE_OUT 内）
  in_ubam_dir="$(dirname "$in_ubam")"
  in_ubam_basename="$(basename "$in_ubam")"

  # 决定容器内输入路径：如果输入 BAM 在 BASE_OUT 下，复用 /workdir 挂载；否则额外挂载
  if [[ "$in_ubam_dir" == "$BASE_OUT"* ]]; then
    in_ubam_container="/workdir/${in_ubam#${BASE_OUT}/}"
    extra_mount=""
  else
    in_ubam_container="/inubam/${in_ubam_basename}"
    extra_mount="-v ${in_ubam_dir}:/inubam"
  fi

  # pbmm2 选项（HiFi 推荐）
  pbmm2_arg=""
  [[ "${USE_PBMM2_OPT}" == "1" ]] && pbmm2_arg="--pbmm2"

  if ! docker run --gpus "${DOCKER_GPUS_FMT}" --user "${DOCKER_USER}" --rm \
        --workdir /workdir \
        -v "${REF_ROOT}":/refdir \
        -v "${BASE_OUT}":/workdir \
        ${extra_mount} \
        "${IMAGE}" pbrun minimap2 \
          --ref "${REF_IN}" \
          --in-bam "${in_ubam_container}" \
          ${pbmm2_arg} \
          --read-group-sm "${sample}" \
          --read-group-id-prefix "${sample}" \
          --read-group-pl "${RG_PL}" \
          --read-group-lb "${RG_LB}" \
          --out-bam "/workdir/aligned_bam/${sample}.minimap2.bam" \
        >> "${log_file}" 2>&1; then
    log_fail "${sample}" "minimap2" "minimap2_exit_nonzero"
    continue
  fi
  step_end minimap2

  # ---- 步骤 3: 索引（pbrun 应该已经生成，但保险起见 fall-back）----
  if [[ ! -f "${out_bai}" ]]; then
    step_start samtools_index
    if ! samtools index -@ 8 "${out_bam}" >> "${log_file}" 2>&1; then
      log_fail "${sample}" "samtools_index" "index_failed"
      continue
    fi
    step_end samtools_index
  fi

  # ---- 产物校验 ----
  if [[ ! -f "${out_bam}" || ! -f "${out_bai}" ]]; then
    log_fail "${sample}" "minimap2" "bam_or_bai_missing"
    continue
  fi
  if ! samtools quickcheck -q "${out_bam}" >/dev/null 2>&1; then
    log_fail "${sample}" "minimap2" "bam_corrupt_after_minimap2"
    rm -f "${out_bam}" "${out_bai}" 2>/dev/null || true
    continue
  fi

  echo ">>> 完成: ${sample}"
  printf '[%s] ALL DONE  %s\n' "$(date '+%F %T')" "${sample}" >> "${TIMING}"

  # ---- 清理中间合并 BAM ----
  if [[ "${KEEP_MERGED_UBAM}" -ne 1 && -f "${merged_ubam}" ]]; then
    echo ">>> 清理中间合并 BAM: ${merged_ubam}"
    rm -f "${merged_ubam}" 2>/dev/null || true
  fi

  echo ">>> 结束: $(date)"
  echo "---"
done < "${SAMPLE_LIST}"

echo
echo "=========================================="
echo "Step1 完成"
echo "  aligned_bam 目录: ${ALIGNED_DIR}"
echo "  每样本 log      : ${LOG_DIR}/<sample>.log"
echo "  每样本 timing   : ${LOG_DIR}/<sample>.timing.log"
echo "  失败清单        : ${FAIL_TSV}"
echo ""
echo "  进度排查示例:"
echo "    跑完的样本数 : grep -l 'ALL DONE' ${LOG_DIR}/*.timing.log | wc -l"
echo "    没跑完的样本 : grep -L 'ALL DONE' ${LOG_DIR}/*.timing.log"
echo "    各步耗时     : grep 'END' ${LOG_DIR}/*.timing.log"
echo "=========================================="
