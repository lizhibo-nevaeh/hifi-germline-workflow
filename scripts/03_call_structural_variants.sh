#!/usr/bin/env bash
# Public portfolio version: paths, scheduler names, sample IDs, and project-specific defaults are parameterized.
# ============================================================================
# HiFi 流程 Step2b: 结构变异检测 (Sniffles2, SLURM 队列并行)
#
# 功能：
#   每个样本生成一个 Sniffles 子脚本并 sbatch 提交到 compute 队列并行跑。
#   输出 SV VCF (单样本) + SNF (多样本联合分型用)。
#
# 说明：
#   - Sniffles 是 CPU 工具，拆到 SLURM 队列并行（不占 GPU 节点）
#   - 沿用二代 reseq 模式：主脚本 while 循环，逐样本生成子脚本 + sbatch
#   - 每个样本独立任务，可独立重跑/排错
#
# 输入：
#   ${BASE_OUT}/aligned_bam/${sample}.minimap2.bam
#
# 输出：
#   ${BASE_OUT}/sv_vcf/${sample}.sv.vcf.gz  + .tbi
#   ${BASE_OUT}/sv_snf/${sample}.sv.snf
#   ${BASE_OUT}/02_sniffles/sh/${sample}.sniffles.sh   (生成的子脚本)
#   ${BASE_OUT}/log/step2b/${sample}.log + .timing.log
#
# 投递方式（<submit-node>，SLURM）：
#   cd /path/to/hifi-project/sh
#   BASE_OUT=/path/to/hifi-project \
#   SAMPLE_LIST=/path/to/hifi-project/samples.list \
#   QUEUE=compute CORE=16 MEM=32 \
#   bash 03_call_structural_variants.sh
#
#   # 只投部分样本（试跑）:
#   bash 03_call_structural_variants.sh SAMPLE01 SAMPLE02
#
# 进度排查：
#   排队/运行 : squeue -u $USER
#   跑完的样本 : grep -l 'ALL DONE' ${BASE_OUT}/log/step2b/*.timing.log | wc -l
#   失败清单   : ${BASE_OUT}/log/step2b/step2b_failed.tsv
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

# 工具绝对路径
SNIFFLES="${SNIFFLES:-sniffles}"
TABIX="${TABIX:-tabix}"

# SLURM 资源
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-16}"
MEM="${MEM:-32}"
TIME_LIMIT="${TIME_LIMIT:-24:00:00}"

# Sniffles 参数（官方默认）
SNIFFLES_MIN_SUPPORT="${SNIFFLES_MIN_SUPPORT:-auto}"
SNIFFLES_MIN_LEN="${SNIFFLES_MIN_LEN:-35}"

###########################
# 输出目录
###########################
SV_VCF_DIR="${BASE_OUT}/sv_vcf"
SV_SNF_DIR="${BASE_OUT}/sv_snf"
SH_DIR="${BASE_OUT}/02_sniffles/sh"
LOG_DIR="${BASE_OUT}/log/step2b"

mkdir -p "$SV_VCF_DIR" "$SV_SNF_DIR" "$SH_DIR" "$LOG_DIR"

FAIL_TSV="${LOG_DIR}/step2b_failed.tsv"
[[ -f "$FAIL_TSV" ]] || echo -e "sample\treason\ttime" > "$FAIL_TSV"

###########################
# 检查
###########################
[[ -s "${SAMPLE_LIST}" ]] || { echo "[FATAL] 样本列表不存在：${SAMPLE_LIST}" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]]   || { echo "[FATAL] 参考缺失：${REF_FASTA}" >&2; exit 1; }
command -v "${SNIFFLES}" >/dev/null 2>&1 || { echo "[FATAL] sniffles not found: ${SNIFFLES}" >&2; exit 1; }

sed -i 's/\r$//' "${SAMPLE_LIST}" 2>/dev/null || true

# 允许命令行传入指定样本（试跑），否则用 SAMPLE_LIST 全量
if [[ $# -gt 0 ]]; then
  SAMPLES=("$@")
else
  mapfile -t SAMPLES < "${SAMPLE_LIST}"
fi

###########################
# 头信息
###########################
echo "=========================================="
echo "HiFi Step2b: Sniffles2 SV (SLURM 队列)"
echo "BASE_OUT       : ${BASE_OUT}"
echo "样本数         : ${#SAMPLES[@]}"
echo "QUEUE          : ${QUEUE}"
echo "CORE / MEM     : ${CORE} / ${MEM}G"
echo "SNIFFLES_PARAMS: min_support=${SNIFFLES_MIN_SUPPORT} min_len=${SNIFFLES_MIN_LEN}"
echo "SV_VCF_DIR     : ${SV_VCF_DIR}"
echo "SV_SNF_DIR     : ${SV_SNF_DIR}"
echo "=========================================="

###########################
# 逐样本生成子脚本 + sbatch
###########################
submitted_count=0
skipped_count=0
missing_count=0

for sample in "${SAMPLES[@]}"; do
  [[ -z "${sample}" ]] && continue

  in_bam="${ALIGNED_DIR}/${sample}.minimap2.bam"
  out_vcf="${SV_VCF_DIR}/${sample}.sv.vcf.gz"
  out_snf="${SV_SNF_DIR}/${sample}.sv.snf"

  # 断点续跑：VCF+SNF 都存在且 VCF 完整则跳过
  if [[ -f "${out_vcf}" && -f "${out_vcf}.tbi" && -f "${out_snf}" ]]; then
    if gunzip -c "${out_vcf}" | head -10 | grep -q "##fileformat=VCF"; then
      echo "  ✓ ${sample}: 已完成，跳过"
      ((skipped_count++)) || true
      continue
    fi
  fi

  # 输入 BAM 检查
  if [[ ! -f "${in_bam}" ]]; then
    echo "  ✗ ${sample}: 输入 BAM 不存在，跳过"
    echo -e "${sample}\tbam_not_found\t$(date '+%F %T')" >> "$FAIL_TSV"
    ((missing_count++)) || true
    continue
  fi

  # 生成子脚本
  sub_script="${SH_DIR}/${sample}.sniffles.sh"
  cat > "${sub_script}" << EOF
#!/usr/bin/env bash
#SBATCH -p ${QUEUE}
#SBATCH -J sf_${sample}
#SBATCH -c ${CORE}
#SBATCH --mem=${MEM}G
#SBATCH -t ${TIME_LIMIT}
#SBATCH -o ${LOG_DIR}/${sample}.log
#SBATCH -e ${LOG_DIR}/${sample}.log

set -uo pipefail

TIMING="${LOG_DIR}/${sample}.timing.log"
: > "\${TIMING}"

fmt_time() { local s=\$1; printf '%dh%02dm%02ds' \$((s/3600)) \$(((s%3600)/60)) \$((s%60)); }

t0=\$(date +%s)
printf '[%s] START  sniffles\n' "\$(date '+%F %T')" >> "\${TIMING}"
echo "[\$(date '+%F %T')] 开始 ${sample} sniffles..."

if ${SNIFFLES} \\
    --input "${in_bam}" \\
    --vcf "${out_vcf}" \\
    --snf "${out_snf}" \\
    --reference "${REF_FASTA}" \\
    --sample-id "${sample}" \\
    --minsupport "${SNIFFLES_MIN_SUPPORT}" \\
    --minsvlen "${SNIFFLES_MIN_LEN}" \\
    --threads ${CORE}; then

  # 索引
  ${TABIX} -p vcf -f "${out_vcf}"

  t1=\$(date +%s); el=\$((t1-t0))
  printf '[%s] END    sniffles         elapsed=%ss (%s)\n' "\$(date '+%F %T')" "\${el}" "\$(fmt_time \${el})" >> "\${TIMING}"
  printf '[%s] ALL DONE  ${sample}\n' "\$(date '+%F %T')" >> "\${TIMING}"
  echo "[\$(date '+%F %T')] 完成 ${sample}，耗时 \$(fmt_time \${el})"
else
  printf '[%s] FAILED sniffles\n' "\$(date '+%F %T')" >> "\${TIMING}"
  echo -e "${sample}\tsniffles_exit_nonzero\t\$(date '+%F %T')" >> "${FAIL_TSV}"
  echo "[\$(date '+%F %T')] 失败 ${sample}" >&2
  exit 1
fi
EOF

  # 提交
  jobid=$(sbatch "${sub_script}" | awk '{print $NF}')
  echo "  → ${sample}: 已提交 (JobID ${jobid})"
  ((submitted_count++)) || true

done

echo
echo "=========================================="
echo "Step2b 提交完成"
echo "  已提交 : ${submitted_count}"
echo "  已跳过 : ${skipped_count}"
echo "  缺 BAM : ${missing_count}"
echo ""
echo "  查看排队/运行: squeue -u \$USER"
echo "  跑完的样本  : grep -l 'ALL DONE' ${LOG_DIR}/*.timing.log | wc -l"
echo "  失败清单    : ${FAIL_TSV}"
echo ""
echo "  下一步: 多样本 SV 联合分型 step3b (用 sv_snf/*.sv.snf)"
echo "=========================================="
