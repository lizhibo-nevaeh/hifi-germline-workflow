#!/usr/bin/env bash
# Public portfolio version: paths, scheduler names, sample IDs, and project-specific defaults are parameterized.
# ============================================================================
# HiFi 流程 Step3b: SV 多样本联合分型 (Sniffles2 multi-sample)
#
# 功能：
#   把 step2b 输出的多个 .snf 文件联合分型为 cohort SV VCF。
#
# 原理：
#   Sniffles2 的 multi-sample 模式：先每个样本出 .snf（step2b 已做），
#   然后把所有 .snf 一起喂给 sniffles 做 population-level SV 联合分型。
#   等价于 small variant 的 GLnexus，但 SV 层用 sniffles 自己的机制。
#
# 输入：
#   ${BASE_OUT}/sv_snf/*.sv.snf   (step2b 产物)
#
# 输出：
#   ${BASE_OUT}/cohort/${COHORT_NAME}.sv.vcf.gz  + .tbi
#
# 说明：
#   - Sniffles multi-sample 是纯 CPU，不需要 docker，任何节点都能跑
#   - 可在 <submit-node> 直接跑，或投 SLURM 队列
#
# 投递方式（SLURM 队列，推荐）：
#   cd /path/to/hifi-project/sh
#   BASE_OUT=/path/to/hifi-project \
#   COHORT_NAME=cohort1 \
#   QUEUE=compute CORE=16 MEM=32 \
#   bash 05_joint_structural_variants.sh
#
# 或直接前台跑（样本少时）：
#   BASE_OUT=... COHORT_NAME=cohort1 RUN_MODE=direct bash 05_joint_structural_variants.sh
#
# 进度排查：
#   squeue -u $USER
#   ${BASE_OUT}/log/step3/${COHORT_NAME}.sniffles_multi.log
# ============================================================================
set -uo pipefail

###########################
# 环境
###########################
ENV_FILE="${ENV_FILE:-}"
[[ -n "$ENV_FILE" && -s "$ENV_FILE" ]] && source "$ENV_FILE"

###########################
# 配置
###########################
BASE_OUT="${BASE_OUT:-${PWD}/work}"
COHORT_NAME="${COHORT_NAME:-cohort}"

SV_SNF_DIR="${SV_SNF_DIR:-${BASE_OUT}/sv_snf}"
COHORT_DIR="${COHORT_DIR:-${BASE_OUT}/cohort}"
LOG_DIR="${LOG_DIR:-${BASE_OUT}/log/step3}"
SH_DIR="${SH_DIR:-${BASE_OUT}/03_sniffles_multi/sh}"

REF_FASTA="${REF_FASTA:-${BASE_OUT}/reference/GRCh38.fasta}"
REF_ROOT="${REF_ROOT:-$(dirname "${REF_FASTA}")}"

SNIFFLES="${SNIFFLES:-sniffles}"
TABIX="${TABIX:-tabix}"

# 运行模式：slurm（默认，投队列）或 direct（前台直接跑）
RUN_MODE="${RUN_MODE:-slurm}"

# SLURM 资源
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-16}"
MEM="${MEM:-32}"
TIME_LIMIT="${TIME_LIMIT:-24:00:00}"

###########################
# 输出目录
###########################
mkdir -p "$COHORT_DIR" "$LOG_DIR" "$SH_DIR"

###########################
# 检查
###########################
[[ -d "${SV_SNF_DIR}" ]] || { echo "[FATAL] SNF 目录不存在：${SV_SNF_DIR}" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]]  || { echo "[FATAL] 参考缺失：${REF_FASTA}" >&2; exit 1; }
command -v "${SNIFFLES}" >/dev/null 2>&1 || { echo "[FATAL] sniffles not found: ${SNIFFLES}" >&2; exit 1; }

# 收集所有 SNF
mapfile -t snf_files < <(ls "${SV_SNF_DIR}"/*.sv.snf 2>/dev/null | sort)
n_samples=${#snf_files[@]}
[[ ${n_samples} -ge 1 ]] || { echo "[FATAL] 没找到 SNF 文件：${SV_SNF_DIR}/*.sv.snf" >&2; exit 1; }

if [[ ${n_samples} -lt 2 ]]; then
  echo "[WARN] 只有 ${n_samples} 个样本，SV 联合分型意义不大，但仍继续。" >&2
fi

OUT_VCF="${COHORT_DIR}/${COHORT_NAME}.sv.vcf.gz"
LOG_FILE="${LOG_DIR}/${COHORT_NAME}.sniffles_multi.log"
TIMING="${LOG_DIR}/${COHORT_NAME}.sv.timing.log"

# 断点续跑
if [[ -f "${OUT_VCF}" && -f "${OUT_VCF}.tbi" ]]; then
  if gunzip -c "${OUT_VCF}" | head -10 | grep -q "##fileformat=VCF"; then
    echo ">>> ${OUT_VCF} 已存在且完整，跳过"
    exit 0
  fi
fi

###########################
# 头信息
###########################
echo "=========================================="
echo "HiFi Step3b: SV 联合分型 (Sniffles multi-sample)"
echo "COHORT_NAME    : ${COHORT_NAME}"
echo "SV_SNF_DIR     : ${SV_SNF_DIR}"
echo "样本数         : ${n_samples}"
echo "OUT_VCF        : ${OUT_VCF}"
echo "RUN_MODE       : ${RUN_MODE}"
echo "=========================================="
echo "样本列表:"
for f in "${snf_files[@]}"; do echo "  $(basename "$f")"; done
echo ""

###########################
# 构造 sniffles 命令
# sniffles multi-sample: --input 后面跟所有 .snf 文件，--vcf 输出 cohort VCF
###########################
build_and_run() {
  cat > "$1" << EOF
#!/usr/bin/env bash
set -uo pipefail
TIMING="${TIMING}"
: > "\${TIMING}"
fmt_time() { local s=\$1; printf '%dh%02dm%02ds' \$((s/3600)) \$(((s%3600)/60)) \$((s%60)); }

t0=\$(date +%s)
printf '[%s] START  sniffles_multi\n' "\$(date '+%F %T')" >> "\${TIMING}"
echo "[\$(date '+%F %T')] 开始 SV 联合分型 (${n_samples} 个样本)..."

if ${SNIFFLES} \\
    --input ${snf_files[*]} \\
    --vcf "${OUT_VCF}" \\
    --reference "${REF_FASTA}" \\
    --threads ${CORE} \\
    > "${LOG_FILE}" 2>&1; then

  ${TABIX} -p vcf -f "${OUT_VCF}"

  t1=\$(date +%s); el=\$((t1-t0))
  printf '[%s] END    sniffles_multi   elapsed=%ss (%s)\n' "\$(date '+%F %T')" "\${el}" "\$(fmt_time \${el})" >> "\${TIMING}"
  printf '[%s] ALL DONE  ${COHORT_NAME}\n' "\$(date '+%F %T')" >> "\${TIMING}"
  echo "[\$(date '+%F %T')] 完成 SV 联合分型，耗时 \$(fmt_time \${el})"
  echo "  cohort SV 数: \$(gunzip -c "${OUT_VCF}" | grep -vc '^#')"
else
  printf '[%s] FAILED sniffles_multi\n' "\$(date '+%F %T')" >> "\${TIMING}"
  echo "[\$(date '+%F %T')] SV 联合分型失败，查看 ${LOG_FILE}" >&2
  exit 1
fi
EOF
}

sub_script="${SH_DIR}/${COHORT_NAME}.sniffles_multi.sh"
build_and_run "${sub_script}"

if [[ "${RUN_MODE}" == "direct" ]]; then
  echo ">>> 前台直接运行..."
  bash "${sub_script}"
  echo
  echo "=========================================="
  echo "Step3b 完成"
  echo "  cohort SV VCF : ${OUT_VCF}"
  echo "  日志          : ${LOG_FILE}"
  echo "=========================================="
else
  # SLURM 投递：给子脚本加 SBATCH 头
  slurm_script="${SH_DIR}/${COHORT_NAME}.sniffles_multi.slurm.sh"
  {
    echo "#!/usr/bin/env bash"
    echo "#SBATCH -p ${QUEUE}"
    echo "#SBATCH -J svmulti_${COHORT_NAME}"
    echo "#SBATCH -c ${CORE}"
    echo "#SBATCH --mem=${MEM}G"
    echo "#SBATCH -t ${TIME_LIMIT}"
    echo "#SBATCH -o ${LOG_DIR}/${COHORT_NAME}.sv.slurm.log"
    echo "#SBATCH -e ${LOG_DIR}/${COHORT_NAME}.sv.slurm.log"
    echo ""
    tail -n +2 "${sub_script}"   # 去掉子脚本第一行 shebang，其余接在 SBATCH 头后
  } > "${slurm_script}"

  jobid=$(sbatch "${slurm_script}" | awk '{print $NF}')
  echo ">>> 已提交 SLURM 任务 (JobID ${jobid})"
  echo
  echo "=========================================="
  echo "Step3b 已提交"
  echo "  cohort SV VCF : ${OUT_VCF} (跑完后生成)"
  echo "  查看进度      : squeue -u \$USER"
  echo "  日志          : ${LOG_FILE}"
  echo "=========================================="
fi
