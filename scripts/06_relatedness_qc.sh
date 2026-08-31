#!/usr/bin/env bash
# Public portfolio version: paths, scheduler names, sample IDs, and project-specific defaults are parameterized.
# ============================================================================
# HiFi 流程 Step3.5: 样本身份 / 亲缘关系 QC (somalier)
#
# 定位：通用sample identity / relatedness QC，applicable to cohort-level QC（肿瘤/人群/家系/出生缺陷均适用）。
#       目的是抓样本错配、样本重复、污染、意外亲缘 —— 保证下游分析用的
#       是身份正确的样本。
#
# 本步只做 somalier，不做 Mendelian。
#   - somalier 亲缘鉴定：不预设关系，从数据推断样本间亲缘 → 抓样本身份问题
#   - Mendelian consistency validation is intentionally kept outside this workflow repository.
#
# PED 处理（supported workflow）：
#   - 传 PED（-p）时，somalier 自己从家系结构推导每对样本的预期亲缘系数，
#     输出到 pairs.tsv 的 expected_relatedness 列。本脚本直接对比
#     "实测 relatedness vs 预期 expected_relatedness"，吻合则 OK，偏差大则
#     flag（如声称亲子却实测无关 → 样本错配/非亲生）。不自己推导关系。
#   - 不传 PED 时，退回纯阈值判断（数据驱动，抓意外亲缘/样本重复）。
#
# 说明：
#   - somalier 从 BAM extract（recommended，比从 GVCF 更准）
#   - somalier 用 sites.hg38.vcf.gz (带 chr，匹配 GRCh38.fasta)
#   - HiFi 亲缘 QC 用 somalier 是supported workflow（somalier 论文即用 pbmm2+HiFi）
#
# 输入：
#   ${BASE_OUT}/aligned_bam/*.minimap2.bam   (somalier extract 用)
#   ${PED_FILE}  (可选，家系文件，用于 somalier 预期关系判断)
#
# 输出：
#   ${BASE_OUT}/qc/somalier/${COHORT_NAME}.pairs.tsv    样本两两亲缘系数
#   ${BASE_OUT}/qc/somalier/${COHORT_NAME}.samples.tsv  每样本信息
#   ${BASE_OUT}/qc/somalier/${COHORT_NAME}.html         可视化报告
#
# 运行方式（重要）：
#   extract 每样本约 20-30 秒（扫 HiFi BAM 采样）。几百样本串行 = 数小时。
#   登录节点(<submit-node>)会杀掉长时间进程，所以默认把整个串行流程包进一个
#   SLURM 作业里跑（单作业串行，非多作业并行——因 extract 快、无需并行，
#   且串行天然容错：一个样本失败不影响其余，断点续跑可接续）。
#
#   RUN_MODE=slurm  (默认)：自动生成 sbatch 包装器提交到队列，作业内串行跑。
#   RUN_MODE=direct        ：前台直接串行跑（仅建议小数据测试，如 3 样本 trio）。
#
# 容错与失败处理：
#   - extract 断点续跑：已有 .somalier 的样本自动跳过。
#   - extract 失败清单：某样本 BAM 缺失/损坏/提取失败 → 记入
#       ${LOG_DIR}/step3_5_failed.tsv（沿用 step1/2a 的失败清单传统）。
#   - relate 缺失报告：relate 前统计实际成功的 .somalier 数，报告哪些样本
#       因 extract 失败被排除，用现有成功样本继续跑（个别失败不卡死整体）。
#   - 亲缘异常清单：relate 后把"实测 vs 预期不符"的异常样本对单独记入
#       ${SOM_OUT_DIR}/${COHORT_NAME}.anomalies.tsv（供追查样本来源）。
#
# 用法：
#   # 生产（有 PED，投队列串行跑）:
#   BASE_OUT=... SAMPLE_LIST=... COHORT_NAME=... PED_FILE=... \
#   bash 06_relatedness_qc.sh
#
#   # 小数据测试（前台跑）:
#   RUN_MODE=direct BASE_OUT=... SAMPLE_LIST=... COHORT_NAME=... PED_FILE=... \
#   bash 06_relatedness_qc.sh
#
#   # 无 PED（普通队列，纯阈值抓样本身份问题）:
#   BASE_OUT=... SAMPLE_LIST=... COHORT_NAME=cohort1 \
#   bash 06_relatedness_qc.sh
# ============================================================================
set -uo pipefail

###########################
# 配置
###########################
BASE_OUT="${BASE_OUT:-${PWD}/work}"
SAMPLE_LIST="${SAMPLE_LIST:-${BASE_OUT}/samples.list}"
COHORT_NAME="${COHORT_NAME:-cohort}"
ALIGNED_DIR="${ALIGNED_DIR:-${BASE_OUT}/aligned_bam}"

REF_FASTA="${REF_FASTA:-${BASE_OUT}/reference/GRCh38.fasta}"
REF_ROOT="${REF_ROOT:-$(dirname "${REF_FASTA}")}"

# somalier
SOMALIER="${SOMALIER:-somalier}"
SITES_VCF="${SITES_VCF:-${BASE_OUT}/resources/somalier/sites.hg38.vcf.gz}"

# PED 家系文件（可选）。给了就做"实测 vs 预期"符合性判断，不给就纯阈值推断
PED_FILE="${PED_FILE:-}"

# 运行模式：slurm（默认，包进 sbatch 作业串行跑）/ direct（前台串行，测试用）
RUN_MODE="${RUN_MODE:-slurm}"

# SLURM 资源（RUN_MODE=slurm 时用）
QUEUE="${QUEUE:-compute}"
CORE="${CORE:-4}"
MEM="${MEM:-16}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"

###########################
# 输出目录
###########################
SOM_EXTRACT_DIR="${BASE_OUT}/qc/somalier/extracted"
SOM_OUT_DIR="${BASE_OUT}/qc/somalier"
LOG_DIR="${BASE_OUT}/log/step3_5"

mkdir -p "$SOM_EXTRACT_DIR" "$SOM_OUT_DIR" "$LOG_DIR"

# --------------------------------------------------------------------
# 运行模式分发（自我投递）
#   RUN_MODE=slurm 且当前不在 SLURM 作业内 → 生成 sbatch 包装器，把本脚本
#   以 RUN_MODE=direct 重新提交到队列（作业内串行执行工作逻辑）。
#   已在作业内（SLURM_JOB_ID 存在）或 RUN_MODE=direct → 直接往下执行。
# --------------------------------------------------------------------
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
if [[ "${RUN_MODE}" == "slurm" && -z "${SLURM_JOB_ID:-}" ]]; then
  SH_DIR="${BASE_OUT}/qc/somalier/sh"
  mkdir -p "${SH_DIR}"
  WRAP="${SH_DIR}/${COHORT_NAME}.step3_5.slurm.sh"
  {
    echo "#!/usr/bin/env bash"
    echo "#SBATCH -p ${QUEUE}"
    echo "#SBATCH -J som_${COHORT_NAME}"
    echo "#SBATCH -c ${CORE}"
    echo "#SBATCH --mem=${MEM}G"
    echo "#SBATCH -t ${TIME_LIMIT}"
    echo "#SBATCH -o ${LOG_DIR}/${COHORT_NAME}.slurm.log"
    echo "#SBATCH -e ${LOG_DIR}/${COHORT_NAME}.slurm.log"
    echo ""
    # 作业内以 direct 模式跑本脚本，透传所有配置
    echo "RUN_MODE=direct \\"
    echo "BASE_OUT='${BASE_OUT}' \\"
    echo "SAMPLE_LIST='${SAMPLE_LIST}' \\"
    echo "COHORT_NAME='${COHORT_NAME}' \\"
    echo "ALIGNED_DIR='${ALIGNED_DIR}' \\"
    echo "REF_ROOT='${REF_ROOT}' \\"
    echo "REF_FASTA='${REF_FASTA}' \\"
    echo "SOMALIER='${SOMALIER}' \\"
    echo "SITES_VCF='${SITES_VCF}' \\"
    echo "PED_FILE='${PED_FILE}' \\"
    echo "bash '${SCRIPT_PATH}'"
  } > "${WRAP}"

  jobid=$(sbatch "${WRAP}" | awk '{print $NF}')
  echo "=========================================="
  echo "Step3.5 已提交到 SLURM（单作业串行）"
  echo "  JobID     : ${jobid}"
  echo "  队列/资源 : ${QUEUE}  ${CORE}核 ${MEM}G  ${TIME_LIMIT}"
  echo "  包装脚本  : ${WRAP}"
  echo "  作业日志  : ${LOG_DIR}/${COHORT_NAME}.slurm.log"
  echo "  查看进度  : squeue -u \$USER    tail -f ${LOG_DIR}/${COHORT_NAME}.slurm.log"
  echo ""
  echo "  提示：小数据测试可用 RUN_MODE=direct 前台直接跑。"
  echo "=========================================="
  exit 0
fi

# ===== 以下为实际工作逻辑（RUN_MODE=direct 或已在 SLURM 作业内执行）=====

LOG_FILE="${LOG_DIR}/${COHORT_NAME}.log"
: > "${LOG_FILE}"

# extract 失败清单（沿用 step1/2a 传统）
FAIL_TSV="${LOG_DIR}/step3_5_failed.tsv"
[[ -f "$FAIL_TSV" ]] || echo -e "sample\tstep\treason\ttime" > "$FAIL_TSV"

# 亲缘异常清单（step3.5 特有：只记 FAIL 的样本对）
ANOMALY_TSV="${SOM_OUT_DIR}/${COHORT_NAME}.anomalies.tsv"

fmt_time() { local s=$1; printf '%dh%02dm%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }

###########################
# 检查
###########################
[[ -s "${SAMPLE_LIST}" ]] || { echo "[FATAL] 样本列表不存在：${SAMPLE_LIST}" >&2; exit 1; }
command -v "${SOMALIER}" >/dev/null 2>&1 || { echo "[FATAL] somalier not found: ${SOMALIER}" >&2; exit 1; }
[[ -s "${SITES_VCF}" ]]   || { echo "[FATAL] sites VCF 不存在：${SITES_VCF}（先跑 somalier_setup_sites.sh）" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]]   || { echo "[FATAL] 参考缺失：${REF_FASTA}" >&2; exit 1; }
[[ -d "${ALIGNED_DIR}" ]] || { echo "[FATAL] BAM 目录不存在：${ALIGNED_DIR}" >&2; exit 1; }

sed -i 's/\r$//' "${SAMPLE_LIST}" 2>/dev/null || true
mapfile -t samples < "${SAMPLE_LIST}"

echo "=========================================="
echo "HiFi Step3.5: 样本关系 QC"
echo "COHORT_NAME    : ${COHORT_NAME}"
echo "样本数         : ${#samples[@]}"
echo "SITES_VCF      : ${SITES_VCF}"
echo "PED_FILE       : ${PED_FILE:-(无，跳过 Mendelian)}"
echo "=========================================="

###########################
# 1. somalier extract（每个样本从 BAM 提取，recommended）
###########################
echo ""
echo "[1] somalier extract (从 BAM)..."
t0=$(date +%s)

n_ok=0; n_skip=0; n_fail=0
for sample in "${samples[@]}"; do
  [[ -z "${sample}" ]] && continue
  bam="${ALIGNED_DIR}/${sample}.minimap2.bam"

  if [[ ! -f "${bam}" ]]; then
    echo "  ✗ ${sample}: BAM 不存在，跳过：${bam}" | tee -a "${LOG_FILE}" >&2
    echo -e "${sample}\textract\tbam_not_found\t$(date '+%F %T')" >> "$FAIL_TSV"
    ((n_fail++)) || true
    continue
  fi

  # 已提取则跳过（断点续跑）
  if [[ -f "${SOM_EXTRACT_DIR}/${sample}.somalier" ]]; then
    echo "  ✓ ${sample}: 已提取，跳过"
    ((n_skip++)) || true
    continue
  fi

  if ${SOMALIER} extract \
       -d "${SOM_EXTRACT_DIR}/" \
       --sites "${SITES_VCF}" \
       -f "${REF_FASTA}" \
       "${bam}" >> "${LOG_FILE}" 2>&1; then
    echo "  ✓ ${sample}: 提取完成"
    ((n_ok++)) || true
  else
    echo "  ✗ ${sample}: 提取失败，查看 ${LOG_FILE}" | tee -a "${LOG_FILE}" >&2
    echo -e "${sample}\textract\tsomalier_extract_failed\t$(date '+%F %T')" >> "$FAIL_TSV"
    ((n_fail++)) || true
  fi
done

t1=$(date +%s)
echo "  extract 完成: 新提取 ${n_ok} / 跳过 ${n_skip} / 失败 ${n_fail}，耗时 $(fmt_time $((t1-t0)))"
if [[ ${n_fail} -gt 0 ]]; then
  echo "  ⚠ ${n_fail} 个样本 extract 失败，已记入失败清单: ${FAIL_TSV}" | tee -a "${LOG_FILE}"
fi

###########################
# 2. somalier relate（亲缘关系推断）
###########################
echo ""
echo "[2] somalier relate (亲缘推断)..."

# 收集 .somalier 文件
mapfile -t som_files < <(ls "${SOM_EXTRACT_DIR}"/*.somalier 2>/dev/null)

# 缺失样本报告：对照样本清单，找出哪些样本没有 .somalier（extract 失败/未完成）
missing_samples=()
for sample in "${samples[@]}"; do
  [[ -z "${sample}" ]] && continue
  [[ -f "${SOM_EXTRACT_DIR}/${sample}.somalier" ]] || missing_samples+=("${sample}")
done
echo "  参与亲缘分析的样本: ${#som_files[@]} / ${#samples[@]}"
if [[ ${#missing_samples[@]} -gt 0 ]]; then
  echo "  ⚠ 以下 ${#missing_samples[@]} 个样本因 extract 失败/缺失被排除（详见 ${FAIL_TSV}）:" | tee -a "${LOG_FILE}"
  printf '      %s\n' "${missing_samples[@]}"
  echo "    → 用现有 ${#som_files[@]} 个样本继续（个别缺失不影响整体亲缘分析）"
fi
echo ""

if [[ ${#som_files[@]} -lt 2 ]]; then
  echo "  ✗ 成功提取的样本 < 2，无法算亲缘。" | tee -a "${LOG_FILE}" >&2
else
  # 若有 PED，传给 somalier 作为预期关系（会在报告里标注 expected vs observed）
  # 若有 PED，传给 somalier 作为预期关系（somalier 会算出 expected_relatedness 列）
  # somalier 0.3.2 用短参 -p 更稳（长参 --ped 部分版本不认）
  ped_arg=""
  [[ -n "${PED_FILE}" && -f "${PED_FILE}" ]] && ped_arg="-p ${PED_FILE}"

  if ${SOMALIER} relate ${ped_arg} \
       -o "${SOM_OUT_DIR}/${COHORT_NAME}" \
       "${som_files[@]}" >> "${LOG_FILE}" 2>&1; then
    echo "  ✓ relate 完成"
    echo ""

    PAIRS_TSV="${SOM_OUT_DIR}/${COHORT_NAME}.pairs.tsv"
    if [[ -f "${PAIRS_TSV}" ]]; then
      # --------------------------------------------------------------
      # 判读逻辑（PED 感知，supported workflow）
      #
      # somalier relate 传入 -p PED 后，pairs.tsv 会输出 expected_relatedness
      # 列（第 17 列）——somalier 自己从 PED 家系结构推导每对样本的预期亲缘
      # 系数（亲子/同胞=0.5，配偶/无关=0，PED 未指定=-1）。
      # 我们不自己推导关系，直接对比：
      #   relatedness (第 3 列, 实测)  vs  expected_relatedness (第 17 列, 预期)
      #
      # 有 PED（expected>=0）：按"实测 vs 预期"是否吻合判断 → 抓样本错配/非亲生
      # 无 PED（expected=-1）：退回纯阈值判断（一级/二级/无关/异常）
      #
      # 列索引：1=sample_a 2=sample_b 3=relatedness 17=expected_relatedness
      # --------------------------------------------------------------
      # 动态定位列号（防止不同 somalier 版本列顺序变化）
      HDR=$(head -1 "${PAIRS_TSV}")
      COL_REL=$(echo "${HDR}" | tr '\t' '\n' | grep -nx 'relatedness' | cut -d: -f1)
      COL_EXP=$(echo "${HDR}" | tr '\t' '\n' | grep -nx 'expected_relatedness' | cut -d: -f1)
      COL_REL="${COL_REL:-3}"
      COL_EXP="${COL_EXP:-17}"

      anomaly_found=0

      echo "  === 样本两两亲缘关系 QC (${COHORT_NAME}.pairs.tsv) ==="
      if [[ -n "${PED_FILE}" && -f "${PED_FILE}" ]]; then
        echo "  模式: 有 PED → 对比实测 vs 预期（抓样本错配/非亲生）"
      else
        echo "  模式: 无 PED → 纯阈值推断（数据驱动，抓意外亲缘/样本重复）"
      fi
      echo ""

      # --------------------------------------------------------------
      # 生成"同家系样本对"清单，解决 founder（父母）对误报问题。
      # somalier 对 PED 里的 founder 对（如父×母，各自父母都是0）给
      # expected_relatedness = -1（无预期关系），而非 0。若按无 PED 的纯阈值
      # 判断，父母对的负 relatedness 会被误报为"负值异常"。
      # 因此：同属一个 FamilyID 的样本对，即使 expected=-1，也视为"预期无关"
      #      （founder 对），只在实测异常高时才 flag，负值/低值视为正常。
      # SAME_FAM 文件每行 "sampleA\tsampleB"（字典序），供 awk 查表。
      # --------------------------------------------------------------
      SAME_FAM="$(mktemp)"
      if [[ -n "${PED_FILE}" && -f "${PED_FILE}" ]]; then
        # 剥注释/空行/CR，取 FamilyID(列1) + IID(列2)，同 family 两两配对
        sed 's/\r$//' "${PED_FILE}" | awk 'NF>0 && $1 !~ /^#/ {fam[$2]=$1; order[NR]=$2}
          END{
            n=0; for(i in order) samp[++n]=order[i]
            for(i=1;i<=n;i++) for(j=i+1;j<=n;j++){
              a=samp[i]; b=samp[j]
              if(fam[a]==fam[b]){
                if(a<b) print a"\t"b; else print b"\t"a
              }
            }
          }' > "${SAME_FAM}" 2>/dev/null || true
      fi

      # 初始化异常清单 TSV（只记 FAIL 的样本对）
      echo -e "sample_a\tsample_b\tmeasured_rel\texpected_rel\treason" > "${ANOMALY_TSV}"

      # 用 awk 逐对判读；把异常行数通过退出码之外的方式回传（写临时标记文件）
      FLAG_FILE="$(mktemp)"
      # 先把同家系对读进 awk 的关联数组
      tail -n +2 "${PAIRS_TSV}" | awk -F'\t' -v cr="${COL_REL}" -v ce="${COL_EXP}" -v flag="${FLAG_FILE}" -v samefam="${SAME_FAM}" -v anom="${ANOMALY_TSV}" '
        BEGIN{
          # 载入同家系对（key = "a\tb" 字典序）
          while((getline line < samefam) > 0){ fampair[line]=1 }
          close(samefam)
        }
        function samefam_pair(a,b,   k){
          if(a<b) k=a"\t"b; else k=b"\t"a
          return (k in fampair)
        }
        function relclass(r){
          if (r >= 0.35)  return "一级亲缘(亲子/同胞)"
          if (r >= 0.15)  return "二级亲缘(祖孙/叔侄/半同胞)"
          if (r >= 0.08)  return "三级亲缘/远亲"
          if (r <= -0.10) return "负值异常(检查数据质量)"
          return "无关"
        }
        {
          a=$1; b=$2; rel=$cr+0; erel=$ce+0
          observed=relclass(rel)

          if (erel < 0) {
            # expected=-1：somalier 未给预期关系
            if (samefam_pair(a,b)) {
              # 但两样本同属一个 PED 家系 → founder 对（如父×母），预期无关。
              # 只在实测异常高（≥0.20，暗示样本重复/污染）时才 flag；
              # 负值/低值都是正常的（父母本就无血缘）。
              if (rel >= 0.20) {
                printf "    %-12s %-12s rel=%-7.3f (同家系founder,预期无关)  ✗ 实测偏高 → 可能样本重复/污染\n", a, b, rel
                print "ANOMALY" >> flag
                printf "%s\t%s\t%.3f\t%s\t%s\n", a, b, rel, "0(founder)", "声称无关(同家系founder)但实测偏高→可能样本重复/污染" >> anom
              } else {
                printf "    %-12s %-12s rel=%-7.3f (同家系founder,预期无关)  ✓ 符合  [%s]\n", a, b, rel, observed
              }
            } else {
              # 真正无 PED 预期（跨家系或无 PED）：纯阈值判断
              mark="→"
              if (rel <= -0.10) {
                mark="✗"; print "ANOMALY" >> flag
                printf "%s\t%s\t%.3f\t%s\t%s\n", a, b, rel, "NA(无预期)", "无PED预期下负值relatedness→检查数据质量/样本污染" >> anom
              }
              printf "    %-12s %-12s rel=%-7.3f  %s %s\n", a, b, rel, mark, observed
            }
          } else {
            # 有 PED 预期：对比实测 vs 预期
            # 预期一级(0.5 附近)但实测明显偏低 → 错配/非亲生
            # 预期无关(0 附近)但实测明显偏高 → 意外亲缘/样本重复
            ok=1; reason=""
            if (erel >= 0.35) {
              # 预期一级
              if (rel < 0.20) { ok=0; reason="声称一级亲缘但实测偏低 → 可能样本错配/非亲生" }
            } else if (erel <= 0.10) {
              # 预期无关
              if (rel >= 0.20) { ok=0; reason="声称无关但实测偏高 → 可能意外亲缘/样本重复" }
            } else {
              # 预期二级等中间值
              if (rel < erel-0.20 || rel > erel+0.20) { ok=0; reason="实测与预期偏差较大" }
            }
            if (ok) {
              printf "    %-12s %-12s rel=%-7.3f (预期%.2f)  ✓ 符合  [%s]\n", a, b, rel, erel, observed
            } else {
              printf "    %-12s %-12s rel=%-7.3f (预期%.2f)  ✗ %s\n", a, b, rel, erel, reason
              print "ANOMALY" >> flag
              printf "%s\t%s\t%.3f\t%.2f\t%s\n", a, b, rel, erel, reason >> anom
            }
          }
        }'

      if [[ -s "${FLAG_FILE}" ]]; then
        anomaly_found=$(wc -l < "${FLAG_FILE}")
        echo ""
        echo "  ⚠ 发现 ${anomaly_found} 对异常关系，请核查样本身份/配对！" | tee -a "${LOG_FILE}"
        echo "    异常清单: ${ANOMALY_TSV}" | tee -a "${LOG_FILE}"
      else
        echo ""
        echo "  ✓ 所有样本对关系正常，无样本错配/污染信号"
        # 无异常，删掉只有表头的空清单
        rm -f "${ANOMALY_TSV}"
      fi
      rm -f "${FLAG_FILE}" "${SAME_FAM}"
    fi
  else
    echo "  ✗ relate 失败，查看 ${LOG_FILE}" | tee -a "${LOG_FILE}" >&2
  fi
fi


echo ""
echo "=========================================="
echo "Step3.5 完成（somalier 样本身份/亲缘 QC）"
echo "  somalier 报告  : ${SOM_OUT_DIR}/${COHORT_NAME}.html"
echo "  亲缘 pairs     : ${SOM_OUT_DIR}/${COHORT_NAME}.pairs.tsv"
echo "  每样本信息     : ${SOM_OUT_DIR}/${COHORT_NAME}.samples.tsv"
echo "  extract 失败   : ${FAIL_TSV}  (若有失败样本，重跑本脚本会自动补跑)"
if [[ -f "${ANOMALY_TSV}" ]]; then
  echo "  亲缘异常清单   : ${ANOMALY_TSV}  ⚠ 有异常样本对，需核查来源"
fi
echo "  日志           : ${LOG_FILE}"
echo ""
echo "  说明: 本步是通用样本身份 QC（任何队列都跑），用于抓样本错配/污染/意外亲缘。"
echo "        Mendelian 一致性验证（拿金标准 trio 测流程质量）已拆为独立脚本:"
echo "        a separate Mendelian validation workflow"
echo "=========================================="
