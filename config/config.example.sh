#!/usr/bin/env bash
# Example environment configuration.
# Copy to config/config.sh, edit values, then:
#   source config/config.sh

export BASE_OUT="/path/to/hifi-project"
export DATA_DIR="${BASE_OUT}/data"
export SAMPLE_LIST="${BASE_OUT}/samples.list"

export REF_FASTA="/path/to/reference/GRCh38.fasta"
export REF_ROOT="$(dirname "${REF_FASTA}")"

# Optional environment bootstrap script. Leave empty if tools are already on PATH.
export ENV_FILE=""
export HIFI_ENV_BIN=""

# Container images used by the workflow. Pin versions appropriate for your environment.
export IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"
export GLNEXUS_IMAGE="ghcr.io/dnanexus-rnd/glnexus:v1.4.1"

# Tools available on PATH by default; absolute paths may also be supplied.
export BCFTOOLS="bcftools"
export TABIX="tabix"
export SNIFFLES="sniffles"
export SOMALIER="somalier"

# Somalier sites VCF.
export SITES_VCF="${BASE_OUT}/resources/somalier/sites.hg38.vcf.gz"

# Scheduler defaults for CPU modules; adjust to your cluster.
export QUEUE="compute"
export CORE="16"
export MEM="32"
export TIME_LIMIT="24:00:00"

# Cohort label used by joint-calling/QC modules.
export COHORT_NAME="example_cohort"
