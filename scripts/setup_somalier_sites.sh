#!/usr/bin/env bash
# Download a GRCh38 Somalier sites VCF with chr-prefixed contig names.
# Override SOMALIER_DIR or SITES_URL as needed.

set -euo pipefail

SOMALIER_DIR="${SOMALIER_DIR:-${PWD}/resources/somalier}"
SITES_URL="${SITES_URL:-https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz}"

mkdir -p "${SOMALIER_DIR}"
cd "${SOMALIER_DIR}"

echo "Downloading Somalier sites VCF..."
wget -c -O sites.hg38.vcf.gz "${SITES_URL}"

echo "Validating VCF..."
if zcat sites.hg38.vcf.gz 2>/dev/null | head -100 | grep -q '^##fileformat=VCF'; then
  echo "OK: VCF header detected."
else
  echo "ERROR: downloaded file does not look like a valid VCF." >&2
  exit 1
fi

echo "First three records:"
zcat sites.hg38.vcf.gz | grep -v '^#' | head -3 | awk '{print $1, $2, $3}'
echo "Sites VCF: ${SOMALIER_DIR}/sites.hg38.vcf.gz"
