# HiFi Germline Workflow

A modular workflow for PacBio HiFi germline analysis, covering alignment, small-variant and structural-variant calling, cohort-level joint genotyping, and sample identity/relatedness QC.

## Workflow

```text
unaligned HiFi BAM(s)
        |
        v
Parabricks minimap2
        |
        v
    aligned BAM
      /     \
     /       \
DeepVariant  Sniffles2
   |           |
 GVCF+VCF    VCF+SNF
   |           |
   v           v
 GLnexus   Sniffles multi-sample
   |           |
   v           v
cohort small   cohort SV
variant VCF    VCF

aligned BAM
    |
    v
Somalier extract + relate
    |
    v
sample identity / relatedness QC
```

## Repository structure

```text
hifi-germline-workflow/
├── config/
│   ├── config.example.sh
│   ├── samples.example.txt
│   └── trio.example.ped
├── scripts/
│   ├── 01_align_hifi.sh
│   ├── 02_call_small_variants.sh
│   ├── 03_call_structural_variants.sh
│   ├── 04_joint_small_variants.sh
│   ├── 05_joint_structural_variants.sh
│   ├── 06_relatedness_qc.sh
│   └── setup_somalier_sites.sh
├── .gitignore
└── README.md
```

## Inputs

- unaligned PacBio HiFi BAM file(s)
- one sample ID per line in a sample list
- GRCh38 reference FASTA
- optional PED file for expected relationships in Somalier QC

For multiple BAM files per sample, the alignment module supports a layout such as:

```text
data/
├── SAMPLE01/
│   ├── cell1.hifi_reads.bam
│   └── cell2.hifi_reads.bam
└── SAMPLE02/
    └── cell1.hifi_reads.bam
```

## Main software

- Docker
- NVIDIA Parabricks for minimap2 and DeepVariant
- samtools
- bcftools / tabix
- Sniffles2
- GLnexus
- Somalier
- SLURM for the CPU modules by default; direct execution is also supported where implemented

The container tags in `config/config.example.sh` are reproducibility examples rather than claims of current latest versions. Review software versions and resource settings for your environment before running.

## Configuration

```bash
cp config/config.example.sh config/config.sh
cp config/samples.example.txt /path/to/hifi-project/samples.list

# Edit config/config.sh, then:
source config/config.sh
```

Download Somalier sites if needed:

```bash
SOMALIER_DIR="${BASE_OUT}/resources/somalier" \
bash scripts/setup_somalier_sites.sh
```

## Run

### 1. Align HiFi reads

```bash
source config/config.sh
bash scripts/01_align_hifi.sh
```

The module accepts one or multiple unaligned HiFi BAMs per sample and produces an aligned BAM plus index.

### 2. Call small variants

```bash
source config/config.sh
bash scripts/02_call_small_variants.sh
```

Produces per-sample DeepVariant GVCFs and variant-only VCFs.

### 3. Call structural variants

```bash
source config/config.sh
bash scripts/03_call_structural_variants.sh
```

Produces per-sample Sniffles VCFs and SNF files. The default mode submits per-sample jobs through SLURM.

### 4. Joint-call small variants

```bash
source config/config.sh
bash scripts/04_joint_small_variants.sh
```

Combines DeepVariant GVCFs into a cohort VCF with GLnexus.

### 5. Joint-call structural variants

```bash
source config/config.sh
bash scripts/05_joint_structural_variants.sh
```

Combines per-sample SNF files with Sniffles multi-sample calling.

### 6. Sample identity / relatedness QC

Without a pedigree:

```bash
source config/config.sh
RUN_MODE=direct bash scripts/06_relatedness_qc.sh
```

With a pedigree:

```bash
source config/config.sh
PED_FILE=config/trio.example.ped \
RUN_MODE=direct \
bash scripts/06_relatedness_qc.sh
```

For larger cohorts, the QC script can submit itself as a SLURM job by using `RUN_MODE=slurm`.

## Notes

- Reference files, sequencing data, and software installations are not included.
- Paths, scheduler settings, GPU allocation, memory, and thread counts should be adapted to the local environment.
- The workflow contains both single-sample and cohort-level calling modules.
- Mendelian consistency and truth-set benchmarking are intentionally kept outside this repository so that routine germline calling and benchmark evaluation remain separate concerns.
