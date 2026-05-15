# cleandata.smk
#
# Per-sample QC + host removal workflow for paired-end FASTQ in rawdata/.
#
# Steps:
#   1. FastQC on raw reads
#   2. MultiQC on raw FastQC reports
#   3. fastp (default parameters + --detect_adapter_for_pe)
#   4. FastQC on fastp-trimmed reads
#   5. MultiQC on post-fastp FastQC reports
#   6. Bowtie2 host removal vs. T2T-CHM13 (existing prebuilt index)
#
# Inputs:
#   rawdata/{sample}_1.fastq
#   rawdata/{sample}_2.fastq
#
# Outputs:
#   trimmed_data/{sample}_1.fastq         (temporary; auto-deleted after
#   trimmed_data/{sample}_2.fastq          fastqc_post + remove_host succeed)
#   trimmed_data/{sample}.fastp.html / .fastp.json  (kept)
#   cleandata/{sample}/{sample}_1.fastq   (host-removed reads, unpaired-conc out)
#   cleandata/{sample}/{sample}_2.fastq
#   qc/multiqc_raw/multiqc_raw_report.html
#   qc/multiqc_post/multiqc_post_report.html
#
# Note: trimmed fastqs are marked temp() so Snakemake deletes them only
# after every downstream consumer (fastqc_post AND remove_host) has finished
# successfully.  If any downstream step fails, the trimmed data stays so
# you can resume without re-running fastp.
#
# Run from repo root with the `qc` conda env active:
#   snakemake -s scripts/cleandata.smk -c 8 -p
# Dry run:
#   snakemake -s scripts/cleandata.smk -n

import glob
import os


# ---------------------------------------------------------------------------
# Paths / parameters
# ---------------------------------------------------------------------------

RAWDATA       = "rawdata"
TRIMMED_DATA  = "trimmed_data"
CLEAN_DATA    = "cleandata"
QC_DIR        = "qc"
LOG_DIR       = "logs/cleandata"

HOST_REF_DIR      = "/home/jys0914/2026-ysjin-seqdep/reference/human_T2T"
HOST_INDEX_PREFIX = f"{HOST_REF_DIR}/GCF_009914755.1_T2T-CHM13v2.0_genomic"
INDEX_EXTS        = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
HOST_INDEX_FILES  = [f"{HOST_INDEX_PREFIX}.{ext}" for ext in INDEX_EXTS]

FASTQC_THREADS  = 2
FASTP_THREADS   = 4
BOWTIE2_THREADS = 8

READS = ["1", "2"]


# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------

def discover_samples():
    samples = []
    for r1 in sorted(glob.glob(os.path.join(RAWDATA, "*_1.fastq"))):
        r2 = r1.replace("_1.fastq", "_2.fastq")
        if not os.path.exists(r2):
            continue
        samples.append(os.path.basename(r1)[:-len("_1.fastq")])
    return sorted(set(samples))


SAMPLES = discover_samples()
if not SAMPLES:
    raise WorkflowError(
        f"No paired *_1.fastq / *_2.fastq files found in {RAWDATA}/."
    )


wildcard_constraints:
    sample=r"[^/]+",
    r=r"[12]",


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

rule all:
    input:
        f"{QC_DIR}/multiqc_raw/multiqc_raw_report.html",
        f"{QC_DIR}/multiqc_post/multiqc_post_report.html",
        expand(f"{CLEAN_DATA}/{{sample}}/{{sample}}_{{r}}.fastq",
               sample=SAMPLES, r=READS),


# ---------------------------------------------------------------------------
# Step 1: FastQC on raw reads
# ---------------------------------------------------------------------------

rule fastqc_raw:
    input:
        f"{RAWDATA}/{{sample}}_{{r}}.fastq",
    output:
        html=f"{QC_DIR}/fastqc_raw/{{sample}}_{{r}}_fastqc.html",
        zip =f"{QC_DIR}/fastqc_raw/{{sample}}_{{r}}_fastqc.zip",
    log:
        f"{LOG_DIR}/fastqc_raw/{{sample}}_{{r}}.log",
    threads:
        FASTQC_THREADS
    shell:
        r"""
        mkdir -p {QC_DIR}/fastqc_raw $(dirname {log})
        fastqc -t {threads} -o {QC_DIR}/fastqc_raw {input} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Step 2: MultiQC over raw FastQC
# ---------------------------------------------------------------------------

rule multiqc_raw:
    input:
        expand(f"{QC_DIR}/fastqc_raw/{{sample}}_{{r}}_fastqc.zip",
               sample=SAMPLES, r=READS),
    output:
        report=f"{QC_DIR}/multiqc_raw/multiqc_raw_report.html",
    log:
        f"{LOG_DIR}/multiqc_raw.log",
    shell:
        r"""
        mkdir -p {QC_DIR}/multiqc_raw $(dirname {log})
        multiqc {QC_DIR}/fastqc_raw \
            -o {QC_DIR}/multiqc_raw \
            -n multiqc_raw_report --force > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Step 3: fastp (default parameters)
# ---------------------------------------------------------------------------

rule fastp:
    input:
        r1=f"{RAWDATA}/{{sample}}_1.fastq",
        r2=f"{RAWDATA}/{{sample}}_2.fastq",
    output:
        r1  =temp(f"{TRIMMED_DATA}/{{sample}}_1.fastq"),
        r2  =temp(f"{TRIMMED_DATA}/{{sample}}_2.fastq"),
        html=f"{TRIMMED_DATA}/{{sample}}.fastp.html",
        json=f"{TRIMMED_DATA}/{{sample}}.fastp.json",
    log:
        f"{LOG_DIR}/fastp/{{sample}}.log",
    threads:
        FASTP_THREADS
    shell:
        r"""
        mkdir -p {TRIMMED_DATA} $(dirname {log})
        fastp -w {threads} \
            --detect_adapter_for_pe \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Step 4: FastQC on fastp-trimmed reads
# ---------------------------------------------------------------------------

rule fastqc_post:
    input:
        f"{TRIMMED_DATA}/{{sample}}_{{r}}.fastq",
    output:
        html=f"{QC_DIR}/fastqc_post/{{sample}}_{{r}}_fastqc.html",
        zip =f"{QC_DIR}/fastqc_post/{{sample}}_{{r}}_fastqc.zip",
    log:
        f"{LOG_DIR}/fastqc_post/{{sample}}_{{r}}.log",
    threads:
        FASTQC_THREADS
    shell:
        r"""
        mkdir -p {QC_DIR}/fastqc_post $(dirname {log})
        fastqc -t {threads} -o {QC_DIR}/fastqc_post {input} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Step 5: MultiQC over post-fastp FastQC
# ---------------------------------------------------------------------------

rule multiqc_post:
    input:
        expand(f"{QC_DIR}/fastqc_post/{{sample}}_{{r}}_fastqc.zip",
               sample=SAMPLES, r=READS),
    output:
        report=f"{QC_DIR}/multiqc_post/multiqc_post_report.html",
    log:
        f"{LOG_DIR}/multiqc_post.log",
    shell:
        r"""
        mkdir -p {QC_DIR}/multiqc_post $(dirname {log})
        multiqc {QC_DIR}/fastqc_post \
            -o {QC_DIR}/multiqc_post \
            -n multiqc_post_report --force > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Step 6: Bowtie2 host removal vs T2T-CHM13
# ---------------------------------------------------------------------------
#
# Uses the prebuilt index at HOST_INDEX_PREFIX.  --un-conc writes the
# concordantly-unaligned (i.e. non-host) reads, with %% replaced by 1/2.

rule remove_host:
    input:
        r1=f"{TRIMMED_DATA}/{{sample}}_1.fastq",
        r2=f"{TRIMMED_DATA}/{{sample}}_2.fastq",
        idx=HOST_INDEX_FILES,
    output:
        r1=f"{CLEAN_DATA}/{{sample}}/{{sample}}_1.fastq",
        r2=f"{CLEAN_DATA}/{{sample}}/{{sample}}_2.fastq",
    params:
        idx_prefix=HOST_INDEX_PREFIX,
        out_prefix=f"{CLEAN_DATA}/{{sample}}/{{sample}}_%.fastq",
    log:
        f"{LOG_DIR}/bowtie2/{{sample}}.log",
    threads:
        BOWTIE2_THREADS
    shell:
        r"""
        mkdir -p {CLEAN_DATA}/{wildcards.sample} $(dirname {log})
        bowtie2 -x {params.idx_prefix} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            --very-sensitive-local \
            --un-conc {params.out_prefix} \
            -S /dev/null > {log} 2>&1
        """
