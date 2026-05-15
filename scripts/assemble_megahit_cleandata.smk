# Megahit: cleandata/<id>/{id}_{{1,2}}.fastq -> assembly/megahit/<id>/final.contigs.fa
# sbatch scripts/submit_assemble_megahit_cleandata.sbatch
#   snakemake -s scripts/assemble_megahit_cleandata.smk -c 8 -n

import glob
import os

from snakemake.utils import min_version

min_version("7.0")

CLEAN = "cleandata"
OUT = "assembly/megahit"
LOG = "logs/megahit_cleandata"

minlen = int(config.get("min_contig_len", 1000))
threads = int(config.get("megahit_threads", 8))

samples = sorted(
    os.path.basename(d)
    for d in glob.glob(os.path.join(CLEAN, "*"))
    if os.path.isdir(d)
    and os.path.isfile(os.path.join(d, f"{os.path.basename(d)}_1.fastq"))
    and os.path.isfile(os.path.join(d, f"{os.path.basename(d)}_2.fastq"))
)

if not samples:
    raise WorkflowError(f"no paired FASTQ under {CLEAN}/*/")

wildcard_constraints:
    sample=r"[^/]+",

shell.prefix("set -euo pipefail; ")


rule all:
    input:
        expand(f"{OUT}/{{sample}}/final.contigs.fa", sample=samples),


rule megahit:
    input:
        r1=f"{CLEAN}/{{sample}}/{{sample}}_1.fastq",
        r2=f"{CLEAN}/{{sample}}/{{sample}}_2.fastq",
    output:
        f"{OUT}/{{sample}}/final.contigs.fa",
    log:
        f"{LOG}/{{sample}}.log",
    threads:
        threads
    params:
        odir=lambda wc: os.path.join(OUT, wc.sample),
        mlen=minlen,
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        rm -rf "{params.odir}"
        megahit -1 {input.r1} -2 {input.r2} -o "{params.odir}" \
            --min-contig-len {params.mlen} -t {threads} > {log} 2>&1
        """
