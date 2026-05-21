# Metagenomic GUS abundance workflow:
# reads -> bowtie2 -> featureCounts -> prodigal -> diamond vs total/candidate GUS -> abundance normalization

import glob
import os
from pathlib import Path

from snakemake.utils import min_version

min_version("7.0")

# Resolved relative to repo root (run snakemake from project directory).
configfile: "config/gus_abundance.yaml"

shell.prefix("set -euo pipefail; ")

CONFIG = config
DIAMOND = CONFIG.get(
    "diamond",
    {"evalue": 1e-5, "max_target_seqs": 25},
)
READ_DIR = CONFIG.get("read_dir", "cleandata")
ASSEMBLY_DIR = CONFIG.get("assembly_dir", "assembly/megahit")
RESULTS = CONFIG.get("results_dir", "results")
SAMPLE = CONFIG.get("sample")
MINLEN = int(CONFIG.get("min_contig_length", 1000))
THREADS_BOWTIE2 = int(CONFIG.get("threads_bowtie2", CONFIG.get("threads", 8)))
THREADS_FEATURECOUNTS = int(CONFIG.get("threads_featurecounts", CONFIG.get("threads", 8)))
THREADS_PRODIGAL = int(CONFIG.get("threads_prodigal", 8))
THREADS_DIAMOND = int(CONFIG.get("threads_diamond", 8))


def paired_samples():
    read_samples = {
        os.path.basename(d)
        for d in glob.glob(os.path.join(READ_DIR, "*"))
        if os.path.isdir(d)
        and os.path.isfile(os.path.join(d, f"{os.path.basename(d)}_1.fastq"))
        and os.path.isfile(os.path.join(d, f"{os.path.basename(d)}_2.fastq"))
    }
    assembly_samples = {
        os.path.basename(d)
        for d in glob.glob(os.path.join(ASSEMBLY_DIR, "*"))
        if os.path.isdir(d)
        and os.path.isfile(os.path.join(d, "final.contigs.fa"))
    }
    samples = sorted(read_samples & assembly_samples)
    if SAMPLE:
        samples = [s for s in samples if s == SAMPLE]
    if not samples:
        raise WorkflowError(f"no paired FASTQ under {READ_DIR}/*/")
    return samples


SAMPLES = paired_samples()
def assembly_contigs(sample):
    return os.path.join(CONFIG["assembly_dir"], sample, CONFIG["contigs_basename"])


rule all:
    input:
        expand(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_abundance.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_abundance.tsv", sample=SAMPLES),


rule bowtie2_index:
    input:
        contigs=lambda wc: assembly_contigs(wc.sample),
    output:
        idx1=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.1.bt2"),
        idx2=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.2.bt2"),
        idx3=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.3.bt2"),
        idx4=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.4.bt2"),
        idxrev1=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.rev.1.bt2"),
        idxrev2=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.rev.2.bt2"),
    threads: THREADS_BOWTIE2
    shell:
        """
        mkdir -p {RESULTS}/mapping
        bowtie2-build {input} {RESULTS}/mapping/{wildcards.sample}.min{MINLEN} > /dev/null 2>&1
        """


rule bowtie2_align:
    input:
        ref=lambda wc: assembly_contigs(wc.sample),
        idx=rules.bowtie2_index.output,
        r1=lambda wc: f"{READ_DIR}/{wc.sample}/{wc.sample}_1.fastq",
        r2=lambda wc: f"{READ_DIR}/{wc.sample}/{wc.sample}_2.fastq",
    output:
        bam=f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.sorted.bam",
    threads: THREADS_BOWTIE2
    shell:
        """
        bowtie2 -x {RESULTS}/mapping/{wildcards.sample}.min{MINLEN} \
            -1 {input.r1} -2 {input.r2} -p {threads} 2> {RESULTS}/mapping/{wildcards.sample}.min{MINLEN}.bowtie2.log \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule prodigal:
    input:
        contigs=lambda wc: assembly_contigs(wc.sample),
    output:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        ffn=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.ffn",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    threads: THREADS_PRODIGAL
    shell:
        """
        mkdir -p {RESULTS}/prodigal
        prodigal -i {input.contigs} -p meta -q \
            -a {output.faa} -d {output.ffn} -o {output.gff} -f gff
        """


rule featurecounts:
    input:
        bam=f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.sorted.bam",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        counts=temp(f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt"),
    threads: THREADS_FEATURECOUNTS
    shell:
        """
        mkdir -p {RESULTS}/counts
        featureCounts -T {threads} -p -B -C -t CDS -g ID -a {input.gff} -o {output.counts} {input.bam}
        """


rule diamond_makedb_total:
    input:
        "results/reference/gus_279_plus_fmn.faa",
    output:
        "results/reference/gus_279_plus_fmn.dmnd",
    shell:
        """
        diamond makedb --in {input} --db results/reference/gus_279_plus_fmn
        """


rule diamond_makedb_candidate:
    input:
        "results/reference/candidate_gus.faa",
    output:
        "results/reference/candidate_gus.dmnd",
    shell:
        """
        diamond makedb --in {input} --db results/reference/candidate_gus
        """


rule diamond_blastp_total:
    input:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        db="results/reference/gus_279_plus_fmn.dmnd",
    output:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.total_gus.tsv",
    params:
        db_prefix="results/reference/gus_279_plus_fmn",
        evalue=DIAMOND["evalue"],
        max_target=DIAMOND["max_target_seqs"],
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {RESULTS}/search
        diamond blastp --query {input.faa} --db {params.db_prefix} --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} --max-target-seqs {params.max_target} --threads {threads}
        """


rule diamond_blastp_candidate:
    input:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        db="results/reference/candidate_gus.dmnd",
    output:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.candidate_gus.tsv",
    params:
        db_prefix="results/reference/candidate_gus",
        evalue=DIAMOND["evalue"],
        max_target=DIAMOND["max_target_seqs"],
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {RESULTS}/search
        diamond blastp --query {input.faa} --db {params.db_prefix} --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} --max-target-seqs {params.max_target} --threads {threads}
        """


rule summarize_total_hits:
    input:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.total_gus.tsv",
        contigs=lambda wc: assembly_contigs(wc.sample),
        metadata="results/reference/gus_279_plus_fmn_metadata.tsv",
    output:
        best=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_best_hits.tsv"),
        clade=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_clade_abundance.tsv"),
    shell:
        """
        python3 scripts/summarize_gus_hits.py \
            --hits {input.hits} --contigs {input.contigs} --metadata {input.metadata} \
            --best-hits {output.best} --clade-summary {output.clade}
        """


rule summarize_candidate_hits:
    input:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.candidate_gus.tsv",
        contigs=lambda wc: assembly_contigs(wc.sample),
        metadata="results/reference/candidate_gus_metadata.tsv",
    output:
        best=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_best_hits.tsv"),
        clade=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_clade_abundance.tsv"),
    shell:
        """
        python3 scripts/summarize_gus_hits.py \
            --hits {input.hits} --contigs {input.contigs} --metadata {input.metadata} \
            --best-hits {output.best} --clade-summary {output.clade}
        """


rule abundance_total:
    input:
        hits=f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_best_hits.tsv",
        counts=f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt",
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_abundance.tsv",
    shell:
        """
        python3 scripts/calculate_gus_abundance.py \
            --hits {input.hits} --counts {input.counts} --faa {input.faa} --gff {input.gff} \
            --output {output}
        """


rule abundance_candidate:
    input:
        hits=f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_best_hits.tsv",
        counts=f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt",
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_abundance.tsv",
    shell:
        """
        python3 scripts/calculate_gus_abundance.py \
            --hits {input.hits} --counts {input.counts} --faa {input.faa} --gff {input.gff} \
            --output {output}
        """
