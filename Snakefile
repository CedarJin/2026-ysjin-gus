# Run from the repository root, with prodigal, diamond, and Python (+pypdf) on PATH.
# Example: snakemake -c 8 --configfile config/config.yaml
# Dry run: snakemake -n

import os
from pathlib import Path

from snakemake.utils import min_version

min_version("7.0")

configfile: "config/config.yaml"

Path("logs").mkdir(parents=True, exist_ok=True)

shell.prefix("set -euo pipefail; ")

ASSEMBLY_CONTIGS = os.path.join(config["assembly_dir"], config["contigs_basename"])


rule all:
    input:
        expand(
            "results/summary/{sample}.min{minlen}.gus_best_hits.tsv",
            sample=[config["sample"]],
            minlen=[config["min_contig_length"]],
        ),
        expand(
            "results/summary/{sample}.min{minlen}.gus_clade_abundance.tsv",
            sample=[config["sample"]],
            minlen=[config["min_contig_length"]],
        ),


rule extract_gus_reference:
    input:
        pdf=config["reference_pdf"],
    output:
        fasta="results/reference/gus_279.faa",
        metadata="results/reference/gus_279_metadata.tsv",
    params:
        expected=config["expected_gus_records"],
    log:
        "logs/extract_gus_reference.log",
    shell:
        """
        python scripts/extract_gus_reference.py \
            --input {input.pdf} \
            --fasta {output.fasta} \
            --metadata {output.metadata} \
            --expected-count {params.expected} \
            > {log} 2>&1
        """


rule filter_contigs:
    input:
        contigs=ASSEMBLY_CONTIGS,
    output:
        "results/filtered_contigs/{sample}.min{minlen}.fa",
    wildcard_constraints:
        sample="[^/]+",
        minlen=r"\d+",
    log:
        "logs/filter_{sample}_min{minlen}.log",
    shell:
        """
        python scripts/filter_contigs_by_length.py \
            --input {input.contigs} \
            --output {output} \
            --min-length {wildcards.minlen} \
            > {log} 2>&1
        """


rule prodigal:
    input:
        "results/filtered_contigs/{sample}.min{minlen}.fa",
    output:
        faa="results/prodigal/{sample}.min{minlen}.faa",
        ffn="results/prodigal/{sample}.min{minlen}.ffn",
        gff="results/prodigal/{sample}.min{minlen}.gff",
    log:
        "logs/prodigal_{sample}_min{minlen}.log",
    shell:
        """
        prodigal -i {input} -p meta -q \
            -a {output.faa} -d {output.ffn} -o {output.gff} -f gff \
            > {log} 2>&1
        """


rule diamond_makedb:
    input:
        "results/reference/gus_279.faa",
    output:
        "results/reference/gus_279.dmnd",
    log:
        "logs/diamond_makedb.log",
    shell:
        """
        diamond makedb --in {input} --db results/reference/gus_279 > {log} 2>&1
        """


rule diamond_blastp:
    input:
        query="results/prodigal/{sample}.min{minlen}.faa",
        db="results/reference/gus_279.dmnd",
    output:
        "results/search/{sample}.min{minlen}.predicted_proteins_vs_gus.tsv",
    params:
        db_prefix="results/reference/gus_279",
        evalue=config["diamond"]["evalue"],
        max_target=config["diamond"]["max_target_seqs"],
    threads: int(config["threads_diamond"])
    log:
        "logs/diamond_blastp_{sample}_min{minlen}.log",
    shell:
        """
        diamond blastp \
            --query {input.query} \
            --db {params.db_prefix} \
            --out {output} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} \
            --max-target-seqs {params.max_target} \
            --threads {threads} \
            > {log} 2>&1
        """


rule summarize_gus_hits:
    input:
        hits="results/search/{sample}.min{minlen}.predicted_proteins_vs_gus.tsv",
        contigs="results/filtered_contigs/{sample}.min{minlen}.fa",
        metadata="results/reference/gus_279_metadata.tsv",
    output:
        best="results/summary/{sample}.min{minlen}.gus_best_hits.tsv",
        clade="results/summary/{sample}.min{minlen}.gus_clade_abundance.tsv",
    params:
        min_identity=config["summarize"]["min_identity"],
        min_query_cov=config["summarize"]["min_query_cov"],
        min_subject_cov=config["summarize"]["min_subject_cov"],
        max_evalue=config["summarize"]["max_evalue"],
    log:
        "logs/summarize_{sample}_min{minlen}.log",
    shell:
        """
        python scripts/summarize_gus_hits.py \
            --hits {input.hits} \
            --contigs {input.contigs} \
            --metadata {input.metadata} \
            --best-hits {output.best} \
            --clade-summary {output.clade} \
            --min-identity {params.min_identity} \
            --min-query-cov {params.min_query_cov} \
            --min-subject-cov {params.min_subject_cov} \
            --max-evalue {params.max_evalue} \
            > {log} 2>&1
        """
