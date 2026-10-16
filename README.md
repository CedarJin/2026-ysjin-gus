# Active GUS abundance in human stool metagenomes

This project quantifies genes matching the **active β-glucuronidase (GUS) reference set** in 18 human stool metagenomes from BioProject PRJNA816918. Active GUS abundance is reported as **RPKM: reads per kilobase of coding sequence per million mapped reads**.

The working analysis requires **≥95% amino-acid identity and ≥90% query and reference coverage** (`pident95_cov90`). A match identifies a homolog of an experimentally active reference protein; it does not establish that the sample-derived protein has the same activity.

## Samples

- 18 paired-end Illumina runs: `SRR18354890`–`SRR18354907`.
- Subject IDs are `S##`, mapped to run accessions using `reference_from_paper/SraRunTable.csv` and `config/PRJNA816918_runs.tsv`.
- 17 subjects have corresponding enzymatic assays; S02 is excluded from assay comparisons.
- Assay display order: **S11, S19, S18, S09, S07, S15, S05, S03, S20, S16, S13, S04, S08, S14, S10, S17, S21**. This is the figure order, not a strict activity ranking.

## Active GUS references

The active set comprises the first 11 entries of Table S3 in Zhang et al. (2026), *Gut Microbiota Mediate the Metabolism of Colonic Prostaglandins*. Reference identifiers are retained as supplied in the analysis.

| Reference ID | Source label |
| --- | --- |
| `NP_416134.1` | E. coli |
| `WP_003467686.1` | C. perfringens |
| `6MVG_A` | R. gnavus |
| `WP_022501306` | Eubacterium |
| `WP_000966715` | S. agalactiae |
| `WP_012740861` | E. eligens |
| `3CMG_A` | B. fragilis |
| `6NCY_A` | F. saccharivorans |
| `WP_270534885` | B. mass |
| `WP_007841259` | B. dorei |
| `WP_217761416` | P. copri |

Reference sequences and annotations are stored locally in `results/reference/active_gus.faa` and `results/reference/active_gus_metadata.tsv`. Eight of these references are detected at the working thresholds. Source labels describe the reference proteins, not independently confirmed species assignments for each sample-derived gene.

## Analysis workflow

```text
Paired-end reads
  → quality control and fastp filtering
  → host removal against T2T-CHM13
  → MEGAHIT assembly (minimum contig length: 1,000 bp)
  → Prodigal gene prediction in metagenomic mode
  → Bowtie2 mapping to each sample's own assembly
  → featureCounts read counts for predicted CDS
  → DIAMOND protein hits restricted to the 11 active GUS references
  → identity/coverage filtering and one best hit per predicted gene
  → CDS-level RPKM
  → summed RPKM per active GUS reference and subject
```

Active hits are extracted from the existing protein-search results. Filtering uses amino-acid identity ≥95%, alignment length / query length ≥90%, alignment length / reference length ≥90%, and e-value ≤1e-5. Among passing active-reference hits, the highest-bitscore hit is retained for each predicted gene to avoid counting the same CDS more than once.

The shared preprocessing and mapping workflow is in `scripts/cleandata.smk`, `scripts/assemble_megahit_cleandata.smk`, and `scripts/gus_abundance.smk`. Active-reference extraction is implemented in `scripts/compute_active_gus_abundance.py`.

## RPKM calculation

For each detected active GUS coding sequence *g* in sample *s*:

```math
\mathrm{RPKM}_{g,s} = \frac{\mathrm{read\ count}_{g,s} \times 10^9}{\mathrm{CDS\ length}_{g}\;(\mathrm{bp}) \times \mathrm{mapped\ reads}_{s}}
```

The calculation uses:

- **Read count:** the original featureCounts v2.1.1 count for the CDS. The run used `-p -B -C -t CDS -g ID` without `--countReadPairs`, so these are individual-read counts, not fragment counts.
- **Length:** nucleotide CDS length from featureCounts, checked against the existing CDS-length records. Protein length in amino acids is not used for RPKM normalization.
- **Mapped reads:** all reads mapped to that sample's assembly, derived consistently from the Bowtie2 log as twice the input pair count minus unaligned individual mates. This includes mapped reads that were not assigned to a CDS.

The denominator is not the sum of active GUS counts or the number of reads assigned to genes. featureCounts summary totals are retained for auditing; the normalization denominator comes from Bowtie2.

Unrounded CDS-level RPKM values are summed for genes assigned to the same reference, then written to the matrix with six decimal places. A subject's **summed active GUS abundance** is the sum across the detected active references. The matrix is not renormalized to a fixed column sum, and RPKM is not a percentage of the genome.

## Active GUS output files

The current RPKM outputs are under:

```text
results/summary/pident95_cov90/abundance_all_subject/
```

| File | Contents |
| --- | --- |
| `active_gus_abundance_gus_genome_all_subject_sample_id.tsv` | Active GUS reference × subject matrix, in **RPKM**; 8 detected references × 18 subjects |
| `active_gus_rpkm_gene_details.tsv` | CDS-level raw counts, nucleotide lengths, mapped-read denominators, and RPKM |
| `active_gus_rpkm_library_sizes.tsv` | Per-subject library sizes and normalization audit |
| `active_gus_rpkm_methods.txt` | Calculation details and provenance |

The `gus_genome` portion of the matrix filename is retained for compatibility; **the values in this specific matrix are RPKM**. Other legacy abundance tables and figures were not converted with this matrix and should not be assumed to use RPKM. Regenerating older merged tables can overwrite this file, so the mapped-read RPKM calculation must be applied again before downstream use.

For stacked abundance plots, use subject on the x-axis, **GUS abundance in RPKM** on the y-axis, and a separate color for each detected active GUS reference. Assay comparisons use the 17-subject order above.

## Interpretation

These abundances measure metagenomic DNA representation, not protein expression or catalytic activity. Matching the same reference does not guarantee identical sample-derived protein sequences, and an individual reference's abundance cannot be equated with whole-sample enzymatic activity.

Comparisons to the enzymatic assay should distinguish original measurements from values estimated from a figure. Sequence substitutions, predicted truncations, and reference expression tags also require separate interpretation before proposing functional effects.

## Repository layout

```text
config/                 Run accessions and workflow configuration
scripts/                Preprocessing, assembly, mapping, and analysis workflows
reference_from_paper/   Paper reference sequences, Table S3, and sample metadata
results/                Local search, count, and abundance outputs (gitignored)
assembly/               Sample assemblies (gitignored)
cleandata/              Host-removed reads (gitignored)
```

Large input and result files are stored locally and are not included in the GitHub repository. `environment.yml` records the mapping and gene-quantification software environment.
