# Active GUS abundance in human stool metagenomes

This project quantifies genes matching the **active β-glucuronidase (GUS) reference set** in 18 human stool metagenomes from BioProject PRJNA816918. 

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
  → Abundance calculation
```

Active hits are extracted from the existing protein-search results. Filtering uses amino-acid identity ≥95%, alignment length / query length ≥90%, alignment length / reference length ≥90%, and e-value ≤1e-5. Among passing active-reference hits, the highest-bitscore hit is retained for each predicted gene to avoid counting the same CDS more than once.

The shared preprocessing and mapping workflow is in `scripts/cleandata.smk`, `scripts/assemble_megahit_cleandata.smk`, and `scripts/gus_abundance.smk`. Active-reference extraction is implemented in `scripts/compute_active_gus_abundance.py`.

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
