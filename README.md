# Metagenomic GUS abundance (PRJNA816918)

This repository quantifies β-glucuronidase (GUS) genes in 18 human stool metagenomes from BioProject [PRJNA816918](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA816918). Host-removed reads are assembled, genes are predicted, reads are mapped back to each sample’s own assembly, predicted proteins are searched against paper GUS reference sets, and abundances are reported as **genome-wide TPM**.

The working analysis uses identity ≥ 95% and query/subject coverage ≥ 90% (`pident95_cov90`).

## Dataset

| Item | Value |
| --- | --- |
| BioProject | PRJNA816918 |
| Samples | 18 paired-end Illumina WGS runs (`SRR18354890`–`SRR18354907`) |
| Sample IDs | `S` + NCBI **Sample Name** (for example `SRR18354890` → `S13`), from `reference_from_paper/SraRunTable.csv` |
| Accession lists | `config/PRJNA816918_accessions.txt`, `config/PRJNA816918_runs.tsv` |

## GUS reference sets

| Set | Source | Size |
| --- | --- | --- |
| **Total GUS** | Pollet et al. 2017 (279) + Simpson et al. 2022 sequences with `Loop_Class == FMN` (28). No ID overlap. | **307** |
| **Candidate GUS** | Table S3 (18 experimentally discussed proteins) from the preprint (Zhang et al. 2026) | 18 |
| **Active GUS** | First 11 rows of Table S3 from the preprint (Zhang et al. 2026) | 11 |

Active GUS IDs (Table S3 order):

`NP_416134.1` (E. coli), `WP_003467686.1` (C. perfringens), `6MVG_A` (R. gnavus), `WP_022501306` (Eubacterium), `WP_000966715` (S. agalactiae), `WP_012740861` (E. eligens), `3CMG_A` (B. fragilis), `6NCY_A` (F. saccharivorans), `WP_270534885` (B. mass), `WP_007841259` (B. dorei), `WP_217761416` (P. copri).

Active abundance is computed by filtering candidate DIAMOND hits to those 11 IDs, then running the same summarize + TPM steps as total/candidate. That is equivalent to searching an 11-sequence database because the candidate search already kept up to 25 targets per query (`--max-target-seqs 25`).

### Total GUS patch (279 + 28 FMN)

The DIAMOND database (`results/reference/gus_279_plus_fmn.faa` / `.dmnd`) is still the original **408** sequences: Pollet 279 plus all 129 Simpson FMN-binding proteins. The search pipeline had labeled every Simpson sequence as clade `FMN`. In Simpson Supplemental Table 2 (`GUS_Summary`), true `Loop_Class` counts are No_Loop 42, **FMN 28**, Mini-Loop_1 22, Loop_2 14, Loop_1 9, NTL 7, Mini-Loop_2 7.

Working total-GUS abundance (`pident95_cov90` merged matrices) was restricted to the **307**-sequence set **without re-running DIAMOND**:

```bash
python3 scripts/patch_total_gus_pollet279_fmn28.py
```

The patch:

- Keeps Pollet 279 + Simpson `Loop_Class == FMN` (28).
- Filters `results/summary/pident95_cov90/abundance_all_subject/total_gus_abundance_*.tsv`.
- Renormalizes `tpm_gus_db` (each sample sums to 10<sup>6</sup>) and `within_gus_db_percent` (each sample sums to 100%).
- Leaves genome-wide TPM and `of_genome_percent` values unchanged (the denominator is the full gene catalog).
- Writes pre-patch copies to `results/summary/pident95_cov90/abundance_all_subject/_backup_before_pollet279_fmn28/`.

At this filter, **123 / 307** references were detected (112 Pollet + 11 of 28 FMN). Re-run the patch after `merge_abundance_tables.py` if those matrices are regenerated. Candidate and active GUS tables are not patched.

## Pipeline

```text
raw FASTQ
  → QC + fastp + T2T-CHM13 host removal     (scripts/cleandata.smk)
  → MEGAHIT assembly (--min-contig-len 1000) (scripts/assemble_megahit_cleandata.smk)
  → Prodigal (meta) on final.contigs.fa
  → Bowtie2: host-removed reads → same assembly
  → featureCounts on Prodigal CDS
  → DIAMOND blastp vs total and candidate GUS
  → filter hits → best hit per predicted gene
  → genome-wide TPM and GUS-set TPM           (scripts/calculate_gus_abundance.py)
```

Snakemake for the 18-sample abundance workflow:

- Workflow: `scripts/gus_abundance.smk`
- Config: `config/gus_abundance.yaml`
- Submit all samples: `sbatch scripts/submit_gus_abundance.sbatch`
- One-sample smoke test: `sbatch scripts/submit_gus_abundance_one_sample.sbatch`

Default DIAMOND / summarize filters in the Snakefile are identity ≥ 40%, query and subject coverage ≥ 50%, e-value ≤ 1e-5. Stricter profiles are applied afterwards from the same DIAMOND tables (no re-search):

| Profile | Identity | Query cov. | Subject cov. | Output directory |
| --- | --- | --- | --- | --- |
| default (Snakemake) | ≥ 40% | ≥ 50% | ≥ 50% | `results/summary/` |
| `pident90` | ≥ 90% | ≥ 50% | ≥ 50% | `results/summary/pident90/` |
| `pident95` | ≥ 95% | ≥ 50% | ≥ 50% | `results/summary/pident95/` |
| **`pident95_cov90`** | **≥ 95%** | **≥ 90%** | **≥ 90%** | **`results/summary/pident95_cov90/`** |

```bash
bash scripts/regenerate_gus_abundance_pident.sh
python3 scripts/compute_active_gus_abundance.py
python3 scripts/merge_abundance_tables.py
python3 scripts/patch_total_gus_pollet279_fmn28.py
python3 scripts/build_gus_report.py
```

Hit filtering (`scripts/summarize_gus_hits.py`):

- Keep a DIAMOND row only if `pident`, alignment / `qlen` (query coverage), alignment / `slen` (subject coverage), and e-value pass the profile.
- Keep **one best hit per predicted gene** (`qseqid`), highest bitscore.
- Query coverage = alignment length / query protein length (`qlen`).
- Subject coverage = alignment length / reference protein length (`slen`).

---

## Genome-wide TPM

Implemented in `scripts/calculate_gus_abundance.py`. The column is `tpm_genome_wide`.

This is a classic gene-level TPM, except the library is **all predicted CDS on the sample assembly**, not RNA-seq transcripts. A GUS protein’s genome-wide TPM is its share of length-normalized read counts among **every counted gene in that sample**. It is comparable across samples and is **not** renormalized inside the GUS set.

### 1. Read counts

Host-removed paired reads are mapped with Bowtie2 to that sample’s MEGAHIT contigs. `featureCounts` then counts fragments per Prodigal CDS:

```text
featureCounts -p -B -C -t CDS -g ID -a <prodigal.gff> -o <counts> <sorted.bam>
```

| Flag | Meaning |
| --- | --- |
| `-p` | Count paired-end fragments |
| `-B` | Require both ends mapped |
| `-C` | Discard chimeric fragments |
| `-t CDS -g ID` | Features are Prodigal CDS; ID matches the GFF `ID=` field |

`raw_count` for a GUS hit is the featureCounts value of the matching CDS. Prodigal FASTA headers (`k141_…`) are mapped to featureCounts gene IDs (`numeric_…`) via the `ID=` field in the `.faa` header.

`library_size` is the sum of counts over **all** CDS. It is used for CPM, not for genome-wide TPM.

### 2. Gene length

For genome-wide TPM, length is the **CDS length in base pairs** from the Prodigal GFF:

```text
cds_length_bp = end - start + 1
```

Do not use the DIAMOND query length (`qlen`, amino acids) for this metric. `qlen` is used only for GUS-set TPM (`tpm_gus_db`) and for the `rpk` column in the per-gene table.

### 3. Reads per kilobase (RPK)

For every CDS *g* with `raw_count > 0` and a GFF length:

```math
\mathrm{RPK}_{g} = \frac{\mathrm{count}_{g}}{\mathrm{cds\_length\_bp}_{g} / 1000}
```

Genes with zero counts are omitted from the sum. Genes present in counts but missing from the GFF are omitted.

### 4. Genome-wide denominator

```math
\mathrm{RPK}_{\mathrm{genome}} = \sum_{g \in \text{all CDS with count}>0} \mathrm{RPK}_{g}
```

This sum is over the **entire predicted gene catalog of the sample**, not over GUS hits. In a given sample, the sum of `tpm_genome_wide` over *all* genes would be 10<sup>6</sup>. The GUS rows are only a small slice of that million.

### 5. TPM for a GUS gene

For a predicted gene that passed the GUS hit filters:

```math
\mathrm{tpm\_genome\_wide} = \frac{\mathrm{RPK}_{\mathrm{GUS\ gene}}}{\mathrm{RPK}_{\mathrm{genome}}} \times 10^{6}
```

Worked numbers (same operations as the code):

1. Count = 6461 fragments, CDS = 2628 bp  
   RPK = 6461 / (2628 / 1000) = 2458.52
2. Suppose RPK<sub>genome</sub> = 6.63 × 10<sup>7</sup>
3. TPM = 2458.52 / 6.63×10<sup>7</sup> × 10<sup>6</sup> ≈ 37.09

That gene is about 37 of every million length-normalized gene counts in the sample.

### 6. From gene-level TPM to reference-level and clade-level tables

Several predicted genes can hit the same reference (`sseqid`). Cross-sample matrices **sum** `tpm_genome_wide` by `sseqid` within each sample:

```text
results/summary/pident95_cov90/abundance_all_subject/active_gus_abundance_gus_genome_all_subject.tsv
```

Clade tables sum those reference TPMs by clade (`results/reference/active_gus_metadata.tsv`). Sample columns can be renamed from SRR to `S` + Sample Name:

```text
.../active_gus_abundance_gus_genome_all_subject_sample_id.tsv
.../active_gus_clade_abundance_gus_genome_all_subject_sample_id.tsv
```

### 7. Percent of genome

Because TPM already scales to 10<sup>6</sup>:

```math
\text{percent of genome} = \mathrm{tpm\_genome\_wide} / 10{,}000
```

Example: TPM 37.086331 → 0.003709% of genome-wide gene abundance. Files named `*_of_genome_percent_all_subject.tsv` are this conversion.

### What genome-wide TPM is not

| Metric | Length | Denominator | Interprets as |
| --- | --- | --- | --- |
| **`tpm_genome_wide`** | CDS bp (GFF) | RPK of **all** CDS in the sample | Abundance among all genes |
| **`tpm_gus_db`** | protein `qlen` (AA) | RPK of **GUS hits only** in that table | Composition within the GUS set (sums to 10<sup>6</sup> among hits) |
| **`cpm`** | none | total fragment counts | Counts per million fragments, not length-normalized |
| **`rpk` in the TSV** | `qlen` / 1000 | none | Intermediate for `tpm_gus_db`, not for genome-wide TPM |

`tpm_gus_db` answers “of the GUS proteins we called, what fraction is this one?”. `tpm_genome_wide` answers “how abundant is this GUS gene relative to every predicted gene in the metagenome?”. Within-database percent tables (`*_within_gus_db_percent_*.tsv`) are `tpm_gus_db` rescaled so each sample column sums to 100.

If a sample has only one active GUS hit, that hit has `tpm_gus_db = 1e6` (100% of the active set) even if its genome-wide TPM is tiny.

---

## Key outputs (`pident95_cov90`)

Per sample, under `results/summary/pident95_cov90/abundance/`:

- `{SRR}.min1000.total_gus_abundance.tsv`
- `{SRR}.min1000.candidate_gus_abundance.tsv`
- `{SRR}.min1000.active_gus_abundance.tsv`

Each row is one predicted gene (best GUS hit) plus `raw_count`, `cds_length_bp`, `rpk`, `cpm`, `tpm_gus_db`, and `tpm_genome_wide`.

Merged matrices, `results/summary/pident95_cov90/abundance_all_subject/`:

`total_gus_*` matrices are the patched 307-sequence set (detected rows only). Clade is encoded in the `sseqid` suffix (Pollet `mL1mL2`, `mL1`, `mL2`, …), not a separate column.

| File | Values |
| --- | --- |
| `*_gus_genome_all_subject.tsv` | genome-wide TPM by reference × sample |
| `*_of_genome_percent_all_subject.tsv` | TPM / 10,000 |
| `*_within_gus_db_percent_all_subject.tsv` | % composition inside that GUS set |
| `*_sample_id.tsv` | same matrices with `S##` column names |

Figures and HTML/Markdown report: `results/summary/pident95_cov90/report/` (`scripts/build_gus_report.py`).

## How to run

Create the mapping environment from `environment.yml` (Python, Prodigal, DIAMOND, Bowtie2, samtools, featureCounts/subread, Snakemake):

```bash
conda env create -f environment.yml
conda activate gus-map
```

QC/host removal uses a separate `qc` environment (`scripts/submit_cleandata.sbatch`). Typical HPC order:

```bash
sbatch scripts/download_PRJNA816918_ena.sbatch          # or *_sratools.sbatch
sbatch scripts/submit_cleandata.sbatch
sbatch scripts/submit_assemble_megahit_cleandata.sbatch
sbatch scripts/submit_gus_abundance.sbatch
bash scripts/regenerate_gus_abundance_pident.sh
python3 scripts/compute_active_gus_abundance.py
python3 scripts/merge_abundance_tables.py
python3 scripts/patch_total_gus_pollet279_fmn28.py
python3 scripts/build_gus_report.py
```

Samples are not listed in `config/gus_abundance.yaml`. The abundance workflow keeps a sample if both exist:

- `cleandata/<SRR>/<SRR>_1.fastq` and `_2.fastq`
- `assembly/megahit/<SRR>/final.contigs.fa`

## Repository layout

```text
config/                     18-run accessions and gus_abundance.yaml
scripts/                    Snakemake, sbatch, and Python/R analysis
reference_from_paper/       paper sequences, Table S3, SraRunTable.csv
results/                    searches, counts, abundance, report (gitignored)
assembly/                   MEGAHIT contigs (gitignored)
cleandata/                  host-removed reads (gitignored)
```

Large inputs and results (`rawdata/`, `cleandata/`, `assembly/`, `results/`, `envs/`) are gitignored and are recreated from the workflows above.

## References:

Pollet RM, D'Agostino EH, Walton WG, Xu Y, Little MS, Biernat KA, Pellock SJ, Patterson LM, Creekmore BC, Isenberg HN, Bahethi RR, Bhatt AP, Liu J, Gharaibeh RZ, Redinbo MR. An Atlas of β-Glucuronidases in the Human Intestinal Microbiome. Structure. 2017 Jul 5;25(7):967-977.e5. doi: 10.1016/j.str.2017.05.003.

Simpson, J. B., Sekela, J. J., Graboski, A. L., Borlandelli, V. B., Bivins, M. M., Barker, N. K., … Redinbo, Matthew. R. (2022). Metagenomics combined with activity-based proteomics point to gut bacterial enzymes that reactivate mycophenolate. Gut Microbes, 14(1). https://doi.org/10.1080/19490976.2022.2107289

Guodong Zhang, Nan Jing, Guodong Cao et al. Gut Microbiota Mediate the Metabolism of Colonic Prostaglandins, 20 February 2026, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-8856024/v1]

Sarah E. Blecksmith, Karen M. Kalanetra, Cheng-Yu Weng, Christopher Suarez, Irnayuli R. Sitepu, Yirui Tang, Shawn Ehlers Cheang, Sophia Jiang, Karina Cernioglo, Karla Damian-Medina, Jennifer T. Smilowitz, Carlito B. Lebrilla, David A. Mills, Danielle G. Lemay; Fecal microbiomes from healthy adult consumers of fruits and vegetables exhibit fiber- and donor-specific fermentation: “5 a day” is not enough. Food Funct. 2025; 16 (19): 7561–7577. https://doi.org/10.1039/d5fo00947b