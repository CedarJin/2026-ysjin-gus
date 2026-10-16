#!/usr/bin/env python3
"""Rebuild the active GUS sample-ID matrix as mapped-read RPKM.

Run on a compute node. Uses nucleotide CDS lengths and read counts, not
protein lengths or rounded TPM. Only the requested matrix is replaced.
"""
import csv
import io
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/summary/pident95_cov90/abundance_all_subject'
TARGET = OUT / 'active_gus_abundance_gus_genome_all_subject_sample_id.tsv'


def read_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    old = read_rows(TARGET)
    columns = list(old[0])
    references = [r['sseqid'] for r in old]
    assert len(references) == len(set(references))
    runs = read_rows(ROOT / 'config/PRJNA816918_runs.tsv')
    sample_runs = {'S' + r['library_name'].split('_')[0].zfill(2): r['run_accession'] for r in runs}
    assert len(sample_runs) == len(runs)
    values = {ref: {} for ref in references}
    audit = []
    details = []
    for sample in columns[1:]:
        run = sample_runs[sample]
        counts_path = ROOT / f'results/counts/{run}.min1000.featurecounts.txt'
        with counts_path.open() as f:
            command = f.readline()
            assert 'v2.1.1' in command and '"-p"' in command and '--countReadPairs' not in command
            counts = {r['Geneid']: (int(r['Length']), int(r[next(reversed(r))])) for r in csv.DictReader(f, delimiter='\t')}
        summary = read_rows(Path(str(counts_path) + '.summary'))
        totals = {r['Status']: int(r[next(reversed(r))]) for r in summary}
        total = sum(totals.values())
        mapped = total - totals['Unassigned_Unmapped']
        assert mapped > 0 and sum(c for _, c in counts.values()) == totals['Assigned']
        log = (ROOT / f'results/mapping/{run}.min1000.bowtie2.log').read_text()
        pairs = int(re.search(r'(?m)^(\d+) reads; of these:', log)[1])
        unmapped = int(re.search(r'(?m)^        (\d+) \([^\n]+\) aligned 0 times$', log)[1])
        summary_total = total
        summary_mapped = mapped
        total = 2 * pairs
        mapped = total - unmapped
        assert 0 < mapped <= total and totals['Assigned'] <= mapped
        hits = read_rows(ROOT / f'results/summary/pident95_cov90/abundance/{run}.min1000.active_gus_abundance.tsv')
        seen = set()
        aggregated = dict.fromkeys(references, 0.0)
        tpm = dict.fromkeys(references, 0.0)
        for row in hits:
            gene = row['featurecounts_gene_id']
            assert gene not in seen
            seen.add(gene)
            length, count = counts[gene]
            assert length == int(row['cds_length_bp']) and count == int(row['raw_count'])
            assert int(row['library_size']) == totals['Assigned']
            ref = row['sseqid']
            assert ref in aggregated and length > 0
            value = count * 1e9 / (length * mapped)
            assert math.isfinite(value) and value >= 0
            aggregated[ref] += value
            tpm[ref] += float(row['tpm_genome_wide'])
            details.append([sample, run, ref, gene, count, length, mapped, f'{value:.12f}'])
        # Validate sample mapping against the original TPM table on first run.
        backup = TARGET.with_suffix('.TPM_before_RPKM.tsv')
        original = read_rows(backup) if backup.exists() else old
        for row in original:
            assert abs(float(row[sample]) - tpm[row['sseqid']]) <= 1e-5
        for ref in references:
            values[ref][sample] = aggregated[ref]
        audit.append([sample, run, total, unmapped, mapped, totals['Assigned'], len(hits), summary_total, summary_mapped, total - summary_total])
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter='\t', lineterminator='\n')
    writer.writerow(columns)
    for ref in references:
        writer.writerow([ref] + [f'{values[ref][s]:.6f}' for s in columns[1:]])
    backup = TARGET.with_suffix('.TPM_before_RPKM.tsv')
    if not backup.exists():
        backup.write_bytes(TARGET.read_bytes())
    for name, header, rows in [
        ('active_gus_rpkm_library_sizes.tsv', ['sample_id','run_accession','total_reads','unmapped_reads','mapped_reads','assigned_reads','active_gus_CDS_count','featurecounts_summary_total','featurecounts_summary_mapped','bowtie2_minus_summary_reads'], audit),
        ('active_gus_rpkm_gene_details.tsv', ['sample_id','run_accession','sseqid','featurecounts_gene_id','raw_count','cds_length_bp','mapped_reads','rpkm'], details),
    ]:
        with (OUT / name).open('w') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(header)
            w.writerows(rows)
    temporary = TARGET.with_suffix('.tmp')
    temporary.write_text(buffer.getvalue())
    temporary.replace(TARGET)
    print(f'Updated {TARGET}: {len(references)} references x {len(columns)-1} samples; {len(details)} CDS records validated.')
    print(f'Original TPM backup: {backup}')


if __name__ == '__main__':
    main()
