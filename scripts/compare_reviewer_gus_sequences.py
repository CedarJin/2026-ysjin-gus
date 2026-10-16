#!/usr/bin/env python3
"""Compare assembled GUS proteins; global unit-cost edit alignment, no dependencies."""
import csv, hashlib, itertools, json, re
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/summary/pident95_cov90/reviewer_sequence_comparison'

def fasta(path):
    header, parts = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None: yield header, ''.join(parts)
                header, parts = line[1:].strip(), []
            else: parts.append(line.strip())
    if header is not None: yield header, ''.join(parts)

def align(a,b):
    # Needleman-Wunsch / Levenshtein: match 0, substitution/insertion/deletion 1.
    # Track number of optimal paths (capped at 2) to flag ambiguous coordinates.
    n,m=len(a),len(b)
    trace=[bytearray(m+1) for _ in range(n+1)]
    prev=list(range(m+1)); ways=[1]*(m+1)
    for j in range(1,m+1): trace[0][j]=2
    for i in range(1,n+1):
        curr=[i]+[0]*m; cw=[1]+[0]*m; trace[i][0]=1
        for j in range(1,m+1):
            scores=(prev[j-1]+(a[i-1]!=b[j-1]),prev[j]+1,curr[j-1]+1)
            best=min(scores); curr[j]=best; trace[i][j]=scores.index(best)
            cw[j]=min(2,sum(w for s,w in zip(scores,(ways[j-1],ways[j],cw[j-1])) if s==best))
        prev,ways=curr,cw
    aa=[];bb=[];i=n;j=m
    while i or j:
        t=trace[i][j]
        if t==0: i-=1;j-=1;aa.append(a[i]);bb.append(b[j])
        elif t==1: i-=1;aa.append(a[i]);bb.append('-')
        else: j-=1;aa.append('-');bb.append(b[j])
    return ''.join(reversed(aa)),''.join(reversed(bb)),prev[m],ways[m]>1

def write_tsv(name, rows):
    with open(OUT/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    samples={f"S{int(r['Sample Name']):02}":r['Run'] for r in csv.DictReader(open(ROOT/'reference_from_paper/SraRunTable.csv')) if int(r['Sample Name']) in (10,11,14,21)}
    records=[];manifest=[]
    for sample,run in sorted(samples.items()):
        path=ROOT/f'results/summary/pident95_cov90/abundance/{run}.min1000.active_gus_abundance.tsv'
        rows=list(csv.DictReader(open(path),delimiter='\t'))
        wanted={r['qseqid'] for r in rows}
        seqs={h.split()[0]:(h,s) for h,s in fasta(ROOT/f'results/prodigal/{run}.min1000.faa') if h.split()[0] in wanted}
        assert set(seqs)==wanted
        for r in rows:
            h,raw=seqs[r['qseqid']];seq=raw.removesuffix('*')
            assert len(raw)==int(r['qlen']); assert '*' not in seq
            key=sample+'|'+r['qseqid']
            records.append(dict(key=key,sample=sample,run=run,reference=r['sseqid'],seq=seq))
            manifest.append(dict(sample=sample,run=run,reference=r['sseqid'],clade=r['reference_clade'],gene=r['qseqid'],protein_length=len(seq),terminal_stop_removed=raw.endswith('*'),prodigal_partial=re.search(r'partial=(\d+)',h).group(1),tpm=r['tpm_genome_wide'],raw_fragment_count=r['raw_count'],diamond_pident=r['pident'],sha256=hashlib.sha256(seq.encode()).hexdigest(),prodigal_header=h))
    refs={h.split()[0]:s.removesuffix('*') for h,s in fasta(ROOT/'results/reference/active_gus.faa')}
    with open(OUT/'sample_proteins.faa','w') as f:
        for r in records:f.write('>'+r['key']+' reference='+r['reference']+'\n'+r['seq']+'\n')
    write_tsv('sequence_inventory.tsv',manifest)
    pairs=[];variants=[];refrows=[]
    (OUT/'alignments').mkdir(exist_ok=True)
    comparisons=[]
    for ref in sorted({r['reference'] for r in records}):
        group=[r for r in records if r['reference']==ref]
        rr=dict(key='REF|'+ref,sample='REF',reference=ref,seq=refs[ref])
        comparisons.extend(('reference',rr,r) for r in group)
        comparisons.extend(('sample',a,b) for a,b in itertools.combinations(group,2))
    for kind,a,b in comparisons:
        x,y,d,amb=align(a['seq'],b['seq']); assert x.replace('-','')==a['seq'] and y.replace('-','')==b['seq']
        matches=sum(u==v for u,v in zip(x,y));sub=sum(u!=v and u!='-' and v!='-' for u,v in zip(x,y));gaps=sum(u=='-' or v=='-' for u,v in zip(x,y))
        assert sub+gaps==d
        row=dict(reference=a['reference'],sequence_a=a['key'],sequence_b=b['key'],length_a=len(a['seq']),length_b=len(b['seq']),identical=a['seq']==b['seq'],matches=matches,substitutions=sub,gap_residues=gaps,alignment_columns=len(x),identity_all_columns_pct=round(matches/len(x)*100,6),identity_paired_residues_pct=round(matches/(len(x)-gaps)*100,6),edit_distance=d,multiple_optimal_alignments=amb)
        (pairs if kind=='sample' else refrows).append(row)
        ap=bp=0;changes=[]
        for col,(u,v) in enumerate(zip(x,y),1):
            ap+=u!='-';bp+=v!='-'
            if u!=v:
                variants.append(dict(comparison=kind,reference=a['reference'],sequence_a=a['key'],sequence_b=b['key'],alignment_column=col,position_a=ap if u!='-' else '',position_b=bp if v!='-' else '',residue_a=u,residue_b=v,multiple_optimal_alignments=amb))
                changes.append(f'{u}{ap if u!="-" else "after"+str(ap)}{v}')
        name=(a['key']+'__'+b['key']).replace('|','_')
        with open(OUT/'alignments'/f'{name}.txt','w') as f:
            f.write(json.dumps(row,indent=2)+'\nDifferences (A numbering): '+', '.join(changes)+'\n')
            for k in range(0,len(x),80):f.write(f'\nA {x[k:k+80]}\n  '+''.join('|' if u==v else ' ' for u,v in zip(x[k:k+80],y[k:k+80]))+f'\nB {y[k:k+80]}\n')
    write_tsv('pairwise_sample_comparisons.tsv',pairs);write_tsv('reference_comparisons.tsv',refrows);write_tsv('residue_differences.tsv',variants)
    report=['Comparison of S10, S11, S14 and S21 assembled GUS proteins', '', 'Methods: existing pident95_cov90 active-GUS best-hit calls; all retained genes included. Sample mapping from SraRunTable.csv. Prodigal terminal stop symbols removed; no amino acids trimmed. Exact string equality plus global minimum-edit alignments (match 0, mismatch/gap residue 1; diagonal-first tie break). Identity = matches / all alignment columns. Ambiguous optimal alignments flagged. Reference coordinates refer to the supplied FASTA, not verified biological/structure numbering. Gene abundance is DNA-based genome-wide TPM. No reads were remapped and sequence variants are assembly-derived, not independently validated.', '', 'Inventory:']
    report += [f"{r['sample']} {r['reference']} length={r['protein_length']} partial={r['prodigal_partial']} TPM={r['tpm']}" for r in manifest]
    report += ['', 'S11 comparisons:']
    report += [json.dumps(r) for r in pairs if 'S11|' in r['sequence_a'] or 'S11|' in r['sequence_b']]
    report += ['', 'Interpretation limits: sequence differences do not establish altered catalytic activity. Terminal reference differences may reflect constructs or annotation and are not assumed to be biological deletions. Prodigal partial=00 is a prediction, not experimental confirmation of full length. No functional-site annotation or causal attribution is made.']
    (OUT/'report.txt').write_text('\n'.join(report)+'\n')
    print('\n'.join(report))
if __name__=='__main__':main()
