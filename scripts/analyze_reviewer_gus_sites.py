#!/usr/bin/env python3
"""Exploratory site frequencies in cohort and Pollet279 homolog panels; not selection scores."""
import csv,json,os,subprocess,collections,hashlib
from pathlib import Path
from compare_reviewer_gus_sequences import fasta
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/summary/pident95_cov90/reviewer_site_conservation'
PREV=ROOT/'results/summary/pident95_cov90/reviewer_sequence_comparison'
TARGETS={'WP_012740861':'C4Z6Z2','6NCY_A':'A0A174EHD1'}

def dump(name,rows):
 with open(OUT/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fa(path,records):
 with open(path,'w') as f:
  for k,s in records.items():f.write(f'>{k}\n{s}\n')
def run(cmd,stdout=None,stderr=None):
 subprocess.run(list(map(str,cmd)),check=True,stdout=stdout,stderr=stderr)
def main():
 OUT.mkdir(exist_ok=True,parents=True)
 threads=os.environ.get('SLURM_CPUS_PER_TASK','1')
 refs={k:json.load(open(OUT/'sources'/f'{acc}.json'))['sequence']['value'] for k,acc in TARGETS.items()}
 fa(OUT/'native_queries.faa',refs)
 diamond=ROOT/'envs/gus-map/bin/diamond'
 run([diamond,'blastp','--query',OUT/'native_queries.faa','--db',ROOT/'results/reference/gus_279.dmnd','--out',OUT/'pollet279_hits.tsv','--threads',threads,'--sensitive','--max-target-seqs','0','--evalue','1e-5','--outfmt','6','qseqid','sseqid','pident','qcovhsp','scovhsp','bitscore'])
 candidates=collections.defaultdict(dict)
 with open(OUT/'pollet279_hits.tsv') as f:
  for q,s,identity,qcov,scov,score in csv.reader(f,delimiter='\t'):
   if float(identity)>=40 and min(float(qcov),float(scov))>=80:candidates[q][s]=float(identity)
 db={h.split()[0]:s.removesuffix('*') for h,s in fasta(ROOT/'results/reference/gus_279.faa')}
 smap={r['Run']:f"S{int(r['Sample Name']):02}" for r in csv.DictReader(open(ROOT/'reference_from_paper/SraRunTable.csv'))}
 cohort=collections.defaultdict(dict);inventory=[]
 for runid,sample in smap.items():
  rows=[r for r in csv.DictReader(open(ROOT/f'results/summary/pident95_cov90/abundance/{runid}.min1000.active_gus_abundance.tsv'),delimiter='\t') if r['sseqid'] in TARGETS]
  if not rows:continue
  wanted={r['qseqid'] for r in rows}
  seqs={h.split()[0]:s.removesuffix('*') for h,s in fasta(ROOT/f'results/prodigal/{runid}.min1000.faa') if h.split()[0] in wanted}
  assert set(seqs)==wanted
  for r in rows:
   k=sample+'|'+r['qseqid'];cohort[r['sseqid']][k]=seqs[r['qseqid']]
   inventory.append(dict(reference=r['sseqid'],sample=sample,run=runid,gene=r['qseqid'],length=len(seqs[r['qseqid']]),tpm=r['tpm_genome_wide']))
 dump('cohort_inventory.tsv',inventory)
 variants=collections.defaultdict(set)
 for r in csv.DictReader(open(PREV/'residue_differences.tsv'),delimiter='\t'):
  if r['comparison']!='sample':continue
  if r['sequence_a'].startswith('S11|') and r['residue_a']!='-' and r['residue_b']!='-':variants[r['reference']].add(int(r['position_a']))
  if r['sequence_b'].startswith('S11|') and r['residue_a']!='-' and r['residue_b']!='-':variants[r['reference']].add(int(r['position_b']))
 stats=[];sites=[];functional=[];panels=[]
 for ref,acc in TARGETS.items():
  records={'NATIVE':refs[ref],**cohort[ref]};seen=set();ids={}
  for k,ident in sorted(candidates[ref].items()):
   if db[k] in seen:continue
   seen.add(db[k]);key='DB|'+k;records[key]=db[k];ids[key]=ident
   panels.append(dict(reference=ref,sequence_id=k,identity_to_native=ident,length=len(db[k])))
  fa(OUT/f'{ref}.unaligned.faa',records)
  with open(OUT/f'{ref}.aligned.faa','w') as out,open(OUT/f'{ref}.mafft.log','w') as log:
   run(['mafft','--thread',threads,'--localpair','--maxiterate','1000',OUT/f'{ref}.unaligned.faa'],out,log)
  aln={h.split()[0]:s for h,s in fasta(OUT/f'{ref}.aligned.faa')}
  assert len(set(map(len,aln.values())))==1
  for k,s in aln.items():assert s.replace('-','')==records[k]
  s11=next(k for k in aln if k.startswith('S11|'))
  def coord(k):return {i:col for i,col in enumerate((c for c,a in enumerate(aln[k]) if a!='-'),1)}
  c11=coord(s11);cn=coord('NATIVE');col_to_native={v:k for k,v in cn.items()}
  for pos in sorted(variants[ref]):
   col=c11[pos];sites.append(dict(reference=ref,s11_position=pos,s11_residue=aln[s11][col],native_position=col_to_native.get(col,''),native_residue=aln['NATIVE'][col],four_sample_residues=';'.join(k.split('|')[0]+':'+aln[k][col] for k in cohort[ref] if k.split('|')[0] in ('S10','S11','S14','S21'))))
   groups={'cohort_all_detected_genes':list(cohort[ref])}
   uniq={}
   for k in cohort[ref]:uniq.setdefault(records[k],k)
   groups['cohort_unique_proteins']=list(uniq.values())
   for threshold in (40,60,80):groups[f'Pollet279_unique_identity_ge{threshold}']=[k for k,v in ids.items() if v>=threshold]
   for label,keys in groups.items():
    counts=collections.Counter(aln[k][col] for k in keys);valid=sum(v for a,v in counts.items() if a not in '-X');ordered=sorted(((n,a) for a,n in counts.items() if a not in '-X'),reverse=True)
    stats.append(dict(reference=ref,s11_position=pos,native_position=col_to_native.get(col,''),s11_residue=aln[s11][col],panel=label,total_sequences=len(keys),nongap_nonX=valid,residue_counts=';'.join(f'{a}:{n}' for a,n in sorted(counts.items())),consensus_residue=ordered[0][1] if ordered else '',consensus_fraction=round(ordered[0][0]/valid,4) if valid else '',s11_residue_fraction=round(counts[aln[s11][col]]/valid,4) if valid else ''))
  # Structure-derived ligand contacts from UniProt; Fs sites explicitly described in Pellock 2019.
  annotations=[]
  d=json.load(open(OUT/'sources'/f'{acc}.json'))
  for f in d.get('features',[]):
   if f['type'] in ('Active site','Binding site'):
    annotations.append((f['location']['start']['value'],f['type'],'UniProt '+acc+' '+json.dumps(f.get('evidences',[]))))
  if ref=='6NCY_A':annotations += [(10,'Adjacent-protomer aglycone contact','Pellock 2019 Fig4'),(377,'Substrate recognition, mutagenesis evidence','Pellock 2019 Fig4'),(378,'YW motif','Pellock 2019 Fig4'),(416,'Catalytic acid/base','Pellock 2019 text')]
  for pos,role,source in annotations:
   col=cn[pos]
   for k in cohort[ref]:
    if k.split('|')[0] in ('S10','S11','S14','S21'):
     functional.append(dict(reference=ref,native_position=pos,native_residue=refs[ref][pos-1],role=role,source=source,sample=k.split('|')[0],sample_residue=aln[k][col],unchanged=aln[k][col]==refs[ref][pos-1]))
 dump('site_coordinate_mapping.tsv',sites);dump('site_frequencies.tsv',stats);dump('functional_sites_four_samples.tsv',functional);dump('pollet_homolog_panel.tsv',panels)
 report=['Exploratory conservation analysis', 'MAFFT 7.505 L-INS-i; native UniProt anchors; cohort includes all 18 metagenomes (including S02 if homolog detected); original pident95_cov90 calls. Pollet279 homologs selected by DIAMOND >=40% identity and >=80% coverage on both sequences, E<=1e-5; exact sequence duplicates removed within reference panel. Subpanels >=60/80% use the same MSA. Reported residue frequencies exclude gaps/X and have no phylogenetic weighting. These are descriptive frequencies in selected panels, not calibrated evolutionary conservation scores; small panels cannot establish broad conservation. Functional annotations are limited and absence from this list does not imply absence of functional importance. S11 and native positions explicitly mapped. No structural distance or functional effect prediction performed.']
 for ref in TARGETS:
  report += [f'{ref}: cohort genes={len(cohort[ref])}, Pollet unique homologs={sum(p["reference"]==ref for p in panels)}']
 report+=['Site mapping:',json.dumps(sites,indent=2),'No observed substitution overlaps the retrieved annotated sites: '+str(not any(not f['unchanged'] for f in functional))]
 (OUT/'methods_and_summary.txt').write_text('\n'.join(report)+'\n');print('\n'.join(report))
if __name__=='__main__':main()
