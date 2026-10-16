#!/usr/bin/env python3
"""Audit every detected active-GUS protein in the 17 assay subjects.
All differences are assembly-predicted, not independently validated variants.
"""
import csv,json,re,os,hashlib,subprocess,collections,itertools
from pathlib import Path
from compare_reviewer_gus_sequences import fasta
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/summary/pident95_cov90/cohort_sequence_audit'
SRC=OUT/'sources'

def rows(p):return list(csv.DictReader(open(p),delimiter='\t'))
def dump(name,data,fields=None):
 with open(OUT/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
def fa(p,records):
 with open(p,'w') as f:
  for k,s in records.items():f.write(f'>{k}\n{s}\n')
def msa(name,records,mode='--localpair'):
 fa(OUT/f'{name}.faa',records)
 with open(OUT/f'{name}.aligned.faa','w') as f,open(OUT/f'{name}.mafft.log','w') as log:
  subprocess.run(['mafft','--thread',os.environ.get('SLURM_CPUS_PER_TASK','1'),mode,'--maxiterate','1000',str(OUT/f'{name}.faa')],check=True,stdout=f,stderr=log)
 a=dict(fasta(OUT/f'{name}.aligned.faa'))
 assert len(set(map(len,a.values())))==1
 for k,s in a.items():assert s.replace('-','')==records[k]
 return a

def coord(s):return {i:c for i,c in enumerate((j for j,a in enumerate(s) if a!='-'),1)}
def compare(a,b):
 pairs=[(x,y) for x,y in zip(a,b) if x!='-' or y!='-'];a=''.join(x for x,y in pairs);b=''.join(y for x,y in pairs)
 paired=[i for i,(x,y) in enumerate(zip(a,b)) if x!='-' and y!='-'];first,last=min(paired),max(paired)
 ap=bp=0;diffs=[];subs=[];missing=[];inserted=[]
 for i,(x,y) in enumerate(zip(a,b)):
  ap+=x!='-';bp+=y!='-'
  if x!=y:
   kind='substitution' if x!='-' and y!='-' else 'reference_residue_missing' if y=='-' else 'sample_insertion'
   location='N_terminal' if i<first else 'C_terminal' if i>last else 'internal'
   if x!='-' and y=='-':missing.append((ap,bp,location))
   if x=='-':inserted.append((ap,bp,location))
   if kind=='substitution':subs.append(f'{x}{ap}{y}')
   diffs.append(dict(reference_position=ap if x!='-' else '',sample_position=bp if y!='-' else '',reference_residue=x,sample_residue=y,type=kind,location=location))
 # Near-start gaps also count as terminal-region deficits if <=10 query residues precede them;
 # this makes the terminal classification robust to ambiguous initiator-M alignment.
 blen=len(b.replace('-',''))
 nmiss=sum(1 for ap,bp,loc in missing if loc=='N_terminal' or bp<=10)
 cmiss=sum(1 for ap,bp,loc in missing if loc=='C_terminal' or blen-bp<=10)
 internal=sum(1 for ap,bp,loc in missing if loc=='internal' and bp>10 and blen-bp>10)
 match=sum(x==y for x,y in zip(a,b));plen=len(paired)
 summary=dict(reference_length=len(a.replace('-','')),protein_length=blen,identical=a.replace('-','')==b.replace('-',''),substitutions=len(subs),paired_residues=plen,identity_paired_pct=round(100*match/plen,6),identity_global_pct=round(100*match/len(a),6),missing_reference_residues=len(missing),inserted_residues=len(inserted),N_terminal_region_missing=nmiss,C_terminal_region_missing=cmiss,internal_missing_residues=internal,substitution_list=';'.join(subs))
 return summary,diffs,a,b

def main():
 OUT.mkdir(exist_ok=True,parents=True);(OUT/'pairwise').mkdir(exist_ok=True)
 phenotype=rows(ROOT/'results/summary/pident95_cov90/gus_index_digitized/subject_predictions.tsv');subjects=[r['subject'] for r in phenotype];activity={r['subject']:float(r['estimated_activity']) for r in phenotype}
 paper={h.split()[0]:s.removesuffix('*') for h,s in fasta(ROOT/'results/reference/active_gus.faa')}
 meta={r['reference_id']:r for r in rows(ROOT/'results/reference/active_gus_metadata.tsv')}
 natives={};native_source={};annotations=collections.defaultdict(list)
 for ref,acc in [('WP_012740861','C4Z6Z2'),('6NCY_A','A0A174EHD1'),('3CMG_A','Q5LIC7'),('WP_270534885','U6RLV2')]:
  d=json.load(open(SRC/f'{acc}.json'));natives[ref]=d['sequence']['value'];native_source[ref]='UniProt '+acc
  for f in d.get('features',[]):
   if f['type'] in ('Active site','Binding site'):
    annotations[ref].append((f['location']['start']['value'],f['type'],f.get('description',''),f'UniProt {acc}; '+json.dumps(f.get('evidences',[]))))
 for ref in ['WP_022501306','WP_007841259','WP_217761416']:
  seq=list(fasta(SRC/f'{ref}.faa'));assert len(seq)==1;natives[ref]=seq[0][1];native_source[ref]=seq[0][0]
 # PDB author numbering -22..0 is the 23-residue N-terminal expression extension.
 pdbrg=json.load(open(SRC/'6MVG_entity.json'))['entity_poly']['pdbx_seq_one_letter_code_can']
 assert pdbrg==paper['6MVG_A'] and pdbrg.startswith('HHHHHHSSGVDLGTENLYFQSNA')
 natives['6MVG_A']=pdbrg[23:];native_source['6MVG_A']='6MVG author-numbered residues 1..754; removed expression extension -22..0; NOT local mismatched UniProt R5TSA0'
 for pos,role in [(10,'Adjacent-protomer aglycone contact'),(377,'Substrate-recognition mutagenesis site'),(378,'YW motif'),(416,'Catalytic acid/base')]:annotations['6NCY_A'].append((pos,'Literature functional site',role,'Pellock et al 2019 https://pmc.ncbi.nlm.nih.gov/articles/PMC7977719/'))
 anc=json.load(open(SRC/'P05804.json'))['sequence']['value'];assert anc[412]=='E' and anc[503]=='E' and anc[565]=='N' and anc[567]=='K'
 alignments=[]
 for mode,name in [('--localpair','native_anchor_local'),('--globalpair','native_anchor_global')]:alignments.append(msa(name,{'ECOLI_P05804':anc,**natives},mode))
 transferred=[]
 for ref in natives:
  for ap,expected,role in [(413,'E','GH2 catalytic acid/base'),(504,'E','GH2 catalytic nucleophile'),(566,'N','GH2 NxK recognition motif'),(568,'K','GH2 NxK recognition motif')]:
   mappings=[]
   for al in alignments:
    col=coord(al['ECOLI_P05804'])[ap];inv={c:i for i,c in coord(al[ref]).items()};mappings.append(inv.get(col))
   agrees=mappings[0] is not None and mappings[0]==mappings[1];pos=mappings[0]
   residue=natives[ref][pos-1] if pos else '-'
   transferred.append(dict(reference=ref,Ecoli_position=ap,expected_residue=expected,native_position=pos or '',native_residue=residue,local_global_agree=agrees,role=role))
   if agrees and residue==expected:annotations[ref].append((pos,'Homology-inferred conserved motif',role,'Ecoli P05804 anchor; local/global MAFFT coordinate agreement; not direct assay validation'))
 dump('conserved_motif_mapping.tsv',transferred)
 runmap={'S'+str(int(r['Sample Name'])).zfill(2):r['Run'] for r in csv.DictReader(open(ROOT/'reference_from_paper/SraRunTable.csv'))}
 rpkm={(r['sample_id'],r['featurecounts_gene_id']):float(r['rpkm']) for r in rows(ROOT/'results/summary/pident95_cov90/abundance_all_subject/active_gus_rpkm_gene_details.tsv')}
 records=collections.defaultdict(dict);inventory=[]
 for s in subjects:
  run=runmap[s];hits=rows(ROOT/f'results/summary/pident95_cov90/abundance/{run}.min1000.active_gus_abundance.tsv');wanted={r['qseqid'] for r in hits}
  seqs={h.split()[0]:(h,seq) for h,seq in fasta(ROOT/f'results/prodigal/{run}.min1000.faa') if h.split()[0] in wanted};assert set(seqs)==wanted
  for hit in hits:
   ref=hit['sseqid'];key=s+'|'+hit['qseqid'];header,raw=seqs[hit['qseqid']];seq=raw.removesuffix('*');assert '*' not in seq and len(raw)==int(hit['qlen'])
   records[ref][key]=seq
   inventory.append(dict(subject=s,run=run,reference=ref,source=meta[ref]['source'],clade=meta[ref]['clade'],gene=hit['qseqid'],sequence_key=key,protein_length=len(seq),prodigal_partial=re.search(r'partial=(\d+)',header).group(1),gene_rpkm=rpkm[s,hit['featurecounts_gene_id']],estimated_activity_nM=activity[s],raw_count=hit['raw_count'],sha256=hashlib.sha256(seq.encode()).hexdigest(),prodigal_header=header))
 lookup={r['sequence_key']:r for r in inventory};gene_comparisons=[];varrows=[];constructs=[];sitechecks=[];sitetable=[];haplotypes=[]
 for ref,seqs in records.items():
  al=msa(ref,{'PAPER_REFERENCE':paper[ref],'NATIVE_REFERENCE':natives[ref],**seqs});assert ref in natives
  cs,cd,_,_=compare(al['NATIVE_REFERENCE'],al['PAPER_REFERENCE']);constructs.append(dict(reference=ref,native_source=native_source[ref],**cs))
  for d in cd:varrows.append(dict(subject='PAPER_REFERENCE',reference=ref,comparison='paper_vs_native',gene='',**d))
  for key,seq in seqs.items():
   inv=lookup[key];combined=dict(inv)
   for kind,anchor in [('paper','PAPER_REFERENCE'),('native','NATIVE_REFERENCE')]:
    sm,df,a,b=compare(al[anchor],al[key]);combined.update({kind+'_'+k:v for k,v in sm.items()})
    for d in df:varrows.append(dict(subject=inv['subject'],reference=ref,comparison=kind,gene=inv['gene'],**d))
    with open(OUT/'pairwise'/f'{key.replace("|","_")}_{kind}.txt','w') as f:
     f.write(json.dumps(sm,indent=2)+'\n')
     for i in range(0,len(a),90):f.write(f'\nREF {a[i:i+90]}\n    '+''.join('|' if x==y else ' ' for x,y in zip(a[i:i+90],b[i:i+90]))+f'\nSEQ {b[i:i+90]}\n')
   gene_comparisons.append(combined)
  ncoord=coord(al['NATIVE_REFERENCE'])
  for pos,typ,role,source in annotations[ref]:
   col=ncoord[pos];native_res=natives[ref][pos-1]
   for key in seqs:
    inv=lookup[key];invcoord={c:p for p,c in coord(al[key]).items()}
    sitechecks.append(dict(subject=inv['subject'],reference=ref,gene=inv['gene'],native_position=pos,native_residue=native_res,sample_position=invcoord.get(col,''),sample_residue=al[key][col],unchanged=al[key][col]==native_res,annotation_type=typ,role=role,source=source))
  # Full native-coordinate site matrix, only sites variable among samples (not sample/reference differences fixed in all subjects).
  for pos,col in ncoord.items():
   vals={key:al[key][col] for key in seqs}
   if len(set(vals.values()))<2:continue
   for key,res in vals.items():
    inv=lookup[key];sitetable.append(dict(subject=inv['subject'],reference=ref,gene=inv['gene'],native_position=pos,native_residue=natives[ref][pos-1],sample_residue=res,gene_rpkm=inv['gene_rpkm'],estimated_activity_nM=activity[inv['subject']]))
  groups=collections.defaultdict(list)
  for key,seq in seqs.items():groups[seq].append(key)
  for i,(seq,keys) in enumerate(sorted(groups.items(),key=lambda item:(-len(item[1]),sorted(item[1]))),1):
   hap=ref+f'_H{i:02}'
   for key in keys:
    inv=lookup[key];haplotypes.append(dict(subject=inv['subject'],reference=ref,haplotype=hap,gene=inv['gene'],group_n=len(keys),gene_rpkm=inv['gene_rpkm'],estimated_activity_nM=activity[inv['subject']],sequence_sha256=inv['sha256']))
 dump('gene_inventory.tsv',inventory);dump('gene_reference_comparisons.tsv',gene_comparisons);dump('residue_differences.tsv',varrows);dump('reference_construct_audit.tsv',constructs);dump('functional_site_checks.tsv',sitechecks);dump('variable_sites_long.tsv',sitetable);dump('protein_haplotypes.tsv',haplotypes)
 fa(OUT/'all_sample_proteins.faa',{k:s for rs in records.values() for k,s in rs.items()});fa(OUT/'native_references.faa',natives)
 # Verify exactly the abundance matrix used in the user's figure, preserving sample inclusion.
 inputmatrix=rows(ROOT/'results/summary/pident95_cov90/assay_order_stacked_rpkm/input_abundance.tsv')
 for s in subjects:
  for r in inputmatrix:
   actual=sum(g['gene_rpkm'] for g in inventory if g['subject']==s and g['reference']==r['sseqid']);assert abs(actual-float(r[s]))<1e-5
 summaries=[]
 for s in subjects:
  genes=[r for r in gene_comparisons if r['subject']==s];total=sum(r['gene_rpkm'] for r in genes)
  motif_changed={r['gene'] for r in sitechecks if r['subject']==s and not r['unchanged'] and r['annotation_type']=='Homology-inferred conserved motif'}
  summaries.append(dict(subject=s,total_RPKM=total,estimated_activity_nM=activity[s],detected_GUS_genes=len(genes),native_substitutions_sum=sum(g['native_substitutions'] for g in genes),abundance_weighted_substitutions_per100aa=sum(g['gene_rpkm']*100*g['native_substitutions']/g['native_paired_residues'] for g in genes)/total,exact_native_RPKM=sum(g['gene_rpkm'] for g in genes if g['native_identical']),exact_native_RPKM_fraction=sum(g['gene_rpkm'] for g in genes if g['native_identical'])/total,terminal_deficit_ge5_RPKM=sum(g['gene_rpkm'] for g in genes if g['native_N_terminal_region_missing']+g['native_C_terminal_region_missing']>=5),terminal_deficit_ge5_RPKM_fraction=sum(g['gene_rpkm'] for g in genes if g['native_N_terminal_region_missing']+g['native_C_terminal_region_missing']>=5)/total,prodigal_partial_genes=sum(g['prodigal_partial']!='00' for g in genes),conserved_motif_changed_genes=len(motif_changed)))
 dump('subject_sequence_summary.tsv',summaries)
 identical=[]
 for a,b in itertools.combinations(haplotypes,2):
  if a['haplotype']==b['haplotype']:
   identical.append(dict(reference=a['reference'],haplotype=a['haplotype'],subject_a=a['subject'],subject_b=b['subject'],RPKM_a=a['gene_rpkm'],RPKM_b=b['gene_rpkm'],activity_a=a['estimated_activity_nM'],activity_b=b['estimated_activity_nM'],activity_fold=max(a['estimated_activity_nM'],b['estimated_activity_nM'])/min(a['estimated_activity_nM'],b['estimated_activity_nM'])))
 dump('identical_protein_subject_pairs.tsv',identical)
 print(f'Finished: {len(subjects)} subjects, {len(inventory)} predicted proteins, {len(records)} detected reference groups. All RPKM cells verified.')
 print('Construct audit:',json.dumps(constructs,indent=2));print('Changed functional sites:',json.dumps([r for r in sitechecks if not r['unchanged']],indent=2))
 (OUT/'methods.txt').write_text('All 17 assay subjects; exact same RPKM matrix as supplied. Existing active-GUS pident95_cov90 best-hit calls only. S02 excluded. All detected CDS retained, terminal stop removed; no other sample sequence editing. MAFFT 7.505 L-INS-i separately per reference including supplied paper construct and a native/construct-core reference. Native-source choices in reference_construct_audit.tsv. Pairwise metrics exclude double-gap columns from the group MSA. Differences to paper and native kept separately; counts of substituted positions and gap residues are descriptive, not functional scores. Near-N-terminal deficits include reference gaps with <=10 sample residues preceding them; near-C-terminal analogous. These are alignment features, not confirmed genetic deletions. Native annotated sites are supplemented by homology-inferred Ecoli P05804 E413/E504/N566/K568 only when local and global native-reference MSAs agree and the expected residue is present in the target reference. Absence of a changed annotated site does not exclude functional effects elsewhere. Sample proteins remain assembly-derived; no new read-level validation. Enzymatic activities are approximate 15-min figure-digitized values, not raw replicate data. Source files and exact hashes retained.\n')
if __name__=='__main__':main()
