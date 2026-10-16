#!/usr/bin/env python3
import urllib.request,json,csv,hashlib,shutil
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root/'results/summary/pident95_cov90/cohort_sequence_audit/sources'
out.mkdir(parents=True,exist_ok=True)
urls={a+'.json':'https://rest.uniprot.org/uniprotkb/'+a+'.json' for a in ['P05804','R5TSA0','C4Z6Z2','Q5LIC7','A0A174EHD1','U6RLV2']}
for a in ['WP_022501306','WP_270534885','WP_007841259','WP_217761416']:
 urls[a+'.faa']='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id='+a+'&rettype=fasta&retmode=text'
manifest=[]
for name,url in urls.items():
 p=out/name
 if not p.exists():
  prev=root/'results/summary/pident95_cov90/reviewer_site_conservation/sources'/name
  if prev.exists():shutil.copyfile(prev,p)
  else:
   try:p.write_bytes(urllib.request.urlopen(url,timeout=40).read())
   except Exception as e:print(name,e);continue
 b=p.read_bytes()
 if name.endswith('.json'):
  d=json.loads(b)
  if 'sequence' not in d:
   print(name,'No sequence in record:',d);continue
  print(name,d.get('organism',{}).get('scientificName'),len(d['sequence']['value']), [(f['type'],f['location']['start'].get('value'),f['location']['end'].get('value'),f.get('description')) for f in d.get('features',[]) if f['type'] in ('Active site','Binding site','Signal','Chain')])
 else:print(name,b[:150].decode(),sum(len(s) for s in b.decode().splitlines() if not s.startswith('>')))
 manifest.append(dict(file=name,url=url,sha256=hashlib.sha256(b).hexdigest()))
with open(out/'source_manifest.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=['file','url','sha256'],delimiter='\t');w.writeheader();w.writerows(manifest)
