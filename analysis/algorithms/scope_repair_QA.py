"""Independent original-source reconstruction of all amended condition categories."""
from pathlib import Path
import collections,csv,datetime,gzip,hashlib,json,re,unicodedata
R=Path(__file__).resolve().parents[1];O=R/'stage42'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def normalize(s):return re.sub(r'\s+',' ',unicodedata.normalize('NFKC',s)).strip().casefold()
p=O/'logs/protocol.json';pr=load(p);assert sha(p)==(O/'logs/protocol_sha256.txt').read_text().strip();old=load(R/'stage39/logs/protocol.json');ext=load(old['external_dictionary']['path']);ps=load(old['external_panel']['path']);green={g['symbol']:g for g in ps['genes'] if str(g['confidence_level'])=='3'};assert len(green)==35
sources={k:Path(v['path']) for k,v in ext['sources'].items()}
for v in ext['sources'].values():assert sha(v['path'])==v['sha256']
terms={};term=None
with sources['hp.obo'].open(encoding='utf-8') as f:
 for line in f:
  s=line.rstrip('\n')
  if s=='[Term]':term={'parents':set(),'alts':set(),'obsolete':False}
  elif s.startswith('['):term=None
  elif term is not None:
   if s.startswith('id: '):term['id']=s[4:];terms[term['id']]=term
   elif s.startswith('name: '):term['name']=s[6:]
   elif s.startswith('is_a: '):term['parents'].add(s[6:].split(' ! ')[0].strip())
   elif s.startswith('alt_id: '):term['alts'].add(s[8:])
   elif s=='is_obsolete: true':term['obsolete']=True
closure={'HP:0000787','HP:0000121'}
while True:
 previous=len(closure);closure|={k for k,t in terms.items() if not t['obsolete'] and t['parents']&closure}
 if len(closure)==previous:break
alts={a:k for k in closure for a in terms[k]['alts']};target_names={normalize(terms[k]['name']) for k in closure}
def table(path,striphash=False):
 with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader((s.lstrip('#') if striphash else s for s in f if s.strip()),delimiter='\t'))
official=table(sources['NCBI_gene_condition_source_id.tsv'],True);gd=table(sources['genes_to_disease.txt']);identities={};names={};allowed={};universe_disease_ids=set()
for gene,g in green.items():
 mim={x for s in g['phenotypes'] for x in re.findall(r'(?<![0-9])([1-6][0-9]{5})(?![0-9])',s)};aa={'OMIM:'+x for x in mim}|{r['disease_id'] for r in gd if r['gene_symbol']==gene and r['association_type']=='MENDELIAN'};im=collections.defaultdict(set);nm=collections.defaultdict(set)
 for r in official:
  if gene not in {x.strip() for x in re.split('[,;|]',r['AssociatedGenes'])}:continue
  ds={'OMIM:'+x for x in re.findall(r'\b[1-6][0-9]{5}\b',r['DiseaseMIM'])};cid=r['ConceptID'].upper()
  if re.fullmatch(r'C[0-9]{7}|CN[0-9]{6}',cid):ds.add('MedGen:'+cid)
  im[cid]|=ds;nm[normalize(r['DiseaseName'])]|=ds;aa|=ds
 identities[gene]=im;names[gene]=nm;allowed[gene]=aa;universe_disease_ids|={'OMIM:'+x for x in mim}|{r['disease_id'] for r in gd if r['gene_symbol']==gene}
annotations=collections.defaultdict(list)
with sources['phenotype.hpoa'].open(encoding='utf-8',newline='') as f:
 for r in csv.DictReader((s for s in f if not s.startswith('#')),delimiter='\t'):
  if r['hpo_id'] in closure and r['database_id'] in universe_disease_ids:annotations[r['database_id']].append(r)
assert sum(map(len,annotations.values()))==len(ext['preserved_raw_target_phenotype_annotation_rows'])
def disease(ids):
 found=[r for d in ids for r in annotations.get(d,[])];positive=negative=False
 for r in found:
  f=r['frequency'].strip();no=r['qualifier'].strip().upper()=='NOT' or f in ('HP:0040285','0','0%') or re.fullmatch(r'0/[0-9]+',f) is not None
  if no:negative=True
  else:positive=True
 return 'mixed_positive_negative' if positive and negative else 'gene_disease_positive_stone_NC_annotation' if positive else 'mapped_named_disease_without_positive_annotation'
def ids(text,mode):
 medpattern=r'\bC[0-9]{7}\b' if mode=='A' else r'\b(?:CN[0-9]{6}|C[0-9]{7})\b'
 return ({m.group().upper() for m in re.finditer(medpattern,text,re.I)},{m.group().upper() for m in re.finditer(r'HP:[0-9]{7}',text,re.I)},{'OMIM:'+m.group(1) for m in re.finditer(r'OMIM[#: ]*([0-9]{6})',text,re.I)})
def permitted(ds,mode):return {x for x in ds if not (mode=='A' and x.startswith('MedGen:CN'))}
def compute(s,g,mode,use_names):
 tags=[];components=[]
 for c in [x.strip() for x in re.split('[|;]',s['ReportedPhenotypeInfo']) if x.strip()] or ['']:
  m,h,o=ids(c,mode);known=o&allowed[g]
  for x in m:known|=permitted(identities[g].get(x,set()),mode)
  prefix=r'^\s*(?:C[0-9]{7}|na)\s*:\s*' if mode=='A' else r'^\s*(?:C[0-9]{7}|CN[0-9]{6}|na)\s*:\s*';name=normalize(re.sub(prefix,'',c,flags=re.I))
  direct=bool(m&{'C0392525','C0027709'} or h&closure or mode=='B' and h&set(alts) or use_names and name in target_names)
  if direct:tag='direct_stone_NC_condition'
  else:
   if use_names:known|=permitted(names[g].get(name,set()),mode)
   tag=disease(known) if known else 'unknown_unresolved_condition'
  tags.append(tag);components.append((c,tag,known,m,h,o))
 m,h,o=ids(s['SubmittedPhenotypeInfo'],mode);known=o&allowed[g];supp=[]
 for x in m:known|=permitted(identities[g].get(x,set()),mode)
 if m&{'C0392525','C0027709'} or h&closure or mode=='B' and h&set(alts):supp.append('direct_stone_NC_condition')
 if known:supp.append(disease(known))
 tags+=supp
 return tags[0] if len(set(tags))==1 else 'mixed_scope_or_unresolved_components',components,supp,known
result=load(O/'logs/repair_results_pending_independent_QA.json');assert sha(result['output']['path'])==result['output']['sha256'];n_variants=n_scvs=n_components=0;allcounts=collections.defaultdict(collections.Counter)
with gzip.open(result['output']['path'],'rt',encoding='utf-8') as f:
 for line in f:
  v=json.loads(line);g=v['gene'];vv=collections.defaultdict(list)
  for s in v['all_source_contributing_SCV_evidence']:
   for mode in ('A','B'):
    for label,use_names in [('with_official_names',True),('identifiers_only',False)]:
     key=mode+'_'+label;tag,cc,supp,known=compute(s,g,mode,use_names);expected=s[key];assert tag==expected['scope'] and supp==expected['supplementary_explicit_ID_scope'] and known==set(expected['supplementary_gene_associated_disease_ids']),(v['VariationID'],s['SCV'],key)
     assert len(cc)==len(expected['component_evidence'])
     for (raw,t,ds,med,hp,om),e in zip(cc,expected['component_evidence']):
      assert (raw,t,ds,med,hp,om)==(e['raw_component'],e['tag'],set(e['gene_associated_disease_ids']),set(e['MedGen_ids']),set(e['HPO_ids']),set(e['OMIM_ids']));n_components+=1
     if s['tier']=='PLP':vv[key].append(tag)
   n_scvs+=1
  fields={'A_with_official_names':'A_corrected_official_names_scope','B_with_official_names':'B_corrected_source_supported_scope','A_identifiers_only':'A_identifiers_only_scope','B_identifiers_only':'B_identifiers_only_scope'}
  for key,field in fields.items():
   tags=set(vv[key]);aggregate=next(iter(tags)) if len(tags)==1 else 'mixed_scope_across_PLP_submissions' if tags else 'no_resolved_PLP_contributing_SCV';assert aggregate==v[field];allcounts[field][aggregate]+=1
  n_variants+=1
assert n_variants==4078 and n_scvs==6831
assert dict(allcounts['A_corrected_official_names_scope'])==result['A_remove_sole_panel_shortcut_counts'] and dict(allcounts['B_corrected_source_supported_scope'])==result['B_corrected_source_supported_counts'] and dict(allcounts['B_identifiers_only_scope'])==result['B_identifiers_only_counts']
record={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS_ALL4078_CORRECTED_SCOPE_ALL6831_CONTRIBUTORS_FROM_ORIGINAL_SOURCES','protocol_sha256':sha(p),'script_sha256':sha(__file__),'all4078_variants_checked':n_variants,'all6831_contributor_records_checked':n_scvs,'all_four_mode_component_checks':n_components,'original_HPO_closure_and_alt_IDs_rebuilt':True,'original_official_AssociatedGenes_and_gene_disease_relations_rebuilt':True,'every_scope_component_identifier_set_and_supplementary_evidence_checked':True,'all_category_counts':{k:dict(v) for k,v in allcounts.items()},'scope':'Independent implementation reconstruction on the same current source dictionaries, not independent clinical/functional validity. No automatic legacy-history merge; absent annotations remain unknown.'}
dump(O/'logs/independent_original_dictionary_scope_repair_QA.json',record);print('STAGE42_INDEPENDENT_ORIGINAL_DICTIONARY_REPAIR_QA_PASS',n_variants,n_scvs,n_components,flush=True)
