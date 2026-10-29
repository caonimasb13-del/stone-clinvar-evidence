"""Post-review scope correction; all original4078variants retained."""
from pathlib import Path
import collections,csv,datetime,gzip,hashlib,json,re,unicodedata
R=Path(__file__).resolve().parents[1];O=R/'stage42'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def norm(s):return ' '.join(unicodedata.normalize('NFKC',s).split()).casefold()
p=O/'logs/protocol.json';pr=json.loads(p.read_text(encoding='utf-8'));assert sha(p)==(O/'logs/protocol_sha256.txt').read_text().strip()
for x in pr['inputs']:assert sha(x['path'])==x['sha256'] and Path(x['path']).stat().st_size==x['bytes']
source=next(Path(x['path']) for x in pr['inputs'] if Path(x['path']).name=='all_latest_PLP_aggregate_condition_scope_and_sources.jsonl.gz')
ep=next(Path(x['path']) for x in pr['inputs'] if Path(x['path']).name=='external_fixed35_gene_condition_phenotype_identity_only.json');ext=json.loads(ep.read_text(encoding='utf-8'))
for s in ext['sources'].values():assert sha(s['path'])==s['sha256']
D=source.parent.parent/'derived_stage42';D.mkdir(exist_ok=True)
lookup={};names={};gene_ids={};old_panel={}
for g in ext['all35_external_gene_records']:
 lm=collections.defaultdict(set);ln=collections.defaultdict(set);gd=set(g['HPO_Mendelian_gene_disease_ids'])|{'OMIM:'+i for i in g['panel_explicit_OMIM_disease_ids']}
 for a in g['official_associated_gene_condition_rows']:
  ids={'OMIM:'+s for s in re.findall(r'\b[1-6]\d{5}\b',a['DiseaseMIM'])};cid=a['ConceptID'].upper()
  if re.fullmatch(r'C\d{7}|CN\d{6}',cid):ids.add('MedGen:'+cid)
  lm[cid].update(ids);ln[norm(a['DiseaseName'])].update(ids);gd.update(ids)
 lookup[g['gene']]=lm;names[g['gene']]=ln;gene_ids[g['gene']]=gd
 old_panel[g['gene']]={norm(x):'OMIM:'+g['panel_explicit_OMIM_disease_ids'][0] for x in g['source_panel_phenotypes']} if len(g['panel_explicit_OMIM_disease_ids'])==1 else {}
ann=collections.defaultdict(list)
for a in ext['preserved_raw_target_phenotype_annotation_rows']:ann[a['database_id']].append(a)
def disease_tag(ids):
 rs=[r for d in ids for r in ann.get(d,[])];pos=neg=False
 for r in rs:
  negative=r['qualifier'].strip().upper()=='NOT' or r['frequency'].strip() in {'HP:0040285','0','0%'} or bool(re.fullmatch('0/\d+',r['frequency'].strip()))
  pos|=not negative;neg|=negative
 return 'mixed_positive_negative' if pos and neg else 'gene_disease_positive_stone_NC_annotation' if pos else 'mapped_named_disease_without_positive_annotation'
hpo_primary={x['id'] for x in ext['all_is_a_descendant_terms']};hpo_alt={a:x['id'] for x in ext['all_is_a_descendant_terms'] for a in x['alt_ids']};stone_names={norm(x['name']) for x in ext['all_is_a_descendant_terms']};root_med={'C0392525','C0027709'}
with Path(ext['sources']['NCBI_ConceptID_history.tsv']['path']).open(encoding='utf-8',newline='') as f:hist={r['Previous ConceptID']:r for r in csv.DictReader(f,delimiter='\t')}
def scope(r,g,mode,allow_names):
 rx=r'\b(?:C\d{7}|CN\d{6})\b' if mode=='B' else r'\bC\d{7}\b';parts=[]
 cs=[c.strip() for c in re.split('[|;]',r['ReportedPhenotypeInfo']) if c.strip()] or ['']
 for c in cs:
  mm={x.upper() for x in re.findall(rx,c,re.I)};hp={x.upper() for x in re.findall('HP:\d{7}',c,re.I)};om={'OMIM:'+x for x in re.findall('OMIM[#: ]*(\d{6})',c,re.I)};ids=om&gene_ids[g];paths=[]
  if ids:paths.append('explicit_gene_associated_OMIM')
  for m in mm:
   found=lookup[g].get(m,set())
   if mode=='A':found={x for x in found if not x.startswith('MedGen:CN')}
   if found:ids|=found;paths.append('explicit_official_gene_associated_MedGen')
  clean=re.sub(r'^\s*(?:(?:C\d{7}|CN\d{6})|na)\s*:\s*','',c,flags=re.I) if mode=='B' else re.sub(r'^\s*(?:C\d{7}|na)\s*:\s*','',c,flags=re.I)
  hn={hpo_alt[x] for x in hp if x in hpo_alt} if mode=='B' else set()
  direct=bool(mm&root_med or hp&hpo_primary or hn or allow_names and norm(clean) in stone_names)
  if direct:tag='direct_stone_NC_condition';paths.append('direct_target_phenotype')
  else:
   nn=names[g].get(norm(clean),set()) if allow_names else set()
   if mode=='A':nn={x for x in nn if not x.startswith('MedGen:CN')}
   if nn:ids|=nn;paths.append('exact_official_DiseaseName')
   tag=disease_tag(ids) if ids else 'unknown_unresolved_condition'
  parts.append({'raw_component':c,'tag':tag,'gene_associated_disease_ids':sorted(ids),'MedGen_ids':sorted(mm),'HPO_ids':sorted(hp),'OMIM_ids':sorted(om),'HPO_alternative_to_primary':{x:hpo_alt[x] for x in hp if mode=='B' and x in hpo_alt},'identity_paths':sorted(set(paths))})
 sm={x.upper() for x in re.findall(rx,r['SubmittedPhenotypeInfo'],re.I)};sh={x.upper() for x in re.findall('HP:\d{7}',r['SubmittedPhenotypeInfo'],re.I)};so={'OMIM:'+x for x in re.findall('OMIM[#: ]*(\d{6})',r['SubmittedPhenotypeInfo'],re.I)};supp=[];si=so&gene_ids[g]
 if sm&root_med or sh&hpo_primary or mode=='B' and any(x in hpo_alt for x in sh):supp.append('direct_stone_NC_condition')
 for m in sm:
  found=lookup[g].get(m,set());si|=({x for x in found if not x.startswith('MedGen:CN')} if mode=='A' else found)
 if si:supp.append(disease_tag(si))
 tags={x['tag'] for x in parts}|set(supp)
 return {'scope':next(iter(tags)) if len(tags)==1 else 'mixed_scope_or_unresolved_components','component_evidence':parts,'supplementary_explicit_ID_scope':supp,'supplementary_gene_associated_disease_ids':sorted(si)}
def variant_category(scvs,key):
 tags={s[key]['scope'] for s in scvs if s['tier']=='PLP'}
 return next(iter(tags)) if len(tags)==1 else 'mixed_scope_across_PLP_submissions' if tags else 'no_resolved_PLP_contributing_SCV'
aud=[];record_paths=collections.Counter();transition=collections.Counter();strict_transition=collections.Counter();hist_hits=[];examples={};gene_counts=collections.defaultdict(collections.Counter);all_sources=0
with gzip.open(source,'rt',encoding='utf-8') as f:
 for line in f:
  v=json.loads(line);g=v['gene'];ss=[]
  for s in v['all_source_contributing_SCV_evidence']:
   all_sources+=1;raw=s['ReportedPhenotypeInfo']+' | '+s['SubmittedPhenotypeInfo'];b=dict(s)
   for mode in ['A','B']:
    b[mode+'_with_official_names']=scope(s,g,mode,True);b[mode+'_identifiers_only']=scope(s,g,mode,False)
   pcs=[re.sub(r'^\s*(?:C\d{7}|na)\s*:\s*','',x,flags=re.I) for x in re.split('[|;]',s['ReportedPhenotypeInfo']) if x.strip()]
   shortcut=[{'raw_component':c,'assumed_disease':old_panel[g][norm(c)],'also_official_name':norm(c) in names[g]} for c in pcs if norm(c) in old_panel[g]]
   cn=sorted(set(re.findall(r'\bCN\d{6}\b',raw,re.I)));hp_alt=sorted(set(re.findall('HP:\d{7}',raw,re.I))&set(hpo_alt));history_ids=sorted({x.upper() for x in re.findall(r'\bC(?:N)?\d+\b',raw,re.I)}&set(hist))
   b['postreview_identity_diagnostics']={'original_sole_panel_OMIM_name_matches':shortcut,'CN_identifiers':cn,'target_HPO_alternative_ids':hp_alt,'legacy_history_identifiers_not_auto_merged':history_ids}
   for label,ok in [('sole_panel_name_path',bool(shortcut)),('sole_panel_name_exclusive_of_official_name',any(not x['also_official_name'] for x in shortcut)),('CN_format',bool(cn)),('HPO_alt',bool(hp_alt)),('legacy_history',bool(history_ids))]:
    if ok:record_paths[label]+=1;examples.setdefault(label,{'VariationID':v['VariationID'],'gene':g,'SCV':s['SCV'],'ReportedPhenotypeInfo':s['ReportedPhenotypeInfo'],'SubmittedPhenotypeInfo':s['SubmittedPhenotypeInfo'],'source_line':s['source_line'],'source_line_sha256':s['source_line_sha256']})
   if history_ids:hist_hits.append({'VariationID':v['VariationID'],'SCV':s['SCV'],'gene':g,'source_line':s['source_line'],'histories':[hist[x] for x in history_ids]})
   ss.append(b)
  a=variant_category(ss,'A_with_official_names');b=variant_category(ss,'B_with_official_names');ast=variant_category(ss,'A_identifiers_only');bst=variant_category(ss,'B_identifiers_only')
  q={**v,'all_source_contributing_SCV_evidence':ss,'original_scope':v['PLP_evidence_scope'],'A_corrected_official_names_scope':a,'B_corrected_source_supported_scope':b,'A_identifiers_only_scope':ast,'B_identifiers_only_scope':bst};aud.append(q)
  transition[(v['PLP_evidence_scope'],a,b)]+=1;strict_transition[(v['identifiers_only_PLP_evidence_scope'],ast,bst)]+=1;gene_counts[g][b]+=1
assert len(aud)==4078
dest=D/'all4078_original_and_corrected_condition_scope_with_sources.jsonl.gz'
with gzip.open(dest,'wt',encoding='utf-8',compresslevel=3,newline='\n') as f:
 for v in aud:f.write(json.dumps(v,ensure_ascii=False,separators=(',',':'))+'\n')
counts=lambda key:dict(collections.Counter(v[key] for v in aud))
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'ALL_POSTREVIEW_SCOPE_CORRECTIONS_COMPUTED_PENDING_INDEPENDENT_RAW_DICTIONARY_QA','protocol_sha256':sha(p),'source_sha256':sha(source),'all4078_retained':True,'all_source_contributing_records':all_sources,'all_PL P_source_records'.replace(' ',''):sum(s['tier']=='PLP' for v in aud for s in v['all_source_contributing_SCV_evidence']),'original_counts':counts('original_scope'),'A_remove_sole_panel_shortcut_counts':counts('A_corrected_official_names_scope'),'B_corrected_source_supported_counts':counts('B_corrected_source_supported_scope'),'B_identifiers_only_counts':counts('B_identifiers_only_scope'),'changed_original_to_A_variants':sum(v['original_scope']!=v['A_corrected_official_names_scope'] for v in aud),'changed_original_to_B_variants':sum(v['original_scope']!=v['B_corrected_source_supported_scope'] for v in aud),'contributor_record_path_occurrences':dict(record_paths),'new_inferential_tests':0,'primary_and_current_dated_results_changed':False,'original_results_overwritten':False,'satisfactory_manuscript_mainline':False,'output':{'path':str(dest),'bytes':dest.stat().st_size,'sha256':sha(dest)},'limits':'Post-result implementation repair and explicit expanded ID coverage. No clinical-error/variant-function/patient-risk conclusion. NCBI history not automatically merged; unresolved remains unknown.'}
dump(O/'logs/repair_results_pending_independent_QA.json',out);dump(O/'metadata/all_observed_legacy_ID_histories.json',hist_hits);dump(O/'metadata/deterministic_first_in_complete_sorted_cohort_mapping_path_examples.json',examples)
for path,ct,fields in [(O/'tables/all4078_original_A_B_scope_transitions.tsv',transition,['original_scope','A_scope','B_scope','variants']),(O/'tables/all4078_direct_identifiers_scope_transitions.tsv',strict_transition,['original_scope','A_scope','B_scope','variants'])]:
 with path.open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(fields)
  for k,n in sorted(ct.items()):w.writerow([*k,n])
with (O/'tables/all35_corrected_scope_counts.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['gene','official_PLP_variants','corrected_scope_counts'])
 for g in sorted(lookup):w.writerow([g,sum(gene_counts[g].values()),json.dumps(dict(gene_counts[g]),ensure_ascii=False)])
print('STAGE42_COMPLETE_ALL4078_POSTREVIEW_REPAIR_PENDING_QA',json.dumps(out,ensure_ascii=False),flush=True)
