"""Fixed descriptive census; no tests, fitted classifier or selected gene promotion."""
from pathlib import Path
import json, gzip, hashlib, datetime, re, unicodedata, collections, csv
ROOT=Path(__file__).resolve().parents[1]; O=ROOT/'stage39'; P=O/'logs/protocol.json'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def norm(s):return ' '.join(unicodedata.normalize('NFKC',s or '').split()).casefold()
def sym(s):return set(filter(None,re.split(r'[;,|\s]+',s.strip())))
def hgnc(s):return set(re.findall(r'HGNC:\d+',s))
def readjgz(p):
 with gzip.open(p,'rt',encoding='utf-8') as f:
  for line in f:yield json.loads(line)
def writejgz(p,rows):
 with gzip.open(p,'wt',encoding='utf-8',compresslevel=3,newline='\n') as f:
  for r in rows:f.write(json.dumps(r,ensure_ascii=False,separators=(',',':'))+'\n')
def tsv(p,rows,fields):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
assert sha(P)==(O/'logs/protocol_sha256.txt').read_text().strip()
protocol=json.loads(P.read_text(encoding='utf-8'));D=Path(protocol['storage'])
complete=json.loads((O/'logs/full_source_extraction_complete.json').read_text(encoding='utf-8'));assert complete['status']=='COMPLETE_ALL_ENDPOINT_AND13MONTH_FULL_SOURCE_EXTRACTION' and complete['protocol_sha256']==sha(P)
for r in complete['endpoint_receipts']+complete['monthly_receipts']:assert Path(r['file']).stat().st_size==r['bytes'] and sha(r['file'])==r['sha256']
ext=json.loads(Path(protocol['external_dictionary']['path']).read_text(encoding='utf-8'));assert sha(protocol['external_dictionary']['path'])==protocol['external_dictionary']['sha256']
external={g['gene']:g for g in ext['all35_external_gene_records']}; h2g={g['HGNC']:g['gene'] for g in external.values()}
major={norm(v):k for k,vals in protocol['major_class_mapping'].items() for v in vals}
five={norm(x):x for x in ['Pathogenic','Likely pathogenic','Uncertain significance','Likely benign','Benign']}
def tier(r):return major.get(norm(r['ClinicalSignificance']))
def condition(r):return (norm(r['SubmittedPhenotypeInfo']),norm(r['ReportedPhenotypeInfo']))
def medgen_ids(r):return frozenset(re.findall(r'\bC\d{7}\b',r['ReportedPhenotypeInfo'],flags=re.I))
def complete_reported_ids(r):
 components=[s.strip() for s in re.split(r'[|;]',r['ReportedPhenotypeInfo']) if s.strip()]
 if not components or any(not re.search(r'\bC\d{7}\b',s,re.I) for s in components):return None
 return frozenset(x.upper() for x in medgen_ids(r))
def version(r):return r['_SCV_version']
def unspecified(r):
 parts=[re.sub(r'^\s*(?:C\d{7}|na)\s*:\s*','',p,flags=re.I) for p in re.split(r'[|;]',r['ReportedPhenotypeInfo'])]
 return bool(parts) and all(norm(p) in ['','not provided','not specified','not provided or not specified','na','-'] for p in parts)
def five_sequence_eligible(rows):
 return all(norm(r['ClinicalSignificance']) in five for r in rows) and not any(version(a)==version(b) and norm(a['ClinicalSignificance'])!=norm(b['ClinicalSignificance']) for a,b in zip(rows,rows[1:]))

# Complete endpoint gene/allele identity audit.
endpoints=[]
for receipt in complete['endpoint_receipts']:
 grouped=collections.defaultdict(list)
 for r in readjgz(receipt['file']):grouped[r['VariationID']].append(r)
 endpoints.append(grouped)
all_vids=set(endpoints[0])|set(endpoints[1]);variant_identities={};stable={};current_aggregate={}
for vid in sorted(all_vids,key=int):
 summaries=[]
 for data in endpoints:
  rows=data.get(vid,[]);hg=set().union(*(hgnc(r['HGNC_ID']) for r in rows));symbols=set().union(*(sym(r['GeneSymbol']) for r in rows));alleles={r['AlleleID'] for r in rows}
  summaries.append({'rows':len(rows),'HGNC_IDs':sorted(hg),'gene_symbols':sorted(symbols),'AlleleIDs':sorted(alleles),'missing_HGNC_rows':sum(not hgnc(r['HGNC_ID']) for r in rows),'source_lines':[r['_source_line'] for r in rows],'assemblies':sorted({r['Assembly'] for r in rows})})
 flags=[]
 if not all(s['rows'] for s in summaries):flags.append('endpoint_absent')
 if any(s['missing_HGNC_rows'] for s in summaries):flags.append('missing_HGNC_association')
 if any(len(s['HGNC_IDs'])!=1 or len(s['gene_symbols'])>1 for s in summaries):flags.append('multigene_or_unresolved_gene')
 if any(len(s['AlleleIDs'])!=1 for s in summaries):flags.append('complex_or_unresolved_allele')
 if summaries[0]['HGNC_IDs']!=summaries[1]['HGNC_IDs'] or summaries[0]['AlleleIDs']!=summaries[1]['AlleleIDs']:flags.append('endpoint_identity_changed')
 gs=set().union(*(set(s['HGNC_IDs']) for s in summaries)); gene=h2g.get(next(iter(gs))) if len(gs)==1 else None
 if gene is None:flags.append('not_unique_external_HGNC')
 rec={'VariationID':vid,'gene':gene,'endpoint_identity':summaries,'flags':flags};variant_identities[vid]=rec
 if not flags:
  stable[vid]=gene
  vals={r['ClinicalSignificance'] for r in endpoints[1][vid]};scvs={b for r in endpoints[1][vid] for b in re.findall(r'SCV\d+',r['SCVsForAggregateGermlineClassification'])}
  current_aggregate[vid]={'official_germline_values':sorted(vals),'major_tiers':sorted({major.get(norm(v),'OUTSIDE_MAJOR') for v in vals}),'source_declared_germline_contributor_bases':sorted(scvs),'rows':[{'line':r['_source_line'],'sha256':r['_source_line_sha256'],'SCVsForAggregateGermlineClassification':r['SCVsForAggregateGermlineClassification']} for r in endpoints[1][vid]]}
writejgz(D/'all_endpoint_variant_identity_and_exclusions.jsonl.gz',variant_identities.values())
del endpoints
print('ENDPOINT_IDENTITY_AUDIT_COMPLETE',len(all_vids),'STABLE_PRIMARY_VARIANTS',len(stable),flush=True)

# Load compact records; all original descriptions/source lines remain in the frozen extraction.
fields=['VariationID','ClinicalSignificance','DateLastEvaluated','SubmittedPhenotypeInfo','ReportedPhenotypeInfo','ReviewStatus','CollectionMethod','Submitter','SCV','SubmittedGeneSymbol','ContributesToAggregateClassification','SomaticClinicalImpact','Oncogenicity','_SCV_base','_SCV_version','_source_line','_source_line_sha256']
month_names=[];months=[];malformed=[];exact_duplicates=0
for receipt in complete['monthly_receipts']:
 name=Path(receipt['file']).name.split('.txt.gz')[0];month_names.append(name);groups=collections.defaultdict(dict)
 for r in readjgz(receipt['file']):
  if not r['_SCV_base']:malformed.append({'month':name,'source_line':r['_source_line'],'sha256':r['_source_line_sha256'],'SCV':r['SCV']});continue
  k=r['_source_line_sha256']
  if k in groups[r['_SCV_base']]:exact_duplicates+=1
  groups[r['_SCV_base']][k]={f:r[f] for f in fields}
 months.append({k:list(v.values()) for k,v in groups.items()})
 assert sum(len(v) for v in months[-1].values())+exact_duplicates>=0
 print('MONTH_COMPACT_RECORDS_LOADED',name,len(groups),flush=True)
assert len(months)==13
universe=set().union(*(set(m) for m in months));audit=[];events=[];pair_events=[];primary=[];no_name=[];id_condition=[];exclusion_counts=collections.Counter();nonmajor_counts=collections.Counter();five_events=[];endpoint_events=[];pairs_eligible=0
def sequence_flags(seq, require_name=True, id_conditions=False):
 flags=[]
 if any(not rows for rows in seq):flags.append('missing_month')
 if any(len(rows)>1 for rows in seq):flags.append('ambiguous_multiple_rows_same_SCV_month')
 # Other identity/condition flags use every available row; ambiguous or absent
 # records never become eligible. Version sequence checks require unique rows.
 rows=[r for rs in seq for r in rs];vids={r['VariationID'] for r in rows}
 if not rows:return flags
 if len(vids)!=1:flags.append('VariationID_changed')
 if len(vids)!=1 or next(iter(vids)) not in stable:flags.append('not_stable_primary_variant')
 if id_conditions:
  ids=[complete_reported_ids(r) for r in rows]
  if any(s is None for s in ids):flags.append('reported_condition_has_unresolved_component')
  elif len(set(ids))!=1:flags.append('reported_MedGen_ID_set_changed')
 elif len({condition(r) for r in rows})!=1:flags.append('raw_condition_pair_changed')
 if require_name and len({norm(r['Submitter']) for r in rows})!=1:flags.append('submitter_name_changed')
 if any(tier(r) is None for r in rows):flags.append('outside_major_tier')
 if any(version(r) is None for r in rows):flags.append('missing_numeric_version')
 elif not any(len(rs)>1 for rs in seq):
  if any(version(b)<version(a) for a,b in zip(rows,rows[1:])):flags.append('version_regression')
  if any(version(b)==version(a) and tier(b)!=tier(a) for a,b in zip(rows,rows[1:])):flags.append('same_version_major_tier_inconsistency')
 return flags
def get_events(base,rows):
 return [{'SCV_base':base,'VariationID':a['VariationID'],'gene':stable[a['VariationID']],'left_month':month_names[i],'right_month':month_names[i+1],'left_tier':tier(a),'right_tier':tier(b),'left_version':version(a),'right_version':version(b),'submitter':a['Submitter'],'left_line':a['_source_line'],'right_line':b['_source_line'],'left_line_sha256':a['_source_line_sha256'],'right_line_sha256':b['_source_line_sha256'],'SubmittedPhenotypeInfo':a['SubmittedPhenotypeInfo'],'ReportedPhenotypeInfo':a['ReportedPhenotypeInfo']} for i,(a,b) in enumerate(zip(rows,rows[1:])) if tier(a)!=tier(b) and version(b)>version(a)]
for base in sorted(universe):
 seq=[m.get(base,[]) for m in months];flags=sequence_flags(seq);exclusion_counts.update(flags)
 for rs in seq:
  for r in rs:
   if tier(r) is None:nonmajor_counts[r['ClinicalSignificance']]+=1
 record={'SCV_base':base,'present_months':sum(bool(x) for x in seq),'flags':flags,'primary_eligible':not flags,'source_refs':[{'month':month_names[i],'rows':[{'VariationID':r['VariationID'],'SCV':r['SCV'],'source_line':r['_source_line'],'sha256':r['_source_line_sha256']} for r in rs]} for i,rs in enumerate(seq) if rs]}
 if not flags:
  rows=[x[0] for x in seq];primary.append(base);ev=get_events(base,rows);events.extend(ev);record.update(VariationID=rows[0]['VariationID'],gene=stable[rows[0]['VariationID']],n_major_updates=len(ev),submitter=rows[0]['Submitter'],baseline_tier=tier(rows[0]),current_tier=tier(rows[-1]),unspecified_condition=unspecified(rows[0]),version_increase_without_major_change_intervals=sum(version(b)>version(a) and tier(a)==tier(b) for a,b in zip(rows,rows[1:])),evaluation_date_changed_without_major_change_intervals=sum(b['DateLastEvaluated']!=a['DateLastEvaluated'] and tier(a)==tier(b) for a,b in zip(rows,rows[1:])))
  if tier(rows[0])!=tier(rows[-1]) and version(rows[-1])>version(rows[0]):endpoint_events.append(base)
  if five_sequence_eligible(rows):
   for i,(a,b) in enumerate(zip(rows,rows[1:])):
    if norm(a['ClinicalSignificance'])!=norm(b['ClinicalSignificance']) and version(b)>version(a):five_events.append({'SCV_base':base,'left_month':month_names[i],'left_class':a['ClinicalSignificance'],'right_class':b['ClinicalSignificance']})
 if not sequence_flags(seq,require_name=False):no_name.append({'SCV_base':base,'events':len(get_events(base,[x[0] for x in seq]))})
 if not sequence_flags(seq,id_conditions=True):id_condition.append({'SCV_base':base,'events':len(get_events(base,[x[0] for x in seq]))})
 for i in range(12):
  pair=seq[i:i+2]
  if not sequence_flags(pair):
   pairs_eligible+=1;a,b=pair[0][0],pair[1][0]
   if tier(a)!=tier(b) and version(b)>version(a):pair_events.append({'SCV_base':base,'VariationID':a['VariationID'],'gene':stable[a['VariationID']],'left_month':month_names[i],'right_month':month_names[i+1],'left_tier':tier(a),'right_tier':tier(b),'left_version':version(a),'right_version':version(b),'left_line':a['_source_line'],'right_line':b['_source_line'],'left_line_sha256':a['_source_line_sha256'],'right_line_sha256':b['_source_line_sha256']})
 audit.append(record)
writejgz(D/'all_selected_SCV_eligibility_and_source_refs.jsonl.gz',audit)
writejgz(D/'all_primary_major_events.jsonl.gz',events)
writejgz(D/'all_five_level_classification_events.jsonl.gz',five_events)
writejgz(D/'all_adjacent_pair_major_events.jsonl.gz',pair_events)
dump(O/'metadata/all_source_nonmajor_germline_terms.json',{'terms_by_selected_record_month_count':dict(nonmajor_counts),'malformed_SCV_rows':malformed,'exact_identical_source_duplicates_collapsed':exact_duplicates})

# Fixed current dated contrast; no arbitrary selection among mixed tied dates.
current=collections.defaultdict(list)
for base,rs in months[-1].items():
 if len(rs)==1 and rs[0]['VariationID'] in stable:current[rs[0]['VariationID']].append(rs[0])
current_release_date=max(i['release_date'][:10] for i in protocol['inputs'] if 'submission_summary' in i['file'])
def date_value(raw):
 for fmt in ['%b %d, %Y','%Y-%m-%d','%m/%d/%Y']:
  try:
   d=datetime.datetime.strptime(raw.strip(),fmt).date().isoformat()
   return d if d<=current_release_date else None
  except ValueError:pass
 return None
trace_by_variant=collections.defaultdict(collections.Counter)
for a in audit:
 vids={r['VariationID'] for ref in a['source_refs'] for r in ref['rows']}
 for vid in vids:
  trace_by_variant[vid]['selected_SCV_bases']+=1
  trace_by_variant[vid]['continuous_primary_eligible_SCV_bases']+=a['primary_eligible']
  for flag in a['flags']:trace_by_variant[vid][flag]+=1
def pair_category(a,b,predicate):
 if not a or not b:return 'NOT_AVAILABLE'
 vals=[predicate(x,y) for x in a for y in b]
 return 'ALL_SAME' if all(vals) else 'ALL_DIFFERENT' if not any(vals) else 'MIXED'
contrasts=[]
for vid in sorted(stable,key=int):
 rows=current.get(vid,[]);valid=[(date_value(r['DateLastEvaluated']),r) for r in rows if tier(r) and date_value(r['DateLastEvaluated'])];status='insufficient_dated_major_assertions';a=b=[]
 if len(valid)>=2:
  lo=min(x[0] for x in valid);hi=max(x[0] for x in valid);a=[r for d,r in valid if d==lo];b=[r for d,r in valid if d==hi]
  if lo==hi:status='all_dated_assertions_same_date'
  elif len({tier(r) for r in a})>1 or len({tier(r) for r in b})>1:status='mixed_major_tiers_at_tied_endpoint_date'
  elif tier(a[0])==tier(b[0]):status='determinate_same_tier'
  else:status='determinate_different_tier'
 record={'VariationID':vid,'gene':stable[vid],'status':status,'n_current_unique_SCVs':len(rows),'n_valid_dated_major_SCVs':len(valid),'earliest_date':min([x[0] for x in valid],default=None),'latest_date':max([x[0] for x in valid],default=None),'earliest_SCVs':[r['SCV'] for r in a],'latest_SCVs':[r['SCV'] for r in b],'earliest_tiers':sorted({tier(r) for r in a}),'latest_tiers':sorted({tier(r) for r in b}), 'all_endpoint_pairs_different_SCV':bool(a and b) and all(x['_SCV_base']!=y['_SCV_base'] for x in a for y in b),'all_endpoint_pairs_same_submitter':bool(a and b) and all(norm(x['Submitter'])==norm(y['Submitter']) for x in a for y in b),'all_endpoint_pairs_same_raw_condition':bool(a and b) and all(condition(x)==condition(y) for x in a for y in b),'SCV_pair_identity':pair_category(a,b,lambda x,y:x['_SCV_base']==y['_SCV_base']),'submitter_pair_identity':pair_category(a,b,lambda x,y:norm(x['Submitter'])==norm(y['Submitter'])),'raw_condition_pair_identity':pair_category(a,b,lambda x,y:condition(x)==condition(y)),'direct_trace_availability_and_exclusions':dict(trace_by_variant[vid]),'primary_direct_update_events':[e for e in events if e['VariationID']==vid]}
 contrasts.append(record)
writejgz(D/'all_current_dated_assertion_contrasts.jsonl.gz',contrasts)

# Condition-level evidence scope. Raw disease associations/qualifiers remain in external source.
annotation=collections.defaultdict(list)
for a in ext['preserved_raw_target_phenotype_annotation_rows']:annotation[a['database_id']].append(a)
def excluded_annotation(a):
 f=a['frequency'].strip();return a['qualifier'].strip().upper()=='NOT' or f in ['HP:0040285','0','0%'] or bool(re.fullmatch(r'0/\d+',f))
def disease_scope(ids):
 rows=[a for d in ids for a in annotation.get(d,[])]
 pos=any(not excluded_annotation(a) for a in rows);neg=any(excluded_annotation(a) for a in rows)
 return 'mixed_positive_negative' if pos and neg else 'gene_disease_positive_stone_NC_annotation' if pos else 'mapped_named_disease_without_positive_annotation'
lookup={};names={};gene_diseases={}
for gene,g in external.items():
 lm=collections.defaultdict(set);ln=collections.defaultdict(set);gd=set(g['HPO_Mendelian_gene_disease_ids'])|{'OMIM:'+x for x in g['panel_explicit_OMIM_disease_ids']}
 for r in g['official_associated_gene_condition_rows']:
  ids={'OMIM:'+x for x in re.findall(r'\b[1-6]\d{5}\b',r['DiseaseMIM'])}
  if re.fullmatch(r'C\d{7}',r['ConceptID'],re.I):ids.add('MedGen:'+r['ConceptID'].upper())
  lm[r['ConceptID'].upper()].update(ids);ln[norm(r['DiseaseName'])].update(ids);gd.update(ids)
 for n in g['source_panel_phenotypes']:
  if len(g['panel_explicit_OMIM_disease_ids'])==1:ln[norm(n)].add('OMIM:'+g['panel_explicit_OMIM_disease_ids'][0])
 lookup[gene]=lm;names[gene]=ln;gene_diseases[gene]=gd
root_medgen={'C0392525','C0027709'};hpo_ids={t['id'] for t in ext['all_is_a_descendant_terms']};stone_names={norm(t['name']) for t in ext['all_is_a_descendant_terms']}
def scope_record(r,gene,allow_names):
 text=r['ReportedPhenotypeInfo'];components=[s.strip() for s in re.split(r'[|;]',text) if s.strip()];parts=[]
 if not components:components=['']
 for comp in components:
  med={x.upper() for x in re.findall(r'\bC\d{7}\b',comp,re.I)};hp={x.upper() for x in re.findall(r'HP:\d{7}',comp,re.I)};om={'OMIM:'+x for x in re.findall(r'OMIM[#: ]*(\d{6})',comp,re.I)}
  cleaned=re.sub(r'^\s*(?:C\d{7}|na)\s*:\s*','',comp,flags=re.I)
  ids=set(om & gene_diseases[gene])
  for c in med:ids.update(lookup[gene].get(c,set()))
  if med&root_medgen or hp&hpo_ids or (allow_names and norm(cleaned) in stone_names):tag='direct_stone_NC_condition'
  else:
   if allow_names:ids.update(names[gene].get(norm(cleaned),set()))
   tag=disease_scope(ids) if ids else 'unknown_unresolved_condition'
  parts.append({'raw_component':comp,'tag':tag,'gene_associated_disease_ids':sorted(ids),'MedGen_ids':sorted(med),'HPO_ids':sorted(hp),'OMIM_ids':sorted(om)})
 # Explicit submitted IDs may supplement identity; they do not replace unresolved reported conditions.
 submitted=r['SubmittedPhenotypeInfo'];submed={x.upper() for x in re.findall(r'\bC\d{7}\b',submitted,re.I)};subhp={x.upper() for x in re.findall(r'HP:\d{7}',submitted,re.I)};subom={'OMIM:'+x for x in re.findall(r'OMIM[#: ]*(\d{6})',submitted,re.I)}
 supplementary=[]
 if submed&root_medgen or subhp&hpo_ids:supplementary.append('direct_stone_NC_condition')
 ids=set(subom&gene_diseases[gene])
 for c in submed:ids.update(lookup[gene].get(c,set()))
 if ids:supplementary.append(disease_scope(ids))
 tags={x['tag'] for x in parts}|set(supplementary)
 tag=next(iter(tags)) if len(tags)==1 else 'mixed_scope_or_unresolved_components'
 return {'scope':tag,'component_evidence':parts,'supplementary_explicit_ID_scope':supplementary}
portability=[];missing_contributors=[];contributing_disagreements=[]
for vid,agg in current_aggregate.items():
 if agg['major_tiers']!=['PLP']:continue
 linked=[]
 for base in agg['source_declared_germline_contributor_bases']:
  rs=months[-1].get(base,[])
  if len(rs)!=1 or rs[0]['VariationID']!=vid:
   missing_contributors.append({'VariationID':vid,'SCV_base':base,'matched_rows':len(rs)});continue
  r=rs[0]; sc=scope_record(r,stable[vid],True);strict=scope_record(r,stable[vid],False)
  linked.append({'SCV_base':base,'SCV':r['SCV'],'tier':tier(r),'raw_class':r['ClinicalSignificance'],'submitter':r['Submitter'],'ReviewStatus':r['ReviewStatus'],'ContributesToAggregateClassification':r['ContributesToAggregateClassification'],'ReportedPhenotypeInfo':r['ReportedPhenotypeInfo'],'SubmittedPhenotypeInfo':r['SubmittedPhenotypeInfo'],'scope_with_external_exact_names':sc,'scope_direct_identifiers_only':strict,'source_line':r['_source_line'],'source_line_sha256':r['_source_line_sha256']})
  if norm(r['ContributesToAggregateClassification']) not in ['yes','true','1','germline']:contributing_disagreements.append({'VariationID':vid,'SCV_base':base,'source_value':r['ContributesToAggregateClassification']})
 plp=[s for s in linked if s['tier']=='PLP'];tags={s['scope_with_external_exact_names']['scope'] for s in plp};strict_tags={s['scope_direct_identifiers_only']['scope'] for s in plp}
 portability.append({'VariationID':vid,'gene':stable[vid],'official_aggregate':agg,'n_linked_source_contributors':len(linked),'n_linked_PLP_SCVs':len(plp),'PLP_evidence_scope':next(iter(tags)) if len(tags)==1 else 'mixed_scope_across_PLP_submissions' if tags else 'no_resolved_PLP_contributing_SCV','identifiers_only_PLP_evidence_scope':next(iter(strict_tags)) if len(strict_tags)==1 else 'mixed_scope_across_PLP_submissions' if strict_tags else 'no_resolved_PLP_contributing_SCV','all_source_contributing_SCV_evidence':linked})
writejgz(D/'all_latest_PLP_aggregate_condition_scope_and_sources.jsonl.gz',portability)
dump(O/'logs/contributor_source_linkage_QA.json',{'missing_or_ambiguous':missing_contributors,'source_contribution_value_disagreements_for_review':contributing_disagreements,'no_own_aggregate_policy_inferred':True})

update_bases={e['SCV_base'] for e in events};update_variants={e['VariationID'] for e in events};transition=collections.Counter((e['left_tier'],e['right_tier']) for e in events)
per_gene=[]
for gene in sorted(external):
 eligible=[a for a in audit if a.get('primary_eligible') and a.get('gene')==gene];ge=[e for e in events if e['gene']==gene];ps=[v for v in portability if v['gene']==gene]
 per_gene.append({'gene':gene,'stable_endpoint_variants':sum(g==gene for g in stable.values()),'continuous_primary_SCVs':len(eligible),'SCVs_with_major_update':len({e['SCV_base'] for e in ge}),'major_update_events':len(ge),'latest_official_PLP_variants':len(ps),'condition_scope_counts':dict(collections.Counter(v['PLP_evidence_scope'] for v in ps))})
tsv(O/'tables/all35_gene_descriptive_census.tsv',[{**g,'condition_scope_counts':json.dumps(g['condition_scope_counts'],ensure_ascii=False)} for g in per_gene],list(per_gene[0]))
submitters=collections.Counter(a['submitter'] for a in audit if a.get('primary_eligible'));es=collections.Counter(e['submitter'] for e in events)
tsv(O/'tables/all_primary_submitter_concentration.tsv',[{'submitter':s,'eligible_SCVs':n,'update_events':es[s],'updated_distinct_SCVs':len({e['SCV_base'] for e in events if e['submitter']==s})} for s,n in submitters.most_common()],['submitter','eligible_SCVs','update_events','updated_distinct_SCVs'])
five_eligible=sum(five_sequence_eligible([m[b][0] for m in months]) for b in primary)
criteria_eligible=[b for b in primary if all(norm(m[b][0]['ReviewStatus']).startswith('criteria provided') or norm(m[b][0]['ReviewStatus']) in ['reviewed by expert panel','practice guideline'] for m in months)]
contribution_values=collections.Counter(r['ContributesToAggregateClassification'] for m in months for rs in m.values() for r in rs)
contribution_eligible=[b for b in primary if all(norm(m[b][0]['ContributesToAggregateClassification']) in ['yes','true','1','germline'] for m in months)]
dump(O/'metadata/all_source_contribution_flags.json',{'source_value_counts_selected_record_months':dict(contribution_values),'recognized_positive_values':['yes','true','1','germline'],'raw_values_preserved':True})
result={'utc':now(),'status':'FROZEN_DESCRIPTIVE_ENDPOINTS_COMPUTED_PENDING_INDEPENDENT_SOURCE_QA','protocol_sha256':sha(P),'selected_endpoint_VariationIDs':len(all_vids),'stable_primary_single_gene_variants':len(stable),'selected_SCV_bases':len(universe),
 'primary':{'eligible_continuous_SCVs':len(primary),'SCVs_with_any_major_version_increase_update':len(update_bases),'fraction':len(update_bases)/len(primary) if primary else None,'distinct_variants_updated':len(update_variants),'major_update_events':len(events),'all6_directed_transitions':{a+'->'+b:transition[(a,b)] for a in ['PLP','VUS','BLB'] for b in ['PLP','VUS','BLB'] if a!=b},'endpoint_only_updated_SCVs':len(endpoint_events),'nonexclusive_exclusion_flags':dict(exclusion_counts),'unspecified_condition_eligible_SCVs':sum(a.get('primary_eligible') and a.get('unspecified_condition') for a in audit),'version_increase_without_major_change_intervals':sum(a.get('version_increase_without_major_change_intervals',0) for a in audit),'evaluation_date_changed_without_major_change_intervals':sum(a.get('evaluation_date_changed_without_major_change_intervals',0) for a in audit)},
 'fixed_sensitivities':{'five_standard_levels_eligible':five_eligible,'five_level_update_events':len(five_events),'without_submitter_constancy_eligible':len(no_name),'without_submitter_constancy_updated_SCVs':sum(x['events']>0 for x in no_name),'complete_reported_MedGen_ID_set_eligible':len(id_condition),'complete_reported_ID_set_updated_SCVs':sum(x['events']>0 for x in id_condition),'adjacent_eligible_SCV_intervals':pairs_eligible,'adjacent_update_events':len(pair_events),'criteria_provided_continuous_SCVs':len(criteria_eligible),'criteria_provided_updated_SCVs':sum(b in update_bases for b in criteria_eligible),'positive_source_contribution_flag_continuous_SCVs':len(contribution_eligible),'positive_source_contribution_flag_updated_SCVs':sum(b in update_bases for b in contribution_eligible)},
 'current_dated_assertion_contrast':{'all_stable_variants':len(stable),'all_statuses':dict(collections.Counter(c['status'] for c in contrasts)),'different_tier_all_pairs_different_SCV':sum(c['status']=='determinate_different_tier' and c['all_endpoint_pairs_different_SCV'] for c in contrasts),'different_tier_all_pairs_same_submitter':sum(c['status']=='determinate_different_tier' and c['all_endpoint_pairs_same_submitter'] for c in contrasts),'different_tier_all_pairs_same_raw_condition':sum(c['status']=='determinate_different_tier' and c['all_endpoint_pairs_same_raw_condition'] for c in contrasts),'different_tier_variants_with_observed_primary_one_year_update':sum(c['status']=='determinate_different_tier' and bool(c['primary_direct_update_events']) for c in contrasts),'different_tier_submitter_identity_categories':dict(collections.Counter(c['submitter_pair_identity'] for c in contrasts if c['status']=='determinate_different_tier')),'different_tier_condition_identity_categories':dict(collections.Counter(c['raw_condition_pair_identity'] for c in contrasts if c['status']=='determinate_different_tier')),'no_cross_window_rate_comparison':True},
 'latest_PLP_condition_scope':{'official_PLP_stable_variants':len(portability),'all_scope_counts':dict(collections.Counter(v['PLP_evidence_scope'] for v in portability)),'direct_ID_only_all_scope_counts':dict(collections.Counter(v['identifiers_only_PLP_evidence_scope'] for v in portability)),'missing_or_ambiguous_source_contributor_links':len(missing_contributors),'condition_is_not_allele_function_or_patient_diagnosis':True},
 'all35_gene_records_retained':len(per_gene),'all_new_statistical_association_tests':0,'satisfactory_manuscript_mainline':False,'independent_QA_complete':False,'limits':'One-year public archive census, current externally pinned vocabulary, no lifetime/patient/population rate, no inferred clinical harm or novel pathogenic gene. Missing condition annotation is unknown. Different estimands and prior-art coverage require review.'}
dump(O/'logs/descriptive_results_pending_independent_QA.json',result)
derived=[{'file':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(D.glob('*')) if p.is_file()]
dump(O/'logs/derived_source_and_result_manifest.json',{'utc':now(),'protocol_sha256':sha(P),'files':derived})
print('STAGE39_DESCRIPTIVE_FIXED_ENDPOINTS_COMPLETE_PENDING_INDEPENDENT_QA',json.dumps(result['primary'],ensure_ascii=False),flush=True)
