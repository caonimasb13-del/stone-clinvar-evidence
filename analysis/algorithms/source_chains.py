"""Source-grounded interpretation guard, all fixed cases, no clinical-error test."""
from pathlib import Path
import collections,csv,datetime,gzip,hashlib,json
R=Path(__file__).resolve().parents[1];O=R/'stage43'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def rows(p):
 with gzip.open(p,'rt',encoding='utf-8') as f:
  for line in f:yield json.loads(line)
def write(p,rs):
 with gzip.open(p,'wt',encoding='utf-8',compresslevel=3,newline='\n') as f:
  for r in rs:f.write(json.dumps(r,ensure_ascii=False,separators=(',',':'))+'\n')
p=O/'logs/protocol.json';pr=load(p);assert sha(p)==(O/'logs/protocol_sha256.txt').read_text().strip()
frozen=load(R/'stage39/logs/protocol.json');D=Path(frozen['storage']);dest=D.parent/'derived_stage43';dest.mkdir(exist_ok=True)
ex=load(R/'stage39/logs/full_source_extraction_complete.json');sens=load(R/'stage39/logs/independent_current_contrast_and_frozen_ID_sensitivity_QA.json');assert sha(sens['sensitivity_output']['file'])==sens['sensitivity_output']['sha256'];ii={r['VariationID']:r for r in rows(sens['sensitivity_output']['file'])}
cases={r['VariationID']:r for r in rows(D/'all_current_dated_assertion_contrasts.jsonl.gz') if r['status']=='determinate_different_tier'};assert len(cases)==1141
history=collections.defaultdict(list);by_current={};receipt_map={Path(i['file']).name:i for i in frozen['inputs']}
for rc in ex['monthly_receipts']:
 assert sha(rc['file'])==rc['sha256'];month=Path(rc['file']).name.split('.txt.gz')[0];original=receipt_map[month+'.txt.gz']
 for r in rows(rc['file']):
  if r['VariationID'] not in cases:continue
  q={k:r[k] for k in ['VariationID','ClinicalSignificance','DateLastEvaluated','SubmittedPhenotypeInfo','ReportedPhenotypeInfo','ReviewStatus','Submitter','SCV','ContributesToAggregateClassification','_SCV_base','_SCV_version','_source_line','_source_line_sha256']};q.update(snapshot=month,original_input_filename=month+'.txt.gz',original_input_sha256=original['sha256'],original_input_url=original['url']);history[r['VariationID']].append(q)
  if month.endswith('2026-10'):by_current[(r['VariationID'],r['SCV'])]=q
out=[];strata=collections.Counter();examples={};table=[]
for vid,c in sorted(cases.items(),key=lambda x:int(x[0])):
 earliest=[by_current[(vid,s)] for s in c['earliest_SCVs']];latest=[by_current[(vid,s)] for s in c['latest_SCVs']];assert earliest and latest
 assert all(a['_SCV_base']!=b['_SCV_base'] for a in earliest for b in latest)
 sensitivity=ii[vid];events=c['primary_direct_update_events'];key=sensitivity.get('complete_reported_ID_set_pair_identity',sensitivity.get('MedGen_ID_pair_identity',sensitivity.get('condition_ID_pair_identity')))
 # Discover only the already fixed independent sensitivity field, not a new endpoint.
 if key is None:
  available=[(k,v) for k,v in sensitivity.items() if isinstance(v,str) and v in {'ALL_SAME','ALL_DIFFERENT','UNRESOLVED_COMPONENT','MIXED'}];assert len(available)==1,(vid,sensitivity);key=available[0][1]
 stratum=key+'|observed_primary_update_'+('yes' if events else 'no');strata[stratum]+=1
 record={'VariationID':vid,'gene':c['gene'],'current_record_URL':'https://www.ncbi.nlm.nih.gov/clinvar/variation/'+vid+'/','snapshot_URL_scope':'Current web page may change; frozen archive URL/hash/source line below is the result provenance.','original_current_contrast':c,'original_fixed_complete_C_ID_sensitivity':sensitivity,'earliest_current_assertion_source_records':earliest,'latest_current_assertion_source_records':latest,'all13_public_snapshot_records_for_variant':history[vid],'all_public_present_SCV_bases_in_window':sorted({r['_SCV_base'] for r in history[vid]}),'observed_primary_major_update_events':events,'task_stratum':stratum,'permitted_label':'CURRENT_DATED_DIFFERENCE_ACROSS_SEPARATE_SCVS','same_SCV_update_from_these_current_endpoint_dates_alone_supported':False,'absence_of_lifetime_update_supported':False,'clinical_error_or_harm_supported':False,'interpretation':'Current endpoints belong to different submissions. Direct major updates are separately traceable in this observed window. Lack of a primary event is not lack of historical reclassification; excluded/missing/condition-changing records retained.'}
 out.append(record);examples.setdefault(stratum,record)
 table.append({'VariationID':vid,'gene':c['gene'],'earliest_evaluation_date':c['earliest_date'],'latest_evaluation_date':c['latest_date'],'earliest_SCVs':'|'.join(c['earliest_SCVs']),'latest_SCVs':'|'.join(c['latest_SCVs']),'raw_condition_identity':c['raw_condition_pair_identity'],'complete_C_ID_identity':key,'observed_primary_major_events':len(events),'same_SCV_update_from_current_dates_supported':False,'lifetime_no_update_claim_supported':False,'task_stratum':stratum})
target=dest/'all1141_dated_contrast_interpretation_source_chains.jsonl.gz';write(target,out);example_file=O/'metadata/deterministic_lowest_numeric_ID_task_stratum_examples.json';dump(example_file,examples)
with (O/'tables/all1141_record_label_interpretation.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(table[0]),delimiter='\t');w.writeheader();w.writerows(table)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'COMPLETE_BOUNDED_ALL1141_INTERPRETATION_TASK_SOURCE_EXPORT_PENDING_SEPARATE_QA','protocol_sha256':sha(p),'all1141_cases_retained':len(out),'all_current_endpoint_SCV_identity_different':True,'all_task_strata':dict(strata),'cases_with_observed_primary_update':sum(bool(x['observed_primary_major_update_events']) for x in out),'no_lifetime_rate_or_clinicalerror_comparison':True,'same_current_date_endpoint_pair_not_itself_longitudinal_SCV_trace':True,'main_cohort_output':{'path':str(target),'bytes':target.stat().st_size,'sha256':sha(target)},'all4078_scope_resource':load(R/'stage42/logs/scientific_review.json')['output'],'case_selection':'LowestnumericVariationIDforeachobservedmetadata stratum, fullcohortprovided; no best-case selection.','value_boundary':'This completeinterpretation use taskdemonstrates source tracing andcorrect label scope, not newgeneralmethod, clinicalbenefit, goldstandard performance orrefutationofhistoricrates. Finitepriorartincrement andportablecompletion remainrequired.','new_association_tests':0,'satisfactory_manuscript_mainline':False};dump(O/'logs/record_interpretation_use_task_pending_QA.json',summary)
print('STAGE43_ALL1141_COMPLETE_USE_TASK_SOURCE_EXPORT',json.dumps(strata),flush=True)
