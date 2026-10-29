"""Apply the frozen post-hoc descriptive protocol to existing verified source chains."""
import collections
import csv
import datetime
import gzip
import hashlib
import itertools
import json
import pathlib
import traceback

P = pathlib.Path(__file__).resolve().parents[1]
N = P/'literature_learning_27pdf_20261006'
O = P/'strengthening1_current_contrast_transparency_20261006'
PROTOCOL = N/'第一补强_当前差异透明度与疾病可比性_事后描述方案_20261006.json'
source = P/'derived_stage43/all1141_dated_contrast_interpretation_source_chains.jsonl.gz'
fixed_cases = ['913','1434','4597','4598','255914']
review_categories = ['CRITERIA_PROVIDED','NO_CRITERIA_PROVIDED','NO_CLASSIFICATION','EXPERT_PANEL','PRACTICE_GUIDELINE','UNKNOWN','MIXED']
contribution_categories = ['YES','NO','UNKNOWN','MIXED']
condition_categories = ['same','different','unresolved']
condition_map = {'ALL_SAME':'same','ALL_DIFFERENT':'different','UNRESOLVED_COMPONENT':'unresolved'}
review_map = {
 'criteria provided, single submitter':'CRITERIA_PROVIDED',
 'criteria provided, conflicting classifications':'CRITERIA_PROVIDED',
 'criteria provided, multiple submitters, no conflicts':'CRITERIA_PROVIDED',
 'no assertion criteria provided':'NO_CRITERIA_PROVIDED',
 'no classification provided':'NO_CLASSIFICATION',
 'no classification for the individual variant':'NO_CLASSIFICATION',
 'reviewed by expert panel':'EXPERT_PANEL', 'practice guideline':'PRACTICE_GUIDELINE'}

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(path,obj): path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def write_table(path,rows,fields=None):
    if fields is None: fields=list(rows[0])
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
def js(value): return json.dumps(value,ensure_ascii=False,separators=(',',':'))
def norm(value): return value.strip().lower() if isinstance(value,str) else None
def one_or_mixed(values):
    values=set(values)
    return sorted(values)[0] if len(values)==1 else ('MIXED' if len(values)>1 else 'UNKNOWN')

if O.exists(): raise RuntimeError('This step already has an output directory; check completion/old process instead of restarting.')
O.mkdir()
dump(O/'START_固定范围_20261006.json',{'started_utc':stamp(),'protocol':str(PROTOCOL),'protocol_sha256':sha(PROTOCOL), 'source':str(source),'new_tests':0,'n_expected_variants':1141,'candidate_selection':'all fixed current contrasts; no filtering'})

try:
    protocol=json.loads(PROTOCOL.read_text(encoding='utf-8'))
    input_contract=[x for x in protocol['input_contract'] if x['path']==str(source)]
    assert len(input_contract)==1 and sha(source)==input_contract[0]['sha256']
    with gzip.open(source,'rt',encoding='utf-8') as f: data=[json.loads(line) for line in f]
    assert len(data)==1141 and len({str(x['VariationID']) for x in data})==1141
    assert set(fixed_cases).issubset({str(x['VariationID']) for x in data})
    main=[]; side_rows=[]; source_objects=[]; events_out=[]; unresolved=[]
    marginal={s:{'review':collections.Counter(),'contribution':collections.Counter()} for s in ['earliest','latest']}
    cell_count=collections.Counter(); observed_conditions=collections.Counter(); event_conditions=collections.Counter()
    raw_review=collections.Counter(); raw_contribution=collections.Counter()
    identity_counts=collections.Counter(); both_category=collections.Counter(); cross_review_contrib=collections.Counter()
    for x in data:
        vid=str(x['VariationID']); original=x['original_current_contrast']
        assert original['status']=='determinate_different_tier' and str(original['VariationID'])==vid
        ci=x['original_fixed_complete_C_ID_sensitivity']
        assert str(ci['VariationID'])==vid
        condition=condition_map[ci['complete_reported_MedGen_ID_pair_category']]
        observed_conditions[condition]+=1
        row={'VariationID':vid,'gene':x['gene'],'current_record_URL':x['current_record_URL'],
             'earliest_evaluation_date':original['earliest_date'],'latest_evaluation_date':original['latest_date'],
             'earliest_original_tiers':js(original['earliest_tiers']),'latest_original_tiers':js(original['latest_tiers']),
             'complete_reported_C_ID_set_comparability':condition,
             'original_condition_category':ci['complete_reported_MedGen_ID_pair_category'],
             'clinical_error_supported':False, 'same_SCV_update_from_current_endpoint_dates_supported':False}
        side_bases={}
        for side in ['earliest','latest']:
            records=x[f'{side}_current_assertion_source_records']
            assert collections.Counter(z['SCV'] for z in records)==collections.Counter(original[f'{side}_SCVs'])
            mapped=[]; contributions=[]; raw_values=[]; contrib_raw_values=[]
            side_bases[side]=set()
            for index,z in enumerate(records):
                assert str(z['VariationID'])==vid
                assert z['_SCV_base']==z['SCV'].rsplit('.',1)[0]
                assert int(z['_SCV_version'])==int(z['SCV'].rsplit('.',1)[1])
                assert all(z[k] for k in ['_source_line_sha256','original_input_sha256','original_input_url'])
                review=z.get('ReviewStatus'); contribution=z.get('ContributesToAggregateClassification')
                rc=review_map.get(norm(review),'UNKNOWN')
                cc={'yes':'YES','no':'NO'}.get(norm(contribution),'UNKNOWN')
                mapped.append(rc);contributions.append(cc);raw_values.append(review);contrib_raw_values.append(contribution)
                raw_review[(side,js(review))]+=1;raw_contribution[(side,js(contribution))]+=1
                side_bases[side].add(z['_SCV_base'])
                source_row={'VariationID':vid,'gene':x['gene'],'side':side,'source_order_0based':index,
                    'SCV':z['SCV'],'SCV_base':z['_SCV_base'],'SCV_version':z['_SCV_version'],
                    'ClinicalSignificance_raw':z.get('ClinicalSignificance'),'DateLastEvaluated_raw':z.get('DateLastEvaluated'),
                    'Submitter_raw':z.get('Submitter'),'ReviewStatus_raw':review,'review_fixed_class':rc,
                    'ContributesToAggregateClassification_raw':contribution,'contribution_fixed_class':cc,
                    'SubmittedPhenotypeInfo_raw':z.get('SubmittedPhenotypeInfo'),'ReportedPhenotypeInfo_raw':z.get('ReportedPhenotypeInfo'),
                    'source_snapshot':z.get('snapshot'),'source_line':z['_source_line'],
                    'source_line_sha256':z['_source_line_sha256'],'source_filename':z['original_input_filename'],
                    'source_file_sha256':z['original_input_sha256'],'source_URL':z['original_input_url']}
                side_rows.append(source_row)
                source_objects.append({'VariationID':vid,'gene':x['gene'],'side':side,'source_order_0based':index,'review_fixed_class':rc,'contribution_fixed_class':cc,'original_source_record':z})
                if rc=='UNKNOWN' or cc=='UNKNOWN': unresolved.append({'VariationID':vid,'side':side,'SCV':z['SCV'],'unknown_review':rc=='UNKNOWN','unknown_contribution':cc=='UNKNOWN','original_source_record':z})
            rsummary=one_or_mixed(mapped); csummary=one_or_mixed(contributions)
            row.update({f'{side}_source_records':len(records),f'{side}_SCVs':js([z['SCV'] for z in records]),
                        f'{side}_review_summary':rsummary,f'{side}_review_class_set':js(sorted(set(mapped))),
                        f'{side}_review_raw_list':js(raw_values),f'{side}_contains_unknown_review':not mapped or 'UNKNOWN' in mapped,
                        f'{side}_contribution_summary':csummary,f'{side}_contribution_class_set':js(sorted(set(contributions))),
                        f'{side}_contribution_raw_list':js(contrib_raw_values),f'{side}_contains_unknown_contribution':not contributions or 'UNKNOWN' in contributions})
            marginal[side]['review'][rsummary]+=1;marginal[side]['contribution'][csummary]+=1
            cross_review_contrib[(side,rsummary,csummary)]+=1
        evs=x['observed_primary_major_update_events']
        row['has_original_strict_window_update']=bool(evs);row['original_strict_window_event_count']=len(evs)
        row['original_updated_SCV_bases']=js(sorted({e['SCV_base'] for e in evs}))
        if evs:event_conditions[condition]+=1
        relations=[]
        for event in evs:
            assert str(event['VariationID'])==vid and int(event['right_version'])>int(event['left_version'])
            base=event['SCV_base']
            first=base in side_bases['earliest'];last=base in side_bases['latest']
            relation=('BOTH_ENDPOINT_SCV_BASE_SETS' if first and last else 'EARLIEST_ENDPOINT_SCV_BASE_SET' if first else 'LATEST_ENDPOINT_SCV_BASE_SET' if last else 'NEITHER_ENDPOINT_SCV_BASE_SET')
            relations.append(relation);identity_counts[relation]+=1
            events_out.append({'VariationID':vid,'condition_comparability':condition,'source_identity_relation':relation,'original_strict_event':event})
        row['event_endpoint_SCV_identity_relations']=js(sorted(set(relations)))
        cell_count[(row['earliest_review_summary'],row['latest_review_summary'],row['earliest_contribution_summary'],row['latest_contribution_summary'],condition,'YES' if evs else 'NO')]+=1
        both_category[(row['earliest_review_summary'],row['latest_review_summary'])]+=1
        main.append(row)
    assert dict(observed_conditions)=={'unresolved':250,'different':659,'same':232}
    assert event_conditions['same']==0 and event_conditions['different']==15 and event_conditions['unresolved']==4
    assert sum(r['has_original_strict_window_update'] for r in main)==19
    for side in marginal:
        assert sum(marginal[side]['review'].values())==1141 and sum(marginal[side]['contribution'].values())==1141
    assert len(side_rows)==sum(r['earliest_source_records']+r['latest_source_records'] for r in main)
    assert all(len(json.loads(r['earliest_SCVs']))==r['earliest_source_records'] for r in main)
    write_table(O/'01_all1141_variant_level_transparency_and_comparability.tsv',main)
    write_table(O/'02_all_endpoint_source_records_with_provenance.tsv',side_rows)
    with (O/'03_all_endpoint_unchanged_source_records.jsonl').open('w',encoding='utf-8') as f:
        for r in source_objects:f.write(js(r)+'\n')
    full_cross=[]
    for er,lr,ec,lc,condition,event in itertools.product(review_categories,review_categories,contribution_categories,contribution_categories,condition_categories,['NO','YES']):
        full_cross.append({'earliest_review':er,'latest_review':lr,'earliest_contribution':ec,'latest_contribution':lc,'condition_comparability':condition,'has_observed_strict_window_event':event,'n_variants':cell_count[(er,lr,ec,lc,condition,event)]})
    assert sum(r['n_variants'] for r in full_cross)==1141
    write_table(O/'04_all_fixed_cross_cells_including_zeros.tsv',full_cross)
    marginals=[]
    for side in ['earliest','latest']:
        for axis,categories in [('review',review_categories),('contribution',contribution_categories)]:
            for category in categories:marginals.append({'side':side,'axis':axis,'category':category,'n_variants':marginal[side][axis][category],'denominator_variants':1141})
    write_table(O/'05_all_side_marginals_including_zeros.tsv',marginals)
    within=[]
    for side,rv,cv in itertools.product(['earliest','latest'],review_categories,contribution_categories):
        within.append({'side':side,'review':rv,'contribution':cv,'n_variants':cross_review_contrib[(side,rv,cv)],'denominator_variants':1141})
    write_table(O/'06_review_by_contribution_full_cells.tsv',within)
    raw_census=[]
    for axis,values in [('ReviewStatus',raw_review),('ContributesToAggregateClassification',raw_contribution)]:
        for (side,value),count in sorted(values.items()):raw_census.append({'axis':axis,'side':side,'raw_value_json':value,'n_source_records':count,'unit':'source-side row; not variants or patients'})
    write_table(O/'07_all_raw_status_and_contribution_values.tsv',raw_census)
    dump(O/'08_all_unknown_source_values.json',{'n_unknown_source_side_rows':len(unresolved),'all_unknown_rows':unresolved,'missing_endpoint_sides':[{'VariationID':r['VariationID'],'side':s} for r in main for s in ['earliest','latest'] if r[f'{s}_source_records']==0], 'mixed_variants':[r for r in main if any(r[f'{s}_{a}_summary']=='MIXED' for s in ['earliest','latest'] for a in ['review','contribution'])], 'all_unresolved_condition_variants':[r['VariationID'] for r in main if r['complete_reported_C_ID_set_comparability']=='unresolved']})
    dump(O/'09_all_original_events_with_endpoint_SCV_identity.json',{'events':events_out,'unit':'original strict event; do not count as independent patients','event_identity_counts':dict(identity_counts)})
    dump(O/'10_all5_fixed_cases_unchanged_selection.json',{'fixed_selection':'existing all5; no new case chosen','cases':[r for vid in fixed_cases for r in main if r['VariationID']==vid]})
    result={'created_utc':stamp(),'status':'COMPLETE_FIXED_POSTHOC_DESCRIPTIVE_STRENGTHENING1',
      'n_variants':1141,'n_endpoint_source_side_rows':len(side_rows),'n_earliest_source_side_rows':sum(r['earliest_source_records'] for r in main),'n_latest_source_side_rows':sum(r['latest_source_records'] for r in main),
      'side_marginals':{s:{axis:dict(counts) for axis,counts in axes.items()} for s,axes in marginal.items()},
      'both_side_review_combinations':[{'earliest':e,'latest':l,'n_variants':count} for (e,l),count in sorted(both_category.items())],
      'condition_counts':dict(observed_conditions),'variants_with_observed_strict_window_events_by_condition':{c:event_conditions[c] for c in condition_categories},
      'n_original_window_event_variants':19,'n_original_window_events_in_these_variants':len(events_out),
      'n_source_unknown_rows':len(unresolved),'n_tied_earliest_sides':sum(r['earliest_source_records']>1 for r in main),'n_tied_latest_sides':sum(r['latest_source_records']>1 for r in main),
      'n_full_zero_retaining_cross_cells':len(full_cross),'new_tests':0,'new_predictions':0,'clinical_truth_provided':False,
      'independent_biological_validation':False,'schema_mapping':['ALL_SAME -> same','ALL_DIFFERENT -> different','UNRESOLVED_COMPONENT -> unresolved'],
      'original_source_SCV_multiplicity_match_all1141':True,'source_raw_records_preserved':True,'original_condition_and_event_counts_reconstructed':True,
      'interpretation':'描述选定当前差异的来源构成；不是实验室错误、患者伤害、真正同病、全总体分类真值或终生更新率。',
      'source_SHA256':input_contract[0]['sha256'],'protocol_SHA256':sha(PROTOCOL),'script_SHA256':sha(pathlib.Path(__file__))}
    dump(O/'11_result_and_validation.json',result)
    outputs=[{'file':x.name,'bytes':x.stat().st_size,'sha256':sha(x)} for x in sorted(O.iterdir()) if x.is_file()]
    dump(O/'COMPLETE_输出完整性.json',{'finished_utc':stamp(),'exit':'success','outputs':outputs,'new_tests':0,'first_strengthening_only':True})
    print(json.dumps(result,ensure_ascii=False,indent=2))
except Exception as exc:
    dump(O/'FAILED_保留记录.json',{'failed_utc':stamp(),'error':repr(exc),'traceback':traceback.format_exc(),'do_not_relaunch_without_old_process_exit_proof':True})
    raise
