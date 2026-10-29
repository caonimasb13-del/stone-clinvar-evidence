"""Reconstruct all current dated contrasts and complete the frozen ID sensitivity."""
from pathlib import Path
import json,gzip,hashlib,datetime,re,unicodedata,collections
R=Path(__file__).resolve().parents[1];O=R/'stage39';P=O/'logs/protocol.json'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def readj(p):
 with gzip.open(p,'rt',encoding='utf-8') as f:
  for line in f:yield json.loads(line)
assert sha(P)==(O/'logs/protocol_sha256.txt').read_text().strip();pro=json.loads(P.read_text(encoding='utf-8'));D=Path(pro['storage'])
res=json.loads((O/'logs/descriptive_results_pending_independent_QA.json').read_text(encoding='utf-8'))
identity=list(readj(D/'all_endpoint_variant_identity_and_exclusions.jsonl.gz'));stable={x['VariationID']:x['gene'] for x in identity if not x['flags']};assert len(stable)==res['stable_primary_single_gene_variants']
current=collections.defaultdict(dict)
for r in readj(D/'submission_summary_2026-10.txt.gz.SCV_closed.jsonl.gz'):
 if r['_SCV_base']:current[r['_SCV_base']][r['_source_line_sha256']]=r
variants=collections.defaultdict(list)
for b,records in current.items():
 if len(records)!=1:continue
 r=next(iter(records.values()))
 if r['VariationID'] in stable:variants[r['VariationID']].append(r)
def normalize(s):return re.sub(r'\s+',' ',unicodedata.normalize('NFKC',s)).strip().casefold()
mt={normalize(s):k for k,ss in pro['major_class_mapping'].items() for s in ss}
def day(s):
 # Independent explicit calendar parser, not the primary strptime function.
 s=s.strip();months={n:i+1 for i,n in enumerate(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])}
 try:
  if re.fullmatch(r'[A-Za-z]{3} \d{1,2}, \d{4}',s):mon,tail=s.split(' ',1);dd,yy=tail.split(', ');d=datetime.date(int(yy),months[mon.title()],int(dd))
  elif re.fullmatch(r'\d{4}-\d{1,2}-\d{1,2}',s):yy,mm,dd=map(int,s.split('-'));d=datetime.date(yy,mm,dd)
  elif re.fullmatch(r'\d{1,2}/\d{1,2}/\d{4}',s):mm,dd,yy=map(int,s.split('/'));d=datetime.date(yy,mm,dd)
  else:return None
  return d.isoformat() if d<=datetime.date(2026,10,1) else None
 except (ValueError,KeyError):return None
def context(r):return (normalize(r['SubmittedPhenotypeInfo']),normalize(r['ReportedPhenotypeInfo']))
def ids(r):
 parts=[p.strip() for p in re.split('[|;]',r['ReportedPhenotypeInfo']) if p.strip()]
 if not parts or any(not re.search(r'\bC[0-9]{7}\b',p,re.I) for p in parts):return None
 return frozenset(c.upper() for p in parts for c in re.findall(r'\bC[0-9]{7}\b',p,re.I))
def category(a,b,key):
 if not a or not b:return 'NOT_AVAILABLE'
 pairs=[(key(x),key(y)) for x in a for y in b]
 if any(x is None or y is None for x,y in pairs):return 'UNRESOLVED_COMPONENT'
 vals=[x==y for x,y in pairs]
 return 'ALL_SAME' if all(vals) else 'ALL_DIFFERENT' if not any(vals) else 'MIXED'
def unspecified(r):
 parts=[re.sub(r'^\s*(?:C[0-9]{7}|na)\s*:\s*','',s,flags=re.I) for s in re.split('[|;]',r['ReportedPhenotypeInfo'])]
 return all(normalize(s) in ['','not provided','not specified','not provided or not specified','na','-'] for s in parts)
expected={r['VariationID']:r for r in readj(D/'all_current_dated_assertion_contrasts.jsonl.gz')};sens=[];counter=collections.Counter()
for vid,gene in stable.items():
 rows=variants.get(vid,[]);valid=[(day(r['DateLastEvaluated']),r) for r in rows if day(r['DateLastEvaluated']) and mt.get(normalize(r['ClinicalSignificance']))];a=b=[];status='insufficient_dated_major_assertions'
 if len(valid)>=2:
  lo=min(d for d,r in valid);hi=max(d for d,r in valid);a=[r for d,r in valid if d==lo];b=[r for d,r in valid if d==hi];at={mt[normalize(r['ClinicalSignificance'])] for r in a};bt={mt[normalize(r['ClinicalSignificance'])] for r in b}
  status='all_dated_assertions_same_date' if lo==hi else 'mixed_major_tiers_at_tied_endpoint_date' if len(at)>1 or len(bt)>1 else 'determinate_same_tier' if at==bt else 'determinate_different_tier'
 e=expected[vid];assert status==e['status'];assert len(rows)==e['n_current_unique_SCVs'] and len(valid)==e['n_valid_dated_major_SCVs'];assert sorted(r['SCV'] for r in a)==sorted(e['earliest_SCVs']) and sorted(r['SCV'] for r in b)==sorted(e['latest_SCVs'])
 assert category(a,b,lambda r:r['_SCV_base'])==e['SCV_pair_identity'];assert category(a,b,lambda r:normalize(r['Submitter']))==e['submitter_pair_identity'];assert category(a,b,context)==e['raw_condition_pair_identity']
 r={'VariationID':vid,'gene':gene,'status':status,'complete_reported_MedGen_ID_pair_category':category(a,b,ids),'earliest_all_unspecified_conditions':bool(a) and all(unspecified(x) for x in a),'latest_all_unspecified_conditions':bool(b) and all(unspecified(x) for x in b),'earliest_MedGen_ID_sets':[sorted(ids(x)) if ids(x) is not None else None for x in a],'latest_MedGen_ID_sets':[sorted(ids(x)) if ids(x) is not None else None for x in b]};sens.append(r);counter[status]+=1
assert dict(counter)==res['current_dated_assertion_contrast']['all_statuses']
dest=D/'all_current_dated_contrast_complete_ID_sensitivity.jsonl.gz';assert not dest.exists()
with gzip.open(dest,'wt',encoding='utf-8',compresslevel=3,newline='\n') as f:
 for r in sens:f.write(json.dumps(r,ensure_ascii=False,separators=(',',':'))+'\n')
different=[r for r in sens if r['status']=='determinate_different_tier']
record={'utc':now(),'status':'PASS_ALL_STABLE_CURRENT_CONTRASTS_AND_FROZEN_CONDITION_ID_SENSITIVITY','protocol_sha256':sha(P),'script_sha256':sha(Path(__file__)),'all_variant_status_date_tie_and_identity_checks':len(stable),'all_fixed_status_counts':dict(counter),'predeclared_complete_MedGen_ID_sensitivity_for_different_tier_contrasts':dict(collections.Counter(r['complete_reported_MedGen_ID_pair_category'] for r in different)),'different_tier_earliest_all_unspecified':sum(r['earliest_all_unspecified_conditions'] for r in different),'different_tier_latest_all_unspecified':sum(r['latest_all_unspecified_conditions'] for r in different),'sensitivity_output':{'file':str(dest),'bytes':dest.stat().st_size,'sha256':sha(dest)},'no_primary_endpoint_selection_or_multiplicity_change':True,'scope':'Complete MedGen ID matching is identity metadata; neither nonidentical names nor different IDs establish allele causality or distinct clinical disease mechanisms. No different observation-window rate ratio.'}
dump(O/'logs/independent_current_contrast_and_frozen_ID_sensitivity_QA.json',record)
print('STAGE39_ALL_CURRENT_CONTRASTS_QA_PASS',len(stable),record['predeclared_complete_MedGen_ID_sensitivity_for_different_tier_contrasts'],flush=True)
