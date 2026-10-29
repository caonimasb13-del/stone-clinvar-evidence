"""Check exact external input identities without writing files or launching models."""
from pathlib import Path
import argparse,hashlib,json,sys
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def digest(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def exact(p,r):
 assert p.is_file() and p.stat().st_size==r['bytes'] and digest(p)==r['sha256'],str(p)
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 for name in ['resources-dir','checkpoints-dir','dictionary-dir','reference-derived-dir','reference-corrected']:ap.add_argument('--'+name,required=True)
 a=ap.parse_args();root=Path(__file__).resolve().parent
 manifest=load(root/'CODE_EXPORT_MANIFEST.json')
 for r in manifest['files']:exact(root/r['path'],r)
 resources=load(root/'external_resources.json')['resources']
 for r in resources:exact(Path(a.resources_dir)/r['path'],r)
 inputs=load(root/'input_identities.json')
 for r in inputs['selected_row_checkpoints']:exact(Path(a.checkpoints_dir)/r['file'],r)
 for r in inputs['dictionary_sources']:exact(Path(a.dictionary_dir)/r['filename'],r)
 refs=load(root/'reference_identities.json')['files']
 for r in refs:exact(Path(a.reference_corrected) if r['directory_argument']=='reference_corrected' else Path(a.reference_derived_dir)/r['file'],r)
 print(json.dumps({'status':'PASS_LOCAL_EXACT_INPUT_IDENTITIES_ONLY','external_resources':len(resources),'checkpoints':15,'dictionary_files':6,'reference_files':len(refs),'models_started':0,'workspace_created':False,'public_inputs_accessible':False},sort_keys=True))
if __name__=='__main__':main()
