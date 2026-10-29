"""Test file-integrity and workspace-write guards with tiny software fixtures only."""
from pathlib import Path
import argparse, ast, copy, hashlib, json, subprocess, sys

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--test-dir',required=True,help='New directory; fixture results are retained')
    args=ap.parse_args()
    root=Path(__file__).resolve().parent
    out=Path(args.test_dir).resolve()
    assert not out.exists(),'Use a new fixture directory; do not overwrite previous results.'
    out.mkdir(parents=True)
    tree=ast.parse((root/'algorithms/source_deletion.py').read_text(encoding='utf-8'))
    values=[v for n in ast.walk(tree) if isinstance(n,ast.Dict) for k,v in zip(n.keys,n.values) if isinstance(k,ast.Constant) and k.value=='export_code_guard_sha256_unchanged']
    assert len(values)==1
    expression=compile(ast.Expression(copy.deepcopy(values[0])),'actual_exported_guard_expression','eval')
    marker=out/'outside_workspace_fixture.txt'
    marker.write_text('Original fixture bytes\n',encoding='utf-8')
    frozen={'export_code_guard_path':str(marker),'export_code_guard_sha256':sha(marker)}
    namespace={'sha':sha,'Path':Path,'frozen':frozen}
    assert eval(expression,namespace) is True
    original=marker.read_bytes()
    marker.write_bytes(original+b'tampered fixture\n')
    assert eval(expression,namespace) is False
    marker.write_bytes(original)
    runs=[]
    for label,allowed in [('allowed_inside_workspace',True),('blocked_outside_workspace',False)]:
        w=out/label
        (w/'scripts').mkdir(parents=True)
        target=w/'allowed.txt' if allowed else marker
        p=w/'scripts'/'fixture.py'
        p.write_text('from pathlib import Path\nPath('+repr(str(target))+').write_text("fixture write",encoding="utf-8")\n',encoding='utf-8')
        command=[sys.executable,'-B','-X','utf8',str(root/'execute_guarded.py'),'--workspace',str(w),str(p)]
        r=subprocess.run(command,cwd=w,capture_output=True,text=True,encoding='utf-8',timeout=30)
        (w/'stdout.txt').write_text(r.stdout,encoding='utf-8')
        (w/'stderr.txt').write_text(r.stderr,encoding='utf-8')
        if allowed:
            assert r.returncode==0 and target.read_text(encoding='utf-8')=='fixture write'
        else:
            assert r.returncode!=0 and 'BLOCKED_WRITE_OUTSIDE_REPLAY_WORKSPACE' in r.stderr
            assert marker.read_bytes()==original
        runs.append({'test':label,'exit_code':r.returncode,'expected_outcome_observed':True})
    result={'status':'PASS_SOFTWARE_GUARD_FIXTURES_ONLY','unchanged_guard_positive':True,'tampered_guard_negative':True,'runs':runs,'clinical_data_used':False,'scientific_models_started':0,'full_pipeline_or_second_host_test':False}
    (out/'COMPLETE_guard_fixtures.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    main()
