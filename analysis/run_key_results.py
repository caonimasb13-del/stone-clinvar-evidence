"""Reproduce the frozen numerical core from verified selected-row checkpoints.

Standard library, Python3.11+. Never downloads, changes original inputs,
rejudges literature, publishes data, or writes to the project progress file.
"""
from pathlib import Path
import argparse, ast, copy, csv, datetime, gzip, hashlib, itertools, json, os, shutil, subprocess, sys, time, traceback

def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(p, o):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(o, ensure_ascii=False, indent=2) + '\n', encoding='utf8')

def bn(p):
    return str(p).replace('\\', '/').rsplit('/', 1)[-1]

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def jrows(p):
    with gzip.open(p, 'rt', encoding='utf8') as f:
        for line in f:
            yield json.loads(line)

def same_jsonl(a, b):
    count = 0
    for left, right in itertools.zip_longest(jrows(a), jrows(b)):
        assert left == right, ('Semantic JSON row mismatch', bn(a), count)
        count += 1
    return count

def table_rows(p):
    with Path(p).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.reader(f, delimiter=',' if str(p).endswith('.csv') else '\t'))

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--checkpoints-dir', required=True)
    ap.add_argument('--dictionary-dir', required=True)
    ap.add_argument('--reference-derived-dir', required=True, help='Preserved original core JSONL results for exact semantic comparison, not generation inputs')
    ap.add_argument('--reference-corrected', required=True)
    ap.add_argument('--workspace', required=True)
    ap.add_argument('--resources-dir', required=True, help='External fixed resource bundle; not included with the code export')
    args = ap.parse_args()
    resources = Path(args.resources_dir).resolve()
    m = Path(__file__).resolve().parent
    w = Path(args.workspace).resolve()
    src = Path(args.checkpoints_dir).resolve()
    dictionary = Path(args.dictionary_dir).resolve()
    assert not w.exists(), 'Use a new workspace. Never automatically repeat a completed step.'
    checkpoint = read(resources / 'frozen/checkpoint15_manifest.json')
    for r in checkpoint:
        p = src / r['file']
        assert p.stat().st_size == r['bytes'] and sha(p) == r['sha256'], p
    adaptation = read(m / 'ADAPTATION_MANIFEST.json')
    for r in adaptation:
        assert sha(m / r['adapted_file']) == r['adapted_sha256']
    w.mkdir(parents=True)
    for d in ['scripts', 'stage39/logs', 'stage39/metadata', 'stage39/tables', 'stage39/reviews', 'stage42/logs', 'stage42/metadata', 'stage42/tables', 'stage43/logs', 'stage43/metadata', 'stage43/tables', 'derived_stage39', 'literature_learning_27pdf_20261006', 'strengthening2_all35_landscape_source_influence_20261006', 'presentation_core/inputs', 'presentation_core/manuscript', 'final_tables', 'logs']:
        (w / d).mkdir(parents=True, exist_ok=True)
    state = {'started_utc': now(), 'status': 'PREPARING_FIXED_CHECKPOINT_REPLAY', 'pid': os.getpid(), 'scope': 'Selected-row checkpoint to numerical core, all1141chains and main source tables; raw15archive extraction not repeated', 'steps': [], 'new_models_or_inferential_tests': 0, 'network_allowed': False}
    dump(w / 'run_progress.json', state)
    timer = time.monotonic()

    def update():
        dump(w / 'run_progress.json', state)

    def execute(name):
        assert time.monotonic() - timer < 3600
        record = {'script': name, 'started_utc': now(), 'script_sha256': sha(w / 'scripts' / name), 'status': 'RUNNING'}
        state['steps'].append(record)
        state['status'] = 'RUNNING_' + name
        update()
        print('STEP_STARTED', name, flush=True)
        with (w / 'logs' / f"{len(state['steps']):02}_{name}.log").open('wb') as log:
            result = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(m / 'execute_guarded.py'), '--workspace', str(w), str(w / 'scripts' / name)], cwd=w, stdout=log, stderr=subprocess.STDOUT, timeout=1200)
        record.update(completed_utc=now(), exit_code=result.returncode, status='PASS' if result.returncode == 0 else 'FAILED')
        update()
        assert result.returncode == 0, ('Step failed; preserve outputs and inspect log', name)
        print('STEP_COMPLETE', name, flush=True)
    try:
        for r in adaptation:
            if r['adapted_file'].endswith('.py'):
                shutil.copyfile(m / r['adapted_file'], w / 'scripts' / Path(r['adapted_file']).name)
        for r in checkpoint:
            shutil.copyfile(src / r['file'], w / 'derived_stage39' / r['file'])
        original = read(resources / 'full_raw_core/frozen/scientific_protocol_original.json')
        pr = copy.deepcopy(original)
        for r in pr['inputs']:
            r['file'] = 'external_raw_sources/' + bn(r['file'])
        ext = read(resources / 'full_raw_core/frozen/external_dictionary_original.json')
        for r in ext['sources'].values():
            p = dictionary / bn(r['path'])
            assert p.is_file() and sha(p) == r['sha256']
            r['path'] = str(p)
        extfile = w / 'stage39/metadata/external_fixed35_gene_condition_phenotype_identity_only.json'
        dump(extfile, ext)
        panel = w / 'stage39/metadata/external_panel_configured.json'
        shutil.copyfile(resources / 'full_raw_core/frozen/external_panel_original.json', panel)
        assert sha(panel) == pr['external_panel']['sha256']
        pr['external_panel']['path'] = str(panel)
        pr['storage'] = str(w / 'derived_stage39')
        pr['external_dictionary'].update(path=str(extfile), sha256=sha(extfile))
        p = w / 'stage39/logs/protocol.json'
        dump(p, pr)
        (p.parent / 'protocol_sha256.txt').write_text(sha(p) + '\n', encoding='utf8')
        ex = read(resources / 'frozen/selected_extraction_receipt_verified_old_replay.json')
        ex['protocol_sha256'] = sha(p)
        for rs in [ex['endpoint_receipts'], ex['monthly_receipts']]:
            for r in rs:
                r['file'] = str(w / 'derived_stage39' / bn(r['file']))
        ex['checkpoint_relocation_note'] = 'Previously completed extraction receipt relocated by exact selected-file hash; no raw-source scan in this run'
        dump(w / 'stage39/logs/full_source_extraction_complete.json', ex)
        execute('census.py')
        execute('current_ID_sensitivity_QA.py')
        comparisons = []
        names = ['all_endpoint_variant_identity_and_exclusions.jsonl.gz', 'all_selected_SCV_eligibility_and_source_refs.jsonl.gz', 'all_primary_major_events.jsonl.gz', 'all_five_level_classification_events.jsonl.gz', 'all_adjacent_pair_major_events.jsonl.gz', 'all_current_dated_assertion_contrasts.jsonl.gz', 'all_latest_PLP_aggregate_condition_scope_and_sources.jsonl.gz', 'all_current_dated_contrast_complete_ID_sensitivity.jsonl.gz']
        for name in names:
            comparisons.append({'file': name, 'rows': same_jsonl(Path(args.reference_derived_dir) / name, w / 'derived_stage39' / name), 'all_scientific_rows_equal': True})
        dump(w / 'core8_jsonl_exact_comparison.json', comparisons)
        p42 = read(resources / 'full_raw_core/frozen/scope_repair_protocol_original.json')
        mapping = {'protocol.json': p, 'descriptive_results_pending_independent_QA.json': w / 'stage39/logs/descriptive_results_pending_independent_QA.json', 'all_latest_PLP_aggregate_condition_scope_and_sources.jsonl.gz': w / 'derived_stage39/all_latest_PLP_aggregate_condition_scope_and_sources.jsonl.gz', 'external_fixed35_gene_condition_phenotype_identity_only.json': extfile}
        historical_ledger_names = {'R2_ledger.json', 'R3_ledger.json'}
        dropped_ledger_inputs = [r for r in p42['inputs'] if bn(r['path']) in historical_ledger_names]
        assert {bn(r['path']) for r in dropped_ledger_inputs} == historical_ledger_names
        p42['inputs'] = [r for r in p42['inputs'] if bn(r['path']) not in historical_ledger_names]
        for r in p42['inputs']:
            path = mapping[bn(r['path'])]
            r.update(path=str(path), bytes=path.stat().st_size, sha256=sha(path))
        p = w / 'stage42/logs/protocol.json'
        dump(p, p42)
        (p.parent / 'protocol_sha256.txt').write_text(sha(p) + '\n', encoding='utf8')
        execute('scope_repair.py')
        execute('scope_repair_QA.py')
        corrected = w / 'derived_stage42/all4078_original_and_corrected_condition_scope_with_sources.jsonl.gz'
        comparisons.append({'file': corrected.name, 'rows': same_jsonl(args.reference_corrected, corrected), 'all_scientific_rows_equal': True})
        dump(w / 'core9_jsonl_exact_comparison.json', comparisons)
        p = w / 'stage43/logs/protocol.json'
        dump(p, read(resources / 'full_raw_core/frozen/release_and_use_task_protocol.json'))
        (p.parent / 'protocol_sha256.txt').write_text(sha(p) + '\n', encoding='utf8')
        repair = read(w / 'stage42/logs/repair_results_pending_independent_QA.json')
        dump(w / 'stage42/logs/scientific_review.json', {'output': repair['output'], 'scope': 'Output pointer required by the historical exporter, not a new independent scientific review'})
        execute('source_chains.py')
        source = w / 'derived_stage43/all1141_dated_contrast_interpretation_source_chains.jsonl.gz'
        qp = read(resources / 'frozen/transparency_protocol_original.json')
        hits = [r for r in qp['input_contract'] if bn(r['path']).endswith('48_all1141_dated_contrast_interpretation_source_chains.jsonl.gz')]
        assert len(hits) == 1
        hits[0].update(path=str(source), bytes=source.stat().st_size, sha256=sha(source))
        dump(w / 'literature_learning_27pdf_20261006/第一补强_当前差异透明度与疾病可比性_事后描述方案_20261006.json', qp)
        execute('transparency.py')
        inputs = w / 'presentation_core/inputs'
        sources = {'27_all35_gene_descriptive_census.tsv': w / 'stage39/tables/all35_gene_descriptive_census.tsv', '29_all_primary_submitter_concentration.tsv': w / 'stage39/tables/all_primary_submitter_concentration.tsv', '33_all35_corrected_scope_counts.tsv': w / 'stage42/tables/all35_corrected_scope_counts.tsv', '40_all_selected_SCV_eligibility_and_source_refs.jsonl.gz': w / 'derived_stage39/all_selected_SCV_eligibility_and_source_refs.jsonl.gz', '41_all_primary_major_events.jsonl.gz': w / 'derived_stage39/all_primary_major_events.jsonl.gz', '47_all4078_original_and_corrected_condition_scope_with_sources.jsonl.gz': corrected}
        for name, path in sources.items():
            assert path.exists(), path
            shutil.copyfile(path, inputs / name)
        static = resources / 'frozen/36_corrected_scope_categories_with_clinical_reading_limits.tsv'
        shutil.copyfile(static, inputs / static.name)
        sources[static.name] = inputs / static.name
        q = read(resources / 'frozen/source_deletion_protocol_original.json')
        s1 = w / 'strengthening1_current_contrast_transparency_20261006'
        for name in ['01_all1141_variant_level_transparency_and_comparability.tsv', '02_all_endpoint_source_records_with_provenance.tsv']:
            sources[name] = s1 / name
        for r in q['inputs']:
            path = sources[bn(r['path'])]
            r.update(path=str(path), bytes=path.stat().st_size, sha256=sha(path))
        out = w / 'strengthening2_all35_landscape_source_influence_20261006'
        qp = out / '第二项补强_全35基因与全部来源影响_检验前固定事后描述协议.json'
        dump(qp, q)
        guard_file = m / 'execute_guarded.py'
        dump(out / 'FREEZE_协议与输入核验.json', {'protocol_sha256': sha(qp), 'export_code_guard_path': str(guard_file), 'export_code_guard_sha256': sha(guard_file), 'binding': 'Supplied code file integrity guard outside the run workspace; no manuscript or reviewer content required; not author approval'})
        execute('source_deletion.py')
        census = read(w / 'stage39/logs/descriptive_results_pending_independent_QA.json')
        task = read(w / 'stage43/logs/record_interpretation_use_task_pending_QA.json')
        cases = read(w / 'stage43/metadata/deterministic_lowest_numeric_ID_task_stratum_examples.json')
        assert {r['VariationID'] for r in cases.values()} == {'913', '1434', '4597', '4598', '255914'}

        def table(name, rows):
            target = w / 'final_tables' / name
            with target.open('w', encoding='utf8', newline='') as f:
                wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
                wr.writeheader()
                wr.writerows(rows)
            return {'file': str(target), 'sha256': sha(target)}
        env = {'original': census, 'repair': repair, 'task': task, 'cases': cases, 'table': table, 'json': json, 'tables': []}
        exec(compile((m / 'algorithms/presentation_tables_verbatim.pyfrag').read_text(encoding='utf8'), 'verbatim_C152_table_generator', 'exec'), env)
        export = {'ST1_all1141_variant_transparency.tsv': s1 / '01_all1141_variant_level_transparency_and_comparability.tsv', 'ST2_all2320_endpoint_source_rows.tsv': s1 / '02_all_endpoint_source_records_with_provenance.tsv', 'ST3_all4704_fixed_cross_cells.tsv': s1 / '04_all_fixed_cross_cells_including_zeros.tsv', 'ST4_all_endpoint_marginals.tsv': s1 / '05_all_side_marginals_including_zeros.tsv', 'Aux_review_by_contribution.tsv': s1 / '06_review_by_contribution_full_cells.tsv', 'Aux_all_raw_status_values.tsv': s1 / '07_all_raw_status_and_contribution_values.tsv', 'ST5_all35_gene_landscape.tsv': out / '01_all35_gene_landscape_corrected_scope_and_separate_units.tsv', 'ST6_all265_historical_sources.tsv': out / '02_all265_historical_submitter_census.tsv', 'ST7_all35_gene_omissions.tsv': out / '03_all35_leave_one_gene_descriptors.tsv', 'ST8_all265_historical_source_omissions.tsv': out / '04_all265_leave_one_historical_submitter_descriptors.tsv', 'ST9_all96_current_source_dependencies.tsv': out / '05_all96_current_source_dependency_deletions.tsv', 'ST10_all1800_direction_cells.tsv': out / '06_all1800_gene_source_direction_cells_including_zeros.tsv', 'Aux_fixed_coverage_figure_labels.tsv': out / '07_fixed_coverage_based_figure_labels.tsv'}
        for name, path in export.items():
            shutil.copyfile(path, w / 'final_tables' / name)
        checks = []
        for expected in sorted((resources / 'expected_original26').iterdir()):
            made = w / 'final_tables' / expected.name
            if made.exists():
                a = table_rows(expected)
                b = table_rows(made)
                assert a == b, ('Exact table mismatch', expected.name)
                checks.append({'file': expected.name, 'status': 'REGENERATED_ALL_HEADERS_AND_CELLS_EXACT', 'rows': len(a) - 1, 'columns': len(a[0]), 'same_bytes': sha(expected) == sha(made), 'output_sha256': sha(made)})
            else:
                checks.append({'file': expected.name, 'status': 'PRESERVED_DOCUMENTARY_OR_MANUAL_INPUT_NOT_REGENERATED', 'source_sha256': sha(expected)})
        dump(w / 'all26_source_table_reproduction_status.json', checks)
        state.update(status='COMPLETE_NUMERICAL_CHECKPOINT_REPRODUCTION', completed_utc=now(), exact_core_JSONL=comparisons, table_status_counts=dict(__import__('collections').Counter((r['status'] for r in checks))), raw_source_extraction_reexecuted=False, later19_context_tables_replayed=False, literature_judgments_independently_revalidated=False, full_article_pipeline_replayed=False)
        update()
        dump(w / 'COMPLETE_checkpoint_reproduction.json', state)
        print(json.dumps({'status': state['status'], 'table_status_counts': state['table_status_counts'], 'core_semantic_rows': sum((r['rows'] for r in comparisons))}, ensure_ascii=False), flush=True)
    except Exception:
        state.update(status='FAILED_RETAINED_NO_AUTOMATIC_REPEAT', failed_utc=now(), traceback=traceback.format_exc())
        update()
        dump(w / 'FAILED_checkpoint_reproduction.json', state)
        raise
if __name__ == '__main__':
    main()
