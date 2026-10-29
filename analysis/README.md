# Frozen numerical entry: external scientific inputs required

This private preparation export separates code from exact scientific inputs. It is not a self-contained public data deposit or a full article pipeline.

## Entry

Python 3.11 or later, standard library only. Supply a fresh workspace that does not exist.

```
python analysis/check_inputs_only.py --resources-dir RESOURCES --checkpoints-dir CHECKPOINTS --dictionary-dir DICTIONARIES --reference-derived-dir REFERENCE_CORE --reference-corrected CORRECTED_JSONL
python analysis/run_key_results.py --resources-dir RESOURCES --checkpoints-dir CHECKPOINTS --dictionary-dir DICTIONARIES --reference-derived-dir REFERENCE_CORE --reference-corrected CORRECTED_JSONL --workspace NEW_WORKSPACE
```

external_resources.json now lists 36 required configuration/comparison resources, 15 checkpoint and six dictionary identities are in input_identities.json, and nine comparison-file identities are in reference_identities.json. These contents remain external and their exact public delivery is unresolved.

## Declared non-numerical repair

The original two internal review ledgers were used only as generic input hash checks. The driver drops exactly those two checks and preserves all four scientific/configuration input entries and every historical rationale field in the Stage42 runtime protocol. The original local protocol and ledgers are preserved separately. No past human or AI review is reclassified.

The old protected manuscript served only as a no-write guard. The exported runtime instead hashes its supplied execute_guarded.py outside the run workspace and reports export_code_guard_sha256_unchanged. No manuscript or ledger content is now needed at runtime. The guarded executor still confines child writes to the new workspace. This is software integrity, not manuscript approval or biological validation.

Only that final guard expression and its two reporting-key literals change in source_deletion.py; inverse AST contracts preserve all numerical computation. Seven other algorithm files and execute_guarded.py retain their bytes. See CODE_ADAPTATION_CONTRACT.json and NON_NUMERICAL_GUARD_CHANGE.json.

## Verification scope

Targeted validation passed: one affected complete descriptive step reproduced all seven table outputs byte for byte; all 30 summary fields matched after excluding only timestamps, protocol relocation hash and the renamed guard status. The actual driver ledger filter retained all four scientific/configuration input objects and every other protocol field. Software fixtures passed the valid guard, tampered guard and actual outside-workspace write-denial checks. These are same-input software checks, not independent biological validation. Subsequently, one complete uninterrupted run of the unchanged exported checkpoint entry passed from a local ZIP containing all 66 exact runtime files. Nine core JSONL outputs (116,505 heterogeneous overlapping rows) and all 23 generated tables matched on the same Windows host. 8 main/child processes were audited; no old project/data-root reads were observed. This is a local numerical-entry test, not a full article pipeline or second-host validation. The tested input ZIP remains local HOLD; see ../metadata/packaged_checkpoint_replay_validation.json and ../metadata/input_access_routes.json. Prior same-host staged replay matched nine core JSONL files (116,505 heterogeneous overlapping rows) and 23 numerical source tables. Three original documentary/manual tables and all 19 later SX tables were not regenerated in that replay.

No source records, dictionary contents, result-table contents, protected manuscript or internal ledger contents are included. Private repository status does not resolve reuse terms or provide public data access. MIT is now assigned to author-controlled software and CC BY4.0 to author-original results/documentation, with upstream provider scopes retained; see ../LICENSE_SCOPE.txt. There is no public Release or DOI.

Run the small, data-free software guard checks in a new directory:

```
python analysis/verify_delivery_guards.py --test-dir NEW_FIXTURE_DIRECTORY
```
