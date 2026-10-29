# Checkpoint numerical code export

This private code candidate separates executable code from the fixed source/configuration and comparison materials. The seven Python algorithm files, table-generator fragment and child write guard are copied byte-for-byte from the previously staged-tested local module. The orchestration driver changes only resource-root location; CODE_ADAPTATION_CONTRACT.json records the inverse-AST equality check. Algorithm code contains expected counts and classification rules; those constants are retained, not relabelled as new data or observations.

Use Python 3.11 or later, standard library only. Every input is supplied separately. No download occurs.

Read-only preflight:

    python -B analysis/check_inputs_only.py --resources-dir FIXED_RESOURCES --checkpoints-dir CHECKPOINTS --dictionary-dir DICTIONARIES --reference-derived-dir REFERENCES --reference-corrected CORRECTED_REFERENCE

The preflight hashes code and the exact declared inputs. It writes no workspace, launches no scientific models and does not establish public input access or numerical correctness of this newly exported driver.

Numerical entry, for a separately authorized computation with exact inputs:

    python -B analysis/run_key_results.py --resources-dir FIXED_RESOURCES --checkpoints-dir CHECKPOINTS --dictionary-dir DICTIONARIES --reference-derived-dir REFERENCES --reference-corrected CORRECTED_REFERENCE --workspace NEW_WORKSPACE

FIXED_RESOURCES must retain the directory structure and exact hashes in external_resources.json. CHECKPOINTS and DICTIONARIES are listed in input_identities.json. References are listed in reference_identities.json. Do not substitute latest sources or loosen hashes. Existing workspaces are not overwritten. The original full raw-source extraction entry is outside this export.

The required resource contents, result tables, raw source-row checkpoints and dictionaries are not uploaded with this code. One legacy resource is a protected manuscript copy used only for a no-write invariant, not numerical generation. That copy and internal reviewer ledgers are intentionally withheld. This creates a real public-reproduction access gap. A later removal/replacement of a non-scientific guard would need its own documented change validation; it is not silently performed here.

The earlier same-host staged run compared nine core files (116,505 parsed rows across overlapping units) and 23 numerical source tables exactly. Three documentary/manual original tables and all19 later SX tables were not regenerated then. That earlier computation is not a new end-to-end execution of this relocated-resource driver, a second-machine test or biological validation. Prior failures, negatives, unknowns and stopped families remain applicable. No open licence, public release or DOI is assigned.
