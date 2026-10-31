# Stone ClinVar evidence

Source-linked interpretation of public ClinVar records in a fixed nephrolithiasis gene panel. The public v1.0.0 release provides the three exact reviewed ZIP packages:

- [Complete bounded numerical reproduction package](https://github.com/caonimasb13-del/stone-clinvar-evidence/releases/download/v1.0.0/Stone_ClinVar_fixed_checkpoint_PORTABLE_PRIVATE.zip):110297439 bytes,122 files and66 runtime file identities.
- [Matching analysis code](https://github.com/caonimasb13-del/stone-clinvar-evidence/releases/download/v1.0.0/Stone_ClinVar_CODE_for_PORTABLE_INPUTS_PRIVATE.zip):105402 bytes,43 files.
- [Preserved derived resource](https://github.com/caonimasb13-del/stone-clinvar-evidence/releases/download/v1.0.0/Stone_ClinVar_derived_resource_LICENSED_PRIVATE.zip):18402823 bytes,66 files,21371 saved variant reports,45 source tables and5 figures.

All three files were downloaded without account credentials and matched their complete SHA256 hashes. Checksums and canonical links are in metadata/release_assets_v1.0.0.json and metadata/public_release_verification.json. Tagged source commit: b670a90c2861dc2407e4a489f47c8cfa327a326a. No DOI has been assigned.

## Run the packaged entry

Extract the complete reproduction ZIP. Use Python3.11 or later and its standard library:

    python -B run_bundle.py check
    python -B run_bundle.py run --workspace NEW_OUTPUT_DIRECTORY_OUTSIDE_THIS_BUNDLE

Use a fresh output folder outside the package. The launcher checks every manifest member and all66 input/reference identities and blocks network access. Three saved NCBI URLs are unversioned; preserve exact file hashes rather than substituting latest downloads. The input checker reports only local identity checks; actual public access was tested separately.

## Verified scope and provenance

One uninterrupted run from a new directory on the same Windows host matched nine JSONL outputs with116505 heterogeneous, potentially overlapping comparison rows and every header/cell of23 generated source tables. Three manual/documentary tables,19 later SX tables and the saved-report interface are preserved rather than regenerated. Eight process audit records had no forbidden outside-root reads or network events. This is a bounded same-host check, not a full article-pipeline replay, second-host test, biological replication or clinical validation.

Seven author-owned JSON protocol/receipt copies normalize74 historical location strings. All non-location content and original upstream hashes are preserved. The remaining59 runtime files and all scientific generating algorithms, numerical driver and write guard retain original bytes. Two analysis identity manifests declare the projected resource identities; metadata/reproduction/path_projection.json explains the transformation. The private original-location map is excluded. Use the paired portable code and input packages.

The PRIVATE suffix and preparation-state descriptions inside the ZIPs identify their byte-frozen pre-publication state. This public release establishes current access. The exact reviewed ZIPs were not repacked. Older dated metadata and Git commits retain historical status; current access is recorded in PUBLICATION_STATUS.json and metadata/public_release_verification.json.

Raw-source locations, fixed versions, hashes and source-to-result mappings are retained in the packages. Selected assertion checkpoints and third-party dictionaries have unchanged bytes. Saved records are distinct from independent functional truth or patient outcomes.

## Licences and disclosure

MIT applies to author-controlled software. CC BY4.0 applies to author-original results and documentation. Third-party components retain provider terms, attribution, versions and source-specific notices; see LICENSE_SCOPE.txt and THIRD_PARTY_NOTICES.md. Author confirmation of no separate additional source-use permission is not a new rights-holder grant.

AI contributed substantively to literature assessment, protocols, code, verification, interpretation, writing and figures. Human source-review history is documented separately. No manuscript was submitted or publication acceptance claimed.
