# PigmentationAtlas first release v1.0.0

This archive contains the source code and selected results of the independently
replayed October 2026 release. The original README below describes the earlier
v71 candidate; the verified release results and limitations are stated in
RELEASE_NOTES_v1.0.0.md. Original PigmentationAtlas code is MIT licensed
(LICENSE); third-party materials have separate terms (THIRD_PARTY_NOTICES.md).

Before running the PLINK-dependent stages, install the PLINK 2 Windows
executable at `tools/plink2/plink2.exe` in the working copy. The archived
workflow used PLINK 2.0 a.7.1; see the provider and license information in
THIRD_PARTY_NOTICES.md. The binary is not included in this archive.

---

# PigmentationAtlas first release candidate v71

This candidate packages the code and validation gate for the independently rebuilt clean workspace. **It is not yet a published or independently archived release.** The original `PigmentationAtlas` project is input only; all computed outputs belong in a separate workspace. The validated clean run and the release package are distinct: packaging alone does not reproduce the results.

## Validated clean-run results (30 September 2026)

| Stage | Verified result |
|---|---:|
| GWAS input / harmonized variants | 21,874,448 / 21,831,407 |
| PLINK clumps / parent loci | 1,325 / 194 |
| Module 11 harmonized apaQTL files / rows | 516 / 14,765,046 |
| Module 12.1 total / eligible / excluded units | 9,915 / 9,814 / 101 |
| Module 12.3 audited unit summaries / failed units | 9,814 / 0 |
| Module 12.4 strong / moderate / suggestive H4 | 30 / 54 / 239 |
| Module 12.6 evaluable / reclassified / unchanged clumps | 130 / 111 / 19 |
| Figure audit | 23/23 table checks and 33/33 figure files |

The 101 exclusions have insufficient shared variants. SuSiE QC excluded 337 nonconverged, one failed, and 722 duplicate clumps from publication consideration; the 265 retained loci include warnings. Do not interpret the 111/130 reclassification rate as a rate among all 1,325 clumps.

The older independent catalogue has 9,471 rows. Of these, 7,609 share the full clump/tissue/phenotype/gene key with the clean catalogue and have identical posteriors and interpretations; 1,862 old-only and 2,205 clean-only records were absent from the other **full** work queue. Clump metadata checked for leading differing IDs retain the same lead variant across projects; a simple clump renumbering is not established. The detailed upstream membership cause remains open. Keep older and clean analyses separate in the manuscript.

## Run and verify

Use Windows PowerShell and R 4.6.0. [RUNBOOK.md](RUNBOOK.md) gives the validated order, source input handling, and commands. To install the corrected scripts into an existing **clean** workspace, preserving any previous workspace copies, run:

```powershell
$pkg = 'C:\Users\User\Desktop\PigmentationAtlas_first_release_candidate_v71\PigmentationAtlas_first_release_candidate'
$w = 'C:\Users\User\Desktop\PigmentationAtlas_clean_phase1_v33'
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$pkg\scripts\install_validated_scripts_v71.ps1" -Workspace $w
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' "$w\scripts\validate_clean_release_v71.R" --workspace $w
if ($LASTEXITCODE -ne 0) { throw 'Release validation failed' }
```

The validation writes `results/release_audit_v71/{audit.tsv,inputs.tsv,R_package_versions.tsv,sessionInfo.txt}` in the clean workspace. A PASS on the existing workspace confirms consistency and availability, not a fresh independent rerun. Check all gates before freezing an archive; record the exact package and input hashes and repository commit in the release record.

## v72 release gates

`scripts/run_fresh_phase1_v72.ps1` starts a new isolated workspace and validates Phase 1. `scripts/hash_release_inputs_v72.ps1` produces SHA-256 comparisons of the original and clean workspace copies of seven large foundational inputs and 46 apaQTL Parquet files. `scripts/lock_environment_v72.R` snapshots the observed Windows R 4.6.0 package versions and their recursive dependencies to an `renv.lock` in the clean workspace. These are new scripts requiring execution on the user's Windows host; no claim is made that the full workflow has yet run from v72.
