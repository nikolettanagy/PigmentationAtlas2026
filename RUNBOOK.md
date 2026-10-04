# Clean-run and release runbook

## Isolation and inputs

Never run computation in or overwrite the original `C:\Users\User\Desktop\PigmentationAtlas`. Create a distinct workspace with `scripts/prepare_clean_workspace.ps1 -SourceRoot ORIGINAL -Workspace NEW`; this copies GWAS, GENCODE GTF, and 1000G reference inputs. Copy the two GTEx skin significant-pair Parquet files and 46 skin apaQTL Parquet files into matching `data/` paths of the workspace. Confirm all 46 names and sizes against the source (a size comparison is a transfer check, not a cryptographic integrity check). Keep the original read-only and record SHA256 hashes before archive.

The workspace must have R 4.6.0, PLINK 2.0 a.7.1 (as used for the validated run), and the packages listed in the validation report. Download PLINK 2 from https://www.cog-genomics.org/plink/2.0/ and place `plink2.exe` in `tools/plink2/` of the working copy; the binary is not distributed with this release. Check the build version against the recorded environment. Preserve `R_package_versions.tsv` and `sessionInfo.txt`. This candidate does not include controlled copies or distribution rights for the large third-party inputs.

## Execution order

Use `Push-Location $w` and set `$env:PIGMENTATIONATLAS_ROOT=$w` for scripts that derive paths from the working directory. Check `$LASTEXITCODE` after each command. The full LD/SuSiE/coloc phases can take many hours; retain their logs and never delete valid per-unit output during resumption.

| Stage | Script or runner | Gate |
|---|---|---|
| 1, S1–S6 | `scripts/run_upstream_phase1_v36.R --workspace $w --rscript $R`; `scripts/audit_upstream_phase1_v37.R $w` | Fixed-count audit PASS |
| 2, S7–S8 | `scripts/run_upstream_expression_v38.R --workspace $w --input-root $original --rscript $R` | 2 tissues, 1,983 summary rows |
| 3, S9–S13 | `scripts/run_upstream_loci_v39.R --workspace $w --rscript $R` | 194 loci, 1,325 clumps |
| 4, S14 | `scripts/upstream/Code S14_v40.R` in full mode; `scripts/recover_large_ld_v42.R --workspace $w` | 1,325 ready; verify oversized 0266/0267 matrices and maps |
| 5, S15–S17 | SuSiE v45 full, PSD retry v46, `scripts/run_susie_qc_v47.R --workspace $w` | 265 publication loci; retain warnings/exclusions |
| 6, S18–S25 | Module11 audit, registry, regions, harmonization preparation, catalogue, final QC | 46 readable Parquet, 516/516 chunks, 14,765,046 rows, PASS |
| 7, S26–S28 | `scripts/upstream/Code S26_v53.R`, S27 `scripts/module12_2_prepare_and_test_single_coloc.R`, S28 `scripts/upstream/Code S28_v55.R`; `scripts/audit_coloc_catalogue_v54.R --workspace $w` | 9,814 eligible, 9,814 unit outputs, zero failed, audit PASS |
| 8, S29–S32 | S29, S29A_v56, S30_v57, S31_v58, S32 in numerical order | 30/54/239 H4 classes, 130 evaluable, 111 reclassified |
| 9, figures | S33–S37, S38_v63–S42_v69, S43 | `scripts/audit_release_figures_v70.R --workspace $w`: 23 checks and 33 files PASS |
| 10, final | `scripts/validate_clean_release_v71.R --workspace $w` | all release checks PASS; save report and checksums |

Stages with known restart corner cases require inspection of their runner logs and status files. The clean run used `Code S28_v55.R` to reconstruct missing aggregate summaries without repeating per-unit coloc. The v46 PSD retry retained 337 `warning_nonconverged` records as labeled. Do not treat a shell exit code of zero alone as evidence that all downstream deliverables were written.

## Remaining conditions before calling it an archived first release

- Run the v71 gate on the clean Windows workspace and review any FAIL with the specific file involved.
- Freeze a versioned code package plus provenance: input and script SHA256 manifest, package versions, PLINK version, output audit, and any redistribution limits.
- Check the final manuscript tables/legends against the clean catalogue; explicitly describe 130/1,325 evaluability and the unresolved old/clean queue membership difference.
- Perform at least a small fresh smoke run from a newly created workspace and record its output; the full clean run has already been completed, but this packaged candidate itself has not been independently rerun.
- Deposit/tag the final immutable archive and record its identifier. No DOI or public archive exists merely by making this package.
