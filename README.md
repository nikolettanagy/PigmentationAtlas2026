
PigmentationAtlas

PigmentationAtlas is an R-based analysis workflow for prioritizing candidate effector genes at human pigmentation loci using GWAS data, skin molecular QTLs, fine-mapping, and Bayesian colocalization.

Version 1.0.0

The independently replayed and archived first release is available on Zenodo:

PigmentationAtlas v1.0.0 — DOI: 10.5281/zenodo.23134292](https://doi.org/10.5281/zenodo.23134292

Use the Zenodo archive when reproducing or citing the v1.0.0 analysis. It contains the validated scripts, `RUNBOOK.md`, `renv.lock`, a worked example, selected results and figures, quality-control reports, and SHA-256 records for 53 external inputs.

The fresh replay produced 1,325 PLINK clumps and 9,814 audited colocalization units with no failed units. Of these, 30 were classified as strong H4 and 54 as moderate H4. Gene reclassification was evaluable for 130 clumps, of which 111 were reclassified. The 111/130 fraction applies only to evaluable clumps.

Reproduction

Follow the RUNBOOK.md inside the Zenodo v1.0.0 archive for the validated execution order and quality-control gates. The workflow requires Windows, R 4.6.0, the package versions recorded in `renv.lock`, and a matching PLINK 2 executable. The required GWAS, GTEx, GENCODE, and 1000 Genomes source inputs are not bundled; their paths and SHA-256 checksums are documented in the archive. Running the supplementary scripts in numerical order alone is not the validated release procedure.

License and data

Original PigmentationAtlas software and documentation are licensed under the MIT License. PLINK 2, R packages, external datasets, and data-derived materials retain their respective terms; see `THIRD_PARTY_NOTICES.md` in the Zenodo archive.

Citation

Nagy, N. (2026). *PigmentationAtlas: An Integrative GWAS–apaQTL Colocalization Pipeline for Human Pigmentation Genetics* (v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23134292
