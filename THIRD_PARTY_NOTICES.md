# License scope and third-party components

The root LICENSE (MIT) covers the original PigmentationAtlas software and
documentation authored by Nikoletta Nagy. It does not relicense third-party
tools, packages, reference data, source datasets, or material derived from
them. Consult the respective source providers' terms before redistributing
or reusing the included example slices and analysis outputs.

PLINK 2 is an external dependency, licensed GPLv3 or later by its authors.
Its executable is intentionally absent from this archive. Obtain the matching
Windows build from https://www.cog-genomics.org/plink/2.0/ and place plink2.exe
at tools/plink2/plink2.exe before running the PLINK-dependent stages.
PLINK source and licensing: https://www.cog-genomics.org/plink/2.0/dev
The validated Windows executable's SHA-256 was
247491bfca7512e070dc99d6565e9fc56f3a52ad5afc01286016271d34c4992f.
Verify the downloaded executable against this value for an exact replay;
the current download may be a later build.

The 53 large GWAS, GTEx, 1000 Genomes, and GENCODE input files are not
distributed here. Their SHA-256 values and paths appear in
results/release_audit_v72/input_sha256.tsv. R packages are installed via
renv.lock under their own licenses.
