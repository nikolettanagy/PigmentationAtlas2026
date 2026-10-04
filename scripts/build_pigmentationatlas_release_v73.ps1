param(
    [Parameter(Mandatory=$true)][string]$PackageRoot,
    [Parameter(Mandatory=$true)][string]$Workspace,
    [Parameter(Mandatory=$true)][string]$OutputDirectory
)

$ErrorActionPreference = 'Stop'
$package = (Resolve-Path -LiteralPath $PackageRoot).Path
$workspace = (Resolve-Path -LiteralPath $Workspace).Path
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
$output = (Resolve-Path -LiteralPath $OutputDirectory).Path
$name = 'PigmentationAtlas_first_release_v1.0.0'
$stage = Join-Path $output $name
$zip = Join-Path $output ($name + '.zip')
if ((Test-Path -LiteralPath $stage) -or (Test-Path -LiteralPath $zip)) {
    throw "Release destination already exists; choose a new OutputDirectory: $output"
}

$required = @(
    'results/release_audit_v72/input_sha256.tsv',
    'results/release_audit_v72/renv.lock',
    'results/release_audit_v71/audit.tsv',
    'results/release_audit_v71/R_package_versions.tsv',
    'results/release_audit_v71/sessionInfo.txt',
    'results/release_figure_audit_v70/checks.tsv',
    'results/release_figure_audit_v70/figure_manifest.tsv',
    'results/module12/reconciliation_v59/audit.tsv',
    'results/module11/harmonization/catalogue/module11_3B_4_final_qc_status.tsv',
    'results/module12/colocalization/work_queue/module12_1_status.tsv',
    'results/module12/colocalization/catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz',
    'results/module12/colocalization/catalogue/summary/module12_3_coloc_catalogue_qc.tsv.gz',
    'results/module12/colocalization/catalogue/summary/module12_3_status.tsv',
    'results/module12/colocalization/prioritization/tables/module12_4_ranked_coloc_catalogue.tsv.gz',
    'results/module12/colocalization/prioritization/module12_4_status.tsv',
    'results/module12/colocalization/biological_annotation/tables/module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz',
    'results/module12/colocalization/biological_annotation/gene_level/module12_5C_gene_level_candidate_prioritization_GENCODEv50.tsv.gz',
    'results/module12/colocalization/biological_annotation/gene_level/module12_5C_best_gene_per_clump_GENCODEv50.tsv.gz',
    'results/module12/reclassification/module12_6_clump_gene_reclassification.tsv.gz',
    'results/module12/reclassification/module12_6_status.tsv'
)
foreach ($relative in $required) {
    $source = Join-Path $workspace ($relative.Replace('/','\'))
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Missing: $source" }
}
$audit = Import-Csv -LiteralPath (Join-Path $workspace 'results/release_audit_v71/audit.tsv') -Delimiter "`t"
$figureAudit = Import-Csv -LiteralPath (Join-Path $workspace 'results/release_figure_audit_v70/checks.tsv') -Delimiter "`t"
$figureFiles = Import-Csv -LiteralPath (Join-Path $workspace 'results/release_figure_audit_v70/figure_manifest.tsv') -Delimiter "`t"
if ($audit.Count -ne 37 -or @($audit | Where-Object result -ne 'PASS').Count -ne 0) { throw 'Release audit is not 37/37 PASS' }
if ($figureAudit.Count -ne 23 -or @($figureAudit | Where-Object status -ne 'PASS').Count -ne 0) { throw 'Figure numeric audit is not 23/23 PASS' }
if ($figureFiles.Count -ne 33 -or @($figureFiles | Where-Object status -ne 'PASS').Count -ne 0) { throw 'Figure file audit is not 33/33 PASS' }

New-Item -ItemType Directory -Path $stage | Out-Null
Get-ChildItem -LiteralPath $package -Force | Copy-Item -Destination $stage -Recurse -Force
foreach ($relative in $required) {
    $source = Join-Path $workspace ($relative.Replace('/','\'))
    $destination = Join-Path $stage ($relative.Replace('/','\'))
    New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
    Copy-Item -LiteralPath $source -Destination $destination
}
Copy-Item -LiteralPath (Join-Path $workspace 'results/release_audit_v72/renv.lock') -Destination (Join-Path $stage 'renv.lock')

$figureSource = Join-Path $workspace 'figures'
if (-not (Test-Path -LiteralPath $figureSource -PathType Container)) { throw "Missing: $figureSource" }
Copy-Item -LiteralPath $figureSource -Destination (Join-Path $stage 'figures') -Recurse

$note = @'
# PigmentationAtlas first release v1.0.0

This archive combines the validated v72 source distribution with selected
publication-level results from the fresh Windows R 4.6.0 replay (October 2026).
The 37-check release audit and 23-check / 33-file figure audit passed.

Run instructions are in RUNBOOK.md. The root renv.lock was restored to an
isolated R library; twelve primary package versions were independently checked.
Install Rtools45 on Windows when source packages must be compiled.

Raw GWAS, 1000 Genomes, GTEx expression, and GTEx apaQTL inputs are not in this
archive. results/release_audit_v72/input_sha256.tsv identifies the 53 checked
input files by SHA-256. Obtain source data under their applicable terms and
verify the hashes before running the full workflow.

The release contains the 9,814-row coloc summary and QC, ranked catalogue,
gene prioritization, reclassification table, and rendered figures. Per-unit
SNP-level files, dense LD matrices, and intermediate apaQTL chunks are not
included; they are regenerated by the workflow from the specified inputs.

The older 9,471-row catalogue is a distinct historical analysis. Its overlap
comparison is supplied as an audit only; it is not the release output.
The 111/130 reclassification fraction applies to evaluable clumps, not all
1,325 clumps. Archive DOI: pending deposit and publication.
'@
Set-Content -LiteralPath (Join-Path $stage 'RELEASE_NOTES_v1.0.0.md') -Value $note -Encoding UTF8

$stageFull = [System.IO.Path]::GetFullPath($stage).TrimEnd('\')
$manifest = Get-ChildItem -LiteralPath $stage -Recurse -File | ForEach-Object {
    [pscustomobject]@{
        relative_path = $_.FullName.Substring($stageFull.Length + 1).Replace('\','/')
        bytes = $_.Length
        sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    }
} | Sort-Object relative_path
$manifest | Export-Csv -LiteralPath (Join-Path $stage 'release_file_sha256.tsv') -Delimiter "`t" -NoTypeInformation -Encoding UTF8

Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.IO.Compression.ZipFile]::CreateFromDirectory($stage, $zip, [System.IO.Compression.CompressionLevel]::Optimal, $false)
$archive = [System.IO.Compression.ZipFile]::OpenRead($zip)
try {
    $entryCount = @($archive.Entries | Where-Object { $_.Name }).Count
    if ($entryCount -ne ($manifest.Count + 1)) { throw "ZIP entry count $entryCount differs from staged files $($manifest.Count + 1)" }
} finally { $archive.Dispose() }
$zipHash = (Get-FileHash -LiteralPath $zip -Algorithm SHA256).Hash.ToLowerInvariant()
$record = Join-Path $output ($name + '.sha256.txt')
Set-Content -LiteralPath $record -Value "$zipHash  $($name).zip" -Encoding ASCII
Write-Host "PASS release bundle: $zip"
Write-Host "Files: $entryCount; ZIP SHA256: $zipHash"
Write-Host "Hash record: $record"
