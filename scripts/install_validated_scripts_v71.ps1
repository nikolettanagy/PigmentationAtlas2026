param([Parameter(Mandatory=$true)][string]$Workspace)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $Workspace).Path
$package = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if ($root -eq $package) { throw 'Workspace and package must be distinct.' }
$map = [ordered]@{
  'Code S26_v53.R'='upstream\Code S26_v53.R'
  'Code S28_v55.R'='upstream\Code S28_v55.R'
  'Code S29A_v56.R'='upstream\Code S29A_v56.R'
  'Code S30_v57.R'='upstream\Code S30_v57.R'
  'Code S31_v58.R'='upstream\Code S31_v58.R'
  'Code S38_v63.R'='Code S38_v63.R'
  'Code S39_v64.R'='Code S39_v64.R'
  'Code S40_v67.R'='Code S40_v67.R'
  'Code S41_v69.R'='Code S41_v69.R'
  'Code S42_v69.R'='Code S42_v69.R'
  'audit_coloc_catalogue_v54.R'='audit_coloc_catalogue_v54.R'
  'audit_release_figures_v70.R'='audit_release_figures_v70.R'
  'reconcile_catalogues_v59.R'='reconcile_catalogues_v59.R'
  'explain_membership_v61.R'='explain_membership_v61.R'
  'reconcile_clump_labels_v62.R'='reconcile_clump_labels_v62.R'
}
$backup = Join-Path $root ('scripts\backup_pre_v71_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
foreach ($name in $map.Keys) {
  $source = Join-Path $PSScriptRoot ('validated\' + $name)
  $target = Join-Path (Join-Path $root 'scripts') $map[$name]
  if (-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Missing package file: $source" }
  New-Item -ItemType Directory -Force -Path (Split-Path $target -Parent) | Out-Null
  if (Test-Path -LiteralPath $target) {
    $saved = Join-Path $backup $map[$name]
    New-Item -ItemType Directory -Force -Path (Split-Path $saved -Parent) | Out-Null
    Copy-Item -LiteralPath $target -Destination $saved
  }
  Copy-Item -LiteralPath $source -Destination $target -Force
  $a=(Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
  $b=(Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
  if ($a -ne $b) { throw "Copy verification failed: $target" }
  Write-Output "INSTALLED $($map[$name]) $a"
}
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'validate_clean_release_v71.R') -Destination (Join-Path $root 'scripts\validate_clean_release_v71.R') -Force
Write-Output "Validated scripts installed in $root; previous workspace versions preserved in $backup"
