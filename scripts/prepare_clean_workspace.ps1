param(
    [Parameter(Mandatory = $true)][string]$SourceRoot,
    [Parameter(Mandatory = $true)][string]$Workspace
)
$ErrorActionPreference = 'Stop'
$source = (Resolve-Path -LiteralPath $SourceRoot).Path
$package = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$workspaceFull = [System.IO.Path]::GetFullPath($Workspace)
if (Test-Path -LiteralPath $workspaceFull) { throw "Workspace already exists: $workspaceFull" }
if ($workspaceFull.StartsWith($source + '\', [System.StringComparison]::OrdinalIgnoreCase) -or
    $workspaceFull -eq $source) { throw 'Workspace cannot be inside the source project.' }
$inputs = @(
    'data\gwas\GCST90691600.h.tsv',
    'data\annotations\gencode\gencode.v50.primary_assembly.annotation.gtf.gz',
    'reference\1000G_GRCh38\1000G_GRCh38_EUR_uniqueID.pgen',
    'reference\1000G_GRCh38\1000G_GRCh38_EUR_uniqueID.pvar.zst',
    'reference\1000G_GRCh38\1000G_GRCh38_EUR_uniqueID.psam'
)
foreach ($rel in $inputs) {
    if (-not (Test-Path -LiteralPath (Join-Path $source $rel) -PathType Leaf)) {
        throw "Missing source input: $rel"
    }
}
New-Item -ItemType Directory -Path $workspaceFull | Out-Null
foreach ($name in @('scripts','config','tools')) {
    Copy-Item -LiteralPath (Join-Path $package $name) -Destination $workspaceFull -Recurse
}
foreach ($rel in $inputs) {
    $target = Join-Path $workspaceFull $rel
    New-Item -ItemType Directory -Path (Split-Path $target -Parent) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $source $rel) -Destination $target
    if ((Get-Item -LiteralPath $target).Length -ne (Get-Item -LiteralPath (Join-Path $source $rel)).Length) {
        throw "Input copy size mismatch: $rel"
    }
    Write-Output "COPIED $rel"
}
New-Item -ItemType Directory -Path (Join-Path $workspaceFull 'data\annotations\processed') -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $workspaceFull 'results') -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $workspaceFull 'logs') -Force | Out-Null
Write-Output "Workspace prepared: $workspaceFull"
Write-Output 'Input copies and computed results are separate from the source project.'
