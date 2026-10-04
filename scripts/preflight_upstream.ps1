param(
    [Parameter(Mandatory = $true)][string]$ProjectRoot,
    [Parameter(Mandatory = $true)][string]$AuditDir,
    [string]$Rscript = 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe'
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $ProjectRoot).Path
if (Test-Path -LiteralPath $AuditDir) {
    throw "AuditDir already exists: $AuditDir"
}
$required = @(
    'data/gwas/GCST90691600.h.tsv',
    'data/annotations/gencode/gencode.v50.primary_assembly.annotation.gtf.gz',
    'reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pgen',
    'reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pvar.zst',
    'reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.psam',
    'data/expression/GTEx/Skin_Sun_Exposed_Lower_leg.v10.eQTLs.signif_pairs.parquet',
    'data/expression/GTEx/Skin_Not_Sun_Exposed_Suprapubic.v10.eQTLs.signif_pairs.parquet'
)
$records = foreach ($rel in $required) {
    $path = Join-Path $root ($rel -replace '/', '\')
    $item = Get-Item -LiteralPath $path -ErrorAction SilentlyContinue
    [pscustomobject]@{
        relative_path = $rel
        present = [bool]$item
        bytes = if ($item) { $item.Length } else { $null }
    }
}
$apa = Join-Path $root 'data\GTEx_v10_apaQTL'
$records += [pscustomobject]@{
    relative_path = 'data/GTEx_v10_apaQTL/'
    present = (Test-Path -LiteralPath $apa -PathType Container)
    bytes = $null
}
$dir = New-Item -ItemType Directory -Path $AuditDir
$records | Export-Csv -LiteralPath (Join-Path $dir.FullName 'input_presence.csv') -NoTypeInformation -Encoding UTF8
$tool = Join-Path $PSScriptRoot '..\tools\plink2\plink2.exe'
if (Test-Path -LiteralPath $Rscript) {
    & $Rscript --version 2>&1 | Out-File -LiteralPath (Join-Path $dir.FullName 'R_version.txt')
}
if (Test-Path -LiteralPath $tool) {
    & $tool --version 2>&1 | Out-File -LiteralPath (Join-Path $dir.FullName 'plink_version.txt')
}
$records | Format-Table -AutoSize
Write-Output "Preflight record: $($dir.FullName)"
