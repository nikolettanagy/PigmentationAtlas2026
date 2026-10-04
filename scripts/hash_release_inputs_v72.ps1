param([Parameter(Mandatory=$true)][string]$Workspace,
      [Parameter(Mandatory=$true)][string]$SourceRoot,
      [Parameter(Mandatory=$true)][string]$AuditDir)
$ErrorActionPreference='Stop'
$w=(Resolve-Path -LiteralPath $Workspace).Path
$src=(Resolve-Path -LiteralPath $SourceRoot).Path
$dest=[IO.Path]::GetFullPath($AuditDir)
if ($dest -eq $src -or $dest.StartsWith($src+'\',[StringComparison]::OrdinalIgnoreCase)) { throw 'Audit output may not be within the original source.' }
$files=@('data\gwas\GCST90691600.h.tsv',
 'data\annotations\gencode\gencode.v50.primary_assembly.annotation.gtf.gz',
 'reference\1000G_GRCh38\1000G_GRCh38_EUR_uniqueID.pgen',
 'reference\1000G_GRCh38\1000G_GRCh38_EUR_uniqueID.pvar.zst',
 'reference\1000G_GRCh38\1000G_GRCh38_EUR_uniqueID.psam',
 'data\expression\GTEx\Skin_Sun_Exposed_Lower_leg.v10.eQTLs.signif_pairs.parquet',
 'data\expression\GTEx\Skin_Not_Sun_Exposed_Suprapubic.v10.eQTLs.signif_pairs.parquet')
$apa='data\GTEx_v10_apaQTL'
$sourceApa=@(Get-ChildItem -LiteralPath (Join-Path $src $apa) -Filter '*.parquet' -File)
if ($sourceApa.Count -ne 46) { throw "Expected 46 apaQTL Parquet files; got $($sourceApa.Count)" }
$files+=@($sourceApa | Sort-Object Name | ForEach-Object { Join-Path $apa $_.Name })
New-Item -ItemType Directory -Path $dest -Force | Out-Null
$manifest=Join-Path $dest 'input_sha256.tsv'
"relative_path`tbytes`tSHA256`tworkspace_sha256`tmatch" | Set-Content -LiteralPath $manifest -Encoding UTF8
foreach($rel in $files) {
  $a=Join-Path $src $rel; $b=Join-Path $w $rel
  if (-not (Test-Path -LiteralPath $a -PathType Leaf) -or -not (Test-Path -LiteralPath $b -PathType Leaf)) { throw "Missing input copy: $rel" }
  $size=(Get-Item -LiteralPath $a).Length
  if ($size -ne (Get-Item -LiteralPath $b).Length) { throw "Byte count differs: $rel" }
  $ha=(Get-FileHash -LiteralPath $a -Algorithm SHA256).Hash
  $hb=(Get-FileHash -LiteralPath $b -Algorithm SHA256).Hash
  $ok=$ha -eq $hb
  "$($rel.Replace('\','/'))`t$size`t$ha`t$hb`t$ok" | Add-Content -LiteralPath $manifest -Encoding UTF8
  Write-Output "$(if($ok){'PASS'}else{'FAIL'}) $rel"
  if (-not $ok) { throw "Input hash differs: $rel" }
}
Write-Output "PASS: $($files.Count) source/workspace SHA256 pairs; $manifest"
