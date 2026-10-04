param([Parameter(Mandatory=$true)][string]$SourceRoot,
      [Parameter(Mandatory=$true)][string]$Workspace,
      [string]$Rscript='C:\Program Files\R\R-4.6.0\bin\Rscript.exe')
$ErrorActionPreference='Stop'
$source=(Resolve-Path -LiteralPath $SourceRoot).Path
$package=(Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) { throw "Rscript missing: $Rscript" }
if (Test-Path -LiteralPath $Workspace) { throw 'Fresh workspace path already exists; choose a new path.' }
& (Join-Path $PSScriptRoot 'prepare_clean_workspace.ps1') -SourceRoot $source -Workspace $Workspace
if (-not $?) { throw 'Workspace preparation failed.' }
$w=(Resolve-Path -LiteralPath $Workspace).Path
$logs=Join-Path $w 'logs\release_v72';New-Item -ItemType Directory -Force -Path $logs | Out-Null
Push-Location $w
try {
  & $Rscript (Join-Path $w 'scripts\run_upstream_phase1_v36.R') --workspace $w --rscript $Rscript
  if ($LASTEXITCODE -ne 0) { throw "Phase 1 failed; see $w\logs\phase1_v36" }
  & $Rscript (Join-Path $w 'scripts\audit_upstream_phase1_v37.R') $w
  if ($LASTEXITCODE -ne 0) { throw 'Phase 1 audit failed.' }
}
finally { Pop-Location }
Write-Output "PASS fresh Phase 1; workspace: $w"
