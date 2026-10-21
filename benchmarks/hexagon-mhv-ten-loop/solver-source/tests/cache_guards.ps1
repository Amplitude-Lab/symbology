param([Parameter(Mandatory=$true)][string]$Source,[Parameter(Mandatory=$true)][string]$Destination)
$ErrorActionPreference='Stop'
$Source=(Resolve-Path -LiteralPath $Source).Path
if(Test-Path -LiteralPath $Destination){throw 'Use a fresh fixture directory'}
if(Get-ChildItem -LiteralPath $Source -Recurse -Force -Attributes ReparsePoint){throw 'Fixture source contains links'}
New-Item -ItemType Directory -Path $Destination | Out-Null
$Destination=(Resolve-Path -LiteralPath $Destination).Path
$exe=Join-Path $PSScriptRoot '../bin/irr_cache.exe'
foreach($name in @('valid','broken','active','linked','foreign','truncated')){Copy-Item -LiteralPath $Source -Destination (Join-Path $Destination $name) -Recurse}
$bad=Get-ChildItem (Join-Path $Destination 'broken') -Filter '*_origin.wxf' | Select-Object -First 1
[IO.File]::AppendAllText($bad.FullName,'corrupt fixture')
[IO.File]::WriteAllText((Join-Path $Destination 'foreign/cache/validation_parameters.txt'),'must not be deleted')
$short=Get-ChildItem (Join-Path $Destination 'truncated/cache') -Filter '*.bin' | Select-Object -First 1
[IO.File]::WriteAllText($short.FullName,'truncated')
$outside=Join-Path $Destination 'outside'
New-Item -ItemType Directory -Path $outside | Out-Null
[IO.File]::WriteAllText((Join-Path $outside 'must_survive.txt'),'sentinel')
New-Item -ItemType Junction -Path (Join-Path $Destination 'linked/cache/external') -Value $outside | Out-Null
$handle=[IO.File]::Open((Join-Path $Destination 'active/.active.lock'),[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
try {
    $before=@(Get-ChildItem (Join-Path $Destination 'valid/cache') -Filter '*.bin').Count
    & $exe $Destination
    if($LASTEXITCODE -or @(Get-ChildItem (Join-Path $Destination 'valid/cache') -Filter '*.bin').Count -ne $before){throw 'Dry run modified files'}
    & $exe $Destination --apply
    if($LASTEXITCODE){throw 'Cleanup failed'}
    if(@(Get-ChildItem (Join-Path $Destination 'valid/cache') -Filter '*.bin').Count){throw 'Completed cache not cleaned'}
    foreach($name in @('broken','active','linked','foreign','truncated')){if(-not @(Get-ChildItem (Join-Path $Destination "$name/cache") -Filter '*.bin').Count){throw "Unsafe deletion: $name"}}
    if(-not (Test-Path (Join-Path $outside 'must_survive.txt'))){throw 'Followed junction'}
} finally {$handle.Dispose()}
Write-Output 'CACHE GUARDS PASS: dry-run, complete, corrupt, active, junction, foreign, truncated; fixtures retained.'
