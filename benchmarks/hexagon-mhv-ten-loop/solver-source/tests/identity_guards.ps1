param([Parameter(Mandatory=$true)][string]$Workspace,[Parameter(Mandatory=$true)][string]$Source,[Parameter(Mandatory=$true)][string]$Destination)
$ErrorActionPreference='Stop'
$Workspace=(Resolve-Path $Workspace).Path
if(Test-Path $Destination){throw 'Use a fresh fixture directory'}
if(Get-ChildItem $Source -Recurse -Force -Attributes ReparsePoint){throw 'Source contains links'}
New-Item -ItemType Directory $Destination | Out-Null
$Destination=(Resolve-Path $Destination).Path
Copy-Item -LiteralPath $Source -Destination "$Destination/different_basis" -Recurse
$basis="$Destination/different_basis/basis_identity.tsv"
$lines=[IO.File]::ReadAllLines($basis)
$key=$lines[0].Split("`t")[0]
$lines[0]=$key+"`t"+('0'*64)
[IO.File]::WriteAllLines($basis,$lines)
# Simulate a completed object from another basis, not a damaged output file.
$receipt="$Destination/different_basis/qq_complete.tsv"
$rows=[IO.File]::ReadAllLines($receipt)
for($i=0;$i -lt $rows.Length;$i++){if($rows[$i].StartsWith("basis_identity.tsv`t")){$rows[$i]="basis_identity.tsv`t"+(Get-FileHash $basis).Hash.ToLowerInvariant()}}
[IO.File]::WriteAllLines($receipt,$rows)
$fec=Join-Path $Workspace 'hexagon/output/irrep_selection/k_recurrence/exact'
$lec=Join-Path $Workspace 'hexagon/output/irrLEC/physical_seed/NMHV'
& "$PSScriptRoot/../bin/irr_conditions.exe" NMHV 4 2 "$Destination/different_basis" "$Destination/rejected" --conditions mrk --samples 32 --fec-root $fec --lec-root $lec --data "$PSScriptRoot/../data" *> "$Destination/rejected.log"
if($LASTEXITCODE -ne 1 -or -not (Select-String "$Destination/rejected.log" -Pattern 'Source recurrence basis changed' -Quiet)){throw 'Mismatched basis was not rejected'}
$bulk=Join-Path $Workspace 'temp/irr_physical_benchmark_20260920/bulk/NMHV_4p2'
Copy-Item $bulk "$Destination/incomplete_boundary" -Recurse
$path="$Destination/incomplete_boundary/inputs.sha256"
$rows=[IO.File]::ReadAllLines($path)
$removed=$false
$filtered=foreach($line in $rows){if(-not $removed -and $line -match 'K_FEC_4_pe_A1'){$removed=$true;continue};$line}
[IO.File]::WriteAllLines($path,$filtered)
& "$PSScriptRoot/../bin/irr_convert.exe" shift NMHV 4 2 "$Destination/incomplete_boundary" "$Destination/bad_shift" --head 4 --fec-root $fec --lec-root $lec *> "$Destination/bad_shift.log"
if($LASTEXITCODE -ne 1 -or -not (Select-String "$Destination/bad_shift.log" -Pattern 'twelve distinct' -Quiet)){throw 'SHIFT accepted incomplete boundary provenance'}
Write-Output 'IDENTITY GUARDS PASS: different declared basis, incomplete legacy boundary, same-split SHIFT.'
