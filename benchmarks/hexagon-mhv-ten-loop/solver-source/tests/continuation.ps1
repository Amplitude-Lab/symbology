param([Parameter(Mandatory=$true)][string]$Workspace,[Parameter(Mandatory=$true)][string]$Destination)
$ErrorActionPreference='Stop'
$Workspace=(Resolve-Path $Workspace).Path
$exe=Join-Path $PSScriptRoot '../bin/irr_conditions.exe'
$fec=Join-Path $Workspace 'hexagon/output/irrep_selection/k_recurrence/exact'
$lec=Join-Path $Workspace 'hexagon/output/irrLEC/physical_seed/NMHV'
$bulk=Join-Path $Workspace 'temp/irr_physical_benchmark_20260920/bulk/NMHV_4p2'
$data=Join-Path $PSScriptRoot '../data'
if(Test-Path $Destination){throw 'Use a fresh test destination'}
New-Item -ItemType Directory $Destination | Out-Null
$Destination=(Resolve-Path $Destination).Path
function RunCase($name,$inputObject,$conditions,$extra=@()) {
    $parent=Join-Path $Destination $name
    & $exe NMHV 4 2 $inputObject $parent --conditions $conditions --samples 64 --fec-root $fec --lec-root $lec --data $data @extra *> (Join-Path $Destination "$name.log")
    if($LASTEXITCODE){throw "Failed case $name"}
    return $parent
}
$a=RunCase 'cache_only' $bulk 'collinear+mrk'
if(Test-Path "$a/collinear/qq_complete.tsv"){throw 'Intermediate was incorrectly exported'}
$b=RunCase 'resume' "$a/collinear" 'mrk'
$c=RunCase 'reverse' $bulk 'mrk+collinear'
$d=RunCase 'repeat' "$a/collinear_mrk" 'mrk'
# Same output directory, different local labels: resuming must not conflict.
$null=RunCase 'cache_only' "$a/collinear" 'mrk'
foreach($target in @("$b/collinear_mrk","$c/mrk_collinear","$d/collinear_mrk_mrk")) {
    foreach($file in Get-ChildItem "$a/collinear_mrk" -Filter '*_origin.wxf*') {
        if((Get-FileHash $file.FullName).Hash -ne (Get-FileHash (Join-Path $target $file.Name)).Hash){throw "Origin differs: $target/$($file.Name)"}
    }
}
$converter=Join-Path $PSScriptRoot '../bin/irr_convert.exe'
& $converter shift NMHV 4 2 "$a/collinear_mrk" "$Destination/shifted" --head 3 --fec-root $fec --lec-root $lec *> "$Destination/shift.log"
if($LASTEXITCODE){throw 'Affine shift failed'}
& $converter shift NMHV 3 3 "$Destination/shifted" "$Destination/roundtrip" --head 4 --fec-root $fec --lec-root $lec *> "$Destination/roundtrip.log"
if($LASTEXITCODE){throw 'Affine shift roundtrip failed'}
& $converter export NMHV 4 2 "$a/collinear_mrk" "$Destination/release" --head 3 --fec-root $fec --lec-root $lec *> "$Destination/release.log"
if($LASTEXITCODE){throw 'Affine release failed'}
foreach($file in Get-ChildItem "$a/collinear_mrk" -Filter '*_origin.wxf*') {
    if((Get-FileHash $file.FullName).Hash -ne (Get-FileHash (Join-Path "$Destination/roundtrip" $file.Name)).Hash){throw 'Shift changed origin'}
}
if((Get-Content "$Destination/shifted/condition_history.tsv") -notmatch 'collinear_mrk'){throw 'Shift lost condition history'}
Write-Output 'CONTINUATION PASS: cache-only, resume, reverse order, unique check, shift, release.'
