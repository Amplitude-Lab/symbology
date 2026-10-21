(* Read-only exact coefficient decoding; no reconstruction or physics fit. *)
$HistoryLength=0;
args=Rest[$ScriptCommandLine];
If[Length[args]!=1,Exit[1]];
jobs=Import[First[args],"RawJSON"];
Do[
 root=job["root"]; vector=job["vector"];
 If[FileExistsQ[vector],Print["refuse overwrite: ",vector];Exit[2]];
 columns=Rest[Import[FileNameJoin[{root,"columns.tsv"}],"TSV"]];
 If[Length[columns]!=6,Exit[3]];
 offset=0; flat={};
 Do[
  {block,first,last,base,file}=row;
  If[base!=offset,Exit[4]];
  If[first last==0,
   values={},
   tensor=Import[FileNameJoin[{root,file}],"WXF"];
   If[Dimensions[tensor]!={first,last},Exit[5]];
   values=Flatten[Normal[tensor]];
   If[!AllTrue[values,MatchQ[#,_Integer|_Rational]&],Exit[6]]
  ];
  flat=Join[flat,values];offset+=first last,
  {row,columns}];
 stream=OpenWrite[vector];
 MapIndexed[WriteString[stream,ToString[First[#2]-1],"\t",
   ToString[Numerator[#1],InputForm],"\t",
   ToString[Denominator[#1],InputForm],"\n"]&,flat];
 Close[stream];
 Print["EXACT_DECODE_PASS root=",root," columns=",Length[flat],
   " nonzero=",Count[flat,x_/;x!=0]],
 {job,jobs}];
Exit[0];
