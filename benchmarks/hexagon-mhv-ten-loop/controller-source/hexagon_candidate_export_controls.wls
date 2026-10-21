(* Independent Wolfram WXF decoder supplies known QQ controls to the exporter.
   The same decoder also checks exported matrices against every TSV entry. *)
$HistoryLength=0;
args=Rest[$ScriptCommandLine];
If[Length[args]!=3,Print["usage: reference|verify ROOT VECTOR_TSV"];Exit[1]];
{mode,root,vector}=args;
columns=Rest[Import[FileNameJoin[{root,"columns.tsv"}],"TSV"]];
If[Length[columns]!=6,Print["invalid columns"];Exit[2]];
offset=0; flat={};
Do[
 {block,first,last,base,file}=row;
 If[base!=offset,Print["offset mismatch"];Exit[3]];
 tensor=Import[FileNameJoin[{root,file}],"WXF"];
 If[Dimensions[tensor]!={first,last},Print["dimensions mismatch"];Exit[4]];
 values=Flatten[Normal[tensor]];
 If[!AllTrue[values,MatchQ[#,_Integer|_Rational]&],Print["nonrational values"];Exit[5]];
 flat=Join[flat,values];offset+=first last;
 ,{row,columns}];
If[mode=="reference",
 If[FileExistsQ[vector],Print["refuse overwrite"];Exit[6]];
 Export[vector,MapIndexed[{First[#2]-1,ToString[Numerator[#1],InputForm],ToString[Denominator[#1],InputForm]}&,flat],"TSV"],
 If[mode!="verify",Exit[7]];
 lines=Select[Import[vector,"Lines"],!StringStartsQ[#,"#"]&];
 expected=Map[Function[line,With[{s=StringSplit[line,"\t"]},If[Length[s]!=3,Exit[8]];{FromDigits[s[[1]]],ToExpression[s[[2]],InputForm]/FromDigits[s[[3]]]}]],lines];
 If[expected[[All,1]]!=Range[0,Length[flat]-1]||expected[[All,2]]!=flat,Print["exact coefficient mismatch"];Exit[9]];
];
Print["INDEPENDENT_WOLFRAM_",ToUpperCase[mode],"_PASS columns=",Length[flat]," nnz=",Count[flat,x_/;x!=0]];
Exit[0];
