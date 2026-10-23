(* Factorized four-loop affine input; reads only analytic/reconstructed lower loops. *)
$HistoryLength=0;base=DirectoryName[DirectoryName[$InputFileName]];r=base<>"/results";p=r<>"/ff_bootstrap/steinmann";out=r<>"/ff_four_loop";start=AbsoluteTime[];
If[!DirectoryQ[out],CreateDirectory[out]];
clean[a_]:=Module[{rr=Most[ArrayRules[SparseArray[a]]]},If[!AllTrue[Last/@rr,MatchQ[#,_Integer|_Rational]&],Print["nonrational tensor"];Exit[2]];SparseArray[Select[rr,Last[#]!=0&],Dimensions[a]]];zero[a_]:=Length[clean[a]["NonzeroValues"]]==0;assert[b_,s_]:=If[!TrueQ[b],Print["FAIL ",s];Exit[2]];
assert[zero[Import[r<>"/ff_freeze_w6/solution_w6.wxf"]-Import[p<>"/affine/solution_w6.wxf"]],"reconstructed F3 coefficient provenance"];
assert[zero[Import[r<>"/ff_freeze_regularity/steinmann_affine/origin.wxf"]-Import[p<>"/affine/solution_w4.wxf"]],"reconstructed F2 coefficient provenance"];
compact=Import[p<>"/affine/final_compact_w6.wxf"];f2=Import[p<>"/affine/reconstructed_w4.wxf"];f1=Import[r<>"/ff_boot_boundary_F41.wxf"];
project[a_,map_,weight_]:=Module[{z=a,d,n=Last[Dimensions[map]]},Do[d=Dimensions[z];z=ArrayReshape[Transpose[map].ArrayReshape[z,{88,Times@@Rest[d]}],Join[{n},Rest[d]]];z=Transpose[z,RotateRight[Range[weight]]],{weight}];clean[z]];
reports={};Do[
 tag=If[type=="DC","double","triple"];n=If[type=="DC",22,29];nf=If[type=="DC",6,9];mapPath=r<>If[type=="DC","/ff_boot_boundary_map.wxf","/ff_boot_triple_map.wxf"];col=p<>If[type=="DC","/collinear","/triple/collinear"];map=Import[mapPath];pr=Import[col<>"/first_w5.wxf"];t=clean[Transpose[pr].compact.map];Clear[pr];steps={};
 Do[b=Import[col<>"/first_w"<>ToString[k]<>"_basis.wxf"];d=Dimensions[b];assert[First[d]==First[Dimensions[t]]&&Last[d]==n,"projected prefix dimensions"];t=clean[ArrayReshape[Transpose[ArrayReshape[b,{d[[1]],d[[2]] n}]].t,{d[[2]],n^(7-k)}]];AppendTo[steps,<|"prefix"->k,"nonzero_words"->Length[t["NonzeroValues"]]|>];Clear[b],{k,5,2,-1}];
 e3=ArrayReshape[t,ConstantArray[n,6]];e1=project[f1,map,2];e2=project[f2,map,4];
 If[type=="TC",prior=Import[r<>"/ff_freeze_w6/triple_frame_expanded.wxf"][[1]];assert[zero[e3-prior],"fresh triple projection equals separately expanded solved affine origin"];Clear[prior]];
 factors=<|"E1"->e1,"F2"->e2,"F3"->e3|>;
 KeyValueMap[Export[out<>"/H8_"<>type<>"_"<>#1<>".wxf",#2,"WXF"]&,factors];
 j=Association@KeyValueMap[#1-><|"dimensions"->Dimensions[#2],"rules"->({First[#],ToString[Last[#],InputForm]}&/@Most[ArrayRules[#2]])|>&,factors];
 Export[out<>"/H8_"<>type<>"_lower_boundary.json",<|"boundary_size"->n,"finite_letters"->Range[nf],"factors"->j|>,"RawJSON"];
 AppendTo[reports,<|"boundary"->type,"alphabet_size"->n,"finite_letters"->nf,"factor_word_counts"->Association@KeyValueMap[#1->Length[#2["NonzeroValues"]]&,factors],"steps"->steps,"seconds_since_start"->AbsoluteTime[]-start,"max_memory_bytes"->MaxMemoryUsed[]|>];Print[Last[reports]];Clear[t,e1,e2,e3,factors,j],{type,{"DC","TC"}}];
sources={r<>"/ff_freeze_w6/solution_w6.wxf",r<>"/ff_freeze_regularity/steinmann_affine/origin.wxf",p<>"/affine/final_compact_w6.wxf",p<>"/affine/reconstructed_w4.wxf",r<>"/ff_boot_boundary_F41.wxf",r<>"/ff_boot_boundary_map.wxf",r<>"/ff_boot_triple_map.wxf"};
Export[out<>"/H8_projection_provenance.json",<|"reports"->reports,"source_hashes"->Association[Table[StringDrop[f,StringLength[base]+1]->FileHash[f,"SHA256","HexString"],{f,sources}]],"known_finite_boundary_targets_read"->False,"published_four_point_symbol_read"->False,"coefficient_provenance_checked"->True,"construction"->"Project saved independently reconstructed F3 compact prefix through DC/TC first_w5 map, then expand only lower-dimensional boundary prefix tensors. F1 analytic; F2 independently reconstructed. Three-loop solution equality checked before reuse of compact checkpoint.","seconds"->AbsoluteTime[]-start,"max_memory_bytes"->MaxMemoryUsed[]|>,"RawJSON"];
