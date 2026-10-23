(* Compact WXF adapter for generated exact H8 JSON factors. *)
$HistoryLength=0;base=DirectoryName[DirectoryName[$InputFileName]];out=base<>"/results/ff_four_loop";
Do[
 j=Import[out<>"/H8_"<>tag<>"_factors.json","RawJSON"];
 tensors=Association@KeyValueMap[Function[{name,a},name->SparseArray[(First[#]->ToExpression[Last[#]])&/@a["rules"],a["dimensions"]]],j["factors"]];
 If[!AllTrue[Join@@(#1["NonzeroValues"]&/@Values[tensors]),MatchQ[#,_Integer|_Rational]&],Print["FAIL rational coefficients"];Exit[2]];
 Export[out<>"/H8_"<>tag<>"_factors.wxf",<|"boundary"->tag,"boundary_size"->j["boundary_size"],"finite_letters"->j["finite_letters"],"factors"->tensors,"spurious_formula"->j["spurious_formula"],"wholly_finite_correction"->j["wholly_finite_correction"],"json_sha256"->FileHash[out<>"/H8_"<>tag<>"_factors.json","SHA256","HexString"]|>,"WXF"];
 Print[tag," ",Dimensions/@tensors," ",Length[#1["NonzeroValues"]]&/@tensors],{tag,{"DC","TC"}}];
