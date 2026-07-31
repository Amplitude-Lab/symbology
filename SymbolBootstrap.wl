(* ::Package:: *)

Get[If[$InputFileName=!="",FileNameJoin[{DirectoryName[$InputFileName],"SparseRREF","SparseRREF.wl"}],FileNameJoin[{Directory[],"SparseRREF","SparseRREF.wl"}]]]


BeginPackage["SymbolBootstrap`"]

DeclareAlphabet::usage="DeclareAlphabet[name, alphabet] registers a new symbol alphabet.";
ResetAlphabet::usage="ResetAlphabet[name, alphabet] overwrites a registered symbol alphabet and clears its conditions and cached tensors.";
ClearAlphabet::usage="ClearAlphabet[name] removes a registered symbol alphabet.";

SetAlphabetExpression::usage="SetAlphabetExpression[name, alphabetExpr] sets the parametrized expressions for an alphabet.";
SetClusterAdjacency::usage="SetClusterAdjacency[name, adjpairs] sets the cluster-adjacent ordered pairs for an alphabet.";
SetExtendedSteinmann::usage="SetExtendedSteinmann[name, nonadjpairs] sets the non-adjacent ordered pairs for extended Steinmann constraints.";
SetFirstEntry::usage="SetFirstEntry[name, firstentry] sets the first-entry letters for an alphabet.";
SetLastEntry::usage="SetLastEntry[name, lastentry] sets the last-entry letters for an alphabet.";
SetLetterTransformation::usage="SetLetterTransformation[name, transName, kineMap] sets one named letter transformation for an alphabet.";
SetAlphabetCondition::usage="SetAlphabetCondition[name, conditionName, content] sets one named alphabet condition.";

GetAlphabetConditionTensor::usage="GetAlphabetConditionTensor[name, conditionName] returns the tensor for a named alphabet condition. GetAlphabetConditionTensor[name, {conditionName1, conditionName2, ...}] combines dlogmat condition tensors.";
GetIntegrabilityTensor::usage="GetIntegrabilityTensor[name] returns the integrability condition tensor.";
GetClusterAdjacencyTensor::usage="GetClusterAdjacencyTensor[name] returns the cluster adjacency condition tensor.";
GetExtendedSteinmannTensor::usage="GetExtendedSteinmannTensor[name] returns the extended Steinmann condition tensor.";
GetFirstEntryTensor::usage="GetFirstEntryTensor[name] returns the first-entry seed tensor.";
GetLastEntryTensor::usage="GetLastEntryTensor[name] returns the last-entry seed tensor.";
GetLetterTransformationTensor::usage="GetLetterTransformationTensor[name, transName] returns one named letter transformation matrix.";

Begin["`Private`"]


$conditionKeys={"Alphabet","Expression","Extended Steinmann","Cluster Adjacency","First Entry","Last Entry","Letter Transformation"};
$resultKeys={"Integrability","Extended Steinmann","Cluster Adjacency","First Entry","Last Entry","Letter Transformation"};


DefaultCondition[alphabet_]:=<|"Alphabet"->alphabet,"Expression"->{},"Extended Steinmann"->{},"Cluster Adjacency"->{},"First Entry"->{},"Last Entry"->{},"Letter Transformation"-><||>|>
DefaultResult[]:=<|"Integrability"->{},"Extended Steinmann"->{},"Cluster Adjacency"->{},"First Entry"->{},"Last Entry"->{},"Letter Transformation"-><||>|>


AlphabetDeclaredQ[name_String]:=AssociationQ[$condition[name]]&&AssociationQ[$result[name]]


RequireAlphabet[name_String]:=If[AlphabetDeclaredQ[name],True,Print["Error: alphabet \"",name,"\" has not been declared"];False]


ConditionSetQ[name_String,key_String]:=KeyExistsQ[$condition[name],key]&&$condition[name][key]=!={}
ResultCachedQ[name_String,key_String]:=KeyExistsQ[$result[name],key]&&$result[name][key]=!={}


SetConditionValue[name_String,key_String,value_]:=($condition[name]=Join[$condition[name],<|key->value|>])
SetResultValue[name_String,key_String,value_]:=($result[name]=Join[$result[name],<|key->value|>])
ClearResult[name_String,key_String]:=If[AlphabetDeclaredQ[name]&&KeyExistsQ[$result[name],key],SetResultValue[name,key,{}]]


DeclareAlphabet[name_String,alphabet_?VectorQ]:=If[AlphabetDeclaredQ[name],Print["Error: alphabet \"",name,"\" is already declared; use ResetAlphabet to overwrite it"];Return[$Failed],
$condition[name]=DefaultCondition[alphabet];$result[name]=DefaultResult[];name]


ResetAlphabet[name_String,alphabet_?VectorQ]:=($condition[name]=DefaultCondition[alphabet];$result[name]=DefaultResult[];name)


ClearAlphabet[name_String]:=If[!AlphabetDeclaredQ[name],Print["Error: alphabet \"",name,"\" has not been declared"];Return[$Failed],
Unset[$condition[name]];Unset[$result[name]];name]


SetAlphabetExpression[name_String,alphabetExpr_?VectorQ]:=(If[!RequireAlphabet[name],Return[$Failed]];
SetConditionValue[name,"Expression",alphabetExpr];ClearResult[name,"Integrability"];SetResultValue[name,"Letter Transformation",<||>];name)


SetClusterAdjacency[name_String,adjpairs_?MatrixQ]:=(If[!RequireAlphabet[name],Return[$Failed]];
SetConditionValue[name,"Cluster Adjacency",adjpairs];ClearResult[name,"Cluster Adjacency"];name)


SetExtendedSteinmann[name_String,nonadjpairs_?MatrixQ]:=(If[!RequireAlphabet[name],Return[$Failed]];
SetConditionValue[name,"Extended Steinmann",nonadjpairs];ClearResult[name,"Extended Steinmann"];name)


SetFirstEntry[name_String,firstentry_?VectorQ]:=(If[!RequireAlphabet[name],Return[$Failed]];
SetConditionValue[name,"First Entry",firstentry];ClearResult[name,"First Entry"];name)


SetLastEntry[name_String,lastentry_?VectorQ]:=(If[!RequireAlphabet[name],Return[$Failed]];
SetConditionValue[name,"Last Entry",lastentry];ClearResult[name,"Last Entry"];name)


SetLetterTransformation[name_String,transName_String,kineMap_]:=(If[!RequireAlphabet[name],Return[$Failed]];
SetConditionValue[name,"Letter Transformation",Join[$condition[name]["Letter Transformation"],<|transName->kineMap|>]];
SetResultValue[name,"Letter Transformation",KeyDrop[$result[name]["Letter Transformation"],transName]];name)


SetAlphabetCondition[name_String,"Expression",content_?VectorQ]:=SetAlphabetExpression[name,content]
SetAlphabetCondition[name_String,"Cluster Adjacency",content_?MatrixQ]:=SetClusterAdjacency[name,content]
SetAlphabetCondition[name_String,"Extended Steinmann",content_?MatrixQ]:=SetExtendedSteinmann[name,content]
SetAlphabetCondition[name_String,"First Entry",content_?VectorQ]:=SetFirstEntry[name,content]
SetAlphabetCondition[name_String,"Last Entry",content_?VectorQ]:=SetLastEntry[name,content]
SetAlphabetCondition[name_String,{"Letter Transformation",transName_String},content_]:=SetLetterTransformation[name,transName,content]
SetAlphabetCondition[name_String,"Letter Transformation",content_]:=(Print["Error: letter transformation condition requires a transformation name; use SetLetterTransformation or {\"Letter Transformation\", transName}"];$Failed)
SetAlphabetCondition[name_String,"Alphabet",content_]:=(Print["Error: use DeclareAlphabet or ResetAlphabet to set the alphabet"];$Failed)
SetAlphabetCondition[name_String,conditionName:{__String},content_]:=(Print["Error: unknown or malformed alphabet condition ",conditionName];$Failed)
SetAlphabetCondition[name_String,conditionName_String,content_]:=(Print["Error: unknown alphabet condition \"",conditionName,"\""];$Failed)


ToList[h_][expr_]:=If[Head[expr]===h,List@@expr,{expr}]


CanonSparseArray[smat_SparseArray]:=With[{data=Select[List@@smat,#1["ExplicitLength"]>0&]},SparseArray[data,{Length[data],Splice[Dimensions[smat][[2;;All]]]}]]


SortSparseArray[smat_SparseArray]:=SparseArray[SortBy[List@@smat,#["ExplicitLength"]&]]


RowRescale[smat_SparseArray]:=With[{rescaleDen=(row|->LCM@@(Denominator/@DeleteDuplicates[row["ExplicitValues"]])*row),rescaleNum=(row|->row/GCD@@DeleteDuplicates[row["ExplicitValues"]])},SparseArray[(row|->If[MatchQ[row["ExplicitValues"],{}|{0}],row,rescaleNum[rescaleDen[row]]])/@List@@smat,Dimensions[smat]]]


Options[SmartArrayReshape]={"Verbose"->True};
SmartArrayReshape[array_?ArrayQ,dims_?VectorQ,OptionsPattern[]]:=Module[{pos,newdims,sol,x},pos=Position[dims,Verbatim[_]];
If[Length[pos]>1,Print["Error: only one unknown dimension is allowed"];Return[$Failed]];If[Length[pos]==0,Return[ArrayReshape[array,dims]]];
newdims=ReplaceAt[dims,Verbatim[_]->x,pos[[1]]];sol=Solve[Times@@Dimensions[array]==Times@@newdims,x];
If[sol==={},Print["Error: unsupported input dimensions"];Return[$Failed],If[TrueQ[OptionValue["Verbose"]],Print["Array reshape to: ",newdims/.sol[[1]]]];ArrayReshape[array,newdims/.sol[[1]]]]]


Options[CombineConditionTensor]=Replace[Options[SparseRREF`SparseRREF],Rule["Threads",_]->Rule["Threads",0],{1}];
CombineConditionTensor[tensors__SparseArray,opts:OptionsPattern[]]:=With[{dims=Dimensions/@{tensors}},
If[!Equal@@(Length/@dims),Print["Error: condition tensors must have the same rank, but got ranks ",Length/@dims];Return[$Failed]];
If[!Equal@@dims[[All,;;2]],Print["Error: condition tensors must have matching first two dimensions, but got ",dims[[All,;;2]]];Return[$Failed]];
If[!Equal@@dims[[1,;;2]],Print["Error: condition tensors must be square in the first two dimensions, but got ",dims[[1,;;2]]];Return[$Failed]];
Module[{n,mat},n=dims[[1,1]];mat=SparseRREF`SparseRREF[Flatten[RowRescale[SortSparseArray[Transpose[Join[tensors,3],{2,3,1}]]],{{1},{2,3}}],opts];If[mat===$Failed,Return[$Failed]];
SmartArrayReshape[Transpose[RowRescale[CanonSparseArray[mat]]],{n,n,_},"Verbose"->False]]]


SetAttributes[SmartMonitor,HoldFirst];
SmartMonitor[expr_,indicator_,verbose_]:=If[$Notebooks&&TrueQ[verbose],Monitor[expr,indicator],expr]


SmartPrintTemporary[expr_,verbose_]:=If[$Notebooks&&TrueQ[verbose],PrintTemporary[expr]]


sqrt/:sqrt[x_]^k_Integer/;EvenQ[k]:=x^(k/2)
sqrt/:sqrt[x_]^k_Integer/;OddQ[k]:=x^Quotient[k,2]sqrt[x]


Options[GenSqrtD]={"Verbose"->True};
GenSqrtD[alphabetExpr_?VectorQ,OptionsPattern[]]:=Module[{alpexpr,dlogexpr,vars,sqrtlist1,sqrtlist2,irrpolys,rules,print,p,i,j,time},time=AbsoluteTime[];vars=Variables[alphabetExpr];print=If[TrueQ[OptionValue["Verbose"]],Print,List];print["Independent variables: ",vars];
alpexpr=alphabetExpr/.Power[a_,k_]/;IntegerQ[k-1/2]:>Power[Factor[a],Floor[k]]sqrt[Factor[a]]//.{sqrt[a_*b_^k_Integer]:>sqrt[a*If[OddQ[k],b,1]] b^Quotient[k,2],sqrt[b_^k_Integer]:>If[OddQ[k],sqrt[b],1]* b^Quotient[k,2]};sqrtlist1=SortBy[Union@Cases[alpexpr,_sqrt,\[Infinity]],ByteCount];
(* Z2 row reduce to get independent sqrt *)
If[sqrtlist1=!={},SmartPrintTemporary[" Z2 row reducing the roots...",OptionValue["Verbose"]];irrpolys=SortBy[DeleteCases[Union@Flatten[sqrtlist1/.{sqrt[a_Times]:>List@@a}/.sqrt->List],_?NumberQ],ByteCount];
With[{irrpolys$pm=DeleteDuplicatesBy[{#,Factor[-#]}&/@irrpolys,Sort]},irrpolys=irrpolys$pm[[All,1]];sqrtlist2=sqrtlist1//.Flatten[MapIndexed[{#1[[1]]->p[#2[[1]]],#1[[2]]->-p[#2[[1]]]}&,irrpolys$pm]]];
With[{mat=NullSpace[(Exponent[#,Array[p,Length[irrpolys]]]&/@sqrtlist2[[All,1]])\[Transpose],Modulus->2]},
rules=Table[With[{terms=DeleteCases[sqrtlist2^vec,1]},terms[[-1]]->Factor[(PowerExpand[Sqrt[Times@@terms/.sqrt->Identity]]Times@@terms[[;;-2]])/(Times@@terms[[;;-2,1]])]]/.Thread[Array[p,Length[irrpolys]]->irrpolys],{vec,mat}]];
alpexpr=Factor[alpexpr/.Dispatch[rules]]/.sqrt[a_]:>sqrt[Expand[a]];sqrtlist1=SortBy[Union@Cases[alpexpr,_sqrt,\[Infinity]],ByteCount];];
If[!AllTrue[alpexpr,RationalExpressionQ[#,Join[vars,sqrtlist1]]&],Print["Error: unsupported input alphabet expression"];Return[$Failed]];
(* Rationalize the denominator *)
alpexpr=Simplify[alpexpr]/.sqrt[a_]:>sqrt[Expand[a]];dlogexpr=Factor[ExpandDenominator[Factor[D[alpexpr,{vars}]/alpexpr//.{Derivative[1][sqrt][x_]:>sqrt[x]/(2x)}]]];
With[{sepsqrt=(If[Head[#]===Times,{Times@@Select[List@@#,Not@*FreeQ[_sqrt]],Times@@Select[List@@#,FreeQ[_sqrt]]},If[FreeQ[#,_sqrt],{1,#},{#,1}]]&),
sqrtconjs=(expr|->With[{sqlist=Union@Cases[{expr},_sqrt,\[Infinity]]},Flatten[Outer[(expr/.Thread[sqlist->(sqlist*{##})])&,Sequence@@Table[{1,-1},Length[sqlist]]]][[2;;]]])},
With[{den2rat=(expr|->With[{dens=sepsqrt[Denominator[expr]]},If[NumberQ[dens[[1]]],Factor[Expand[expr]],Factor[Expand[Numerator[expr]Times@@sqrtconjs[dens[[1]]]]/(Expand[dens[[1]]Times@@sqrtconjs[dens[[1]]]]dens[[2]])]]])},
If[!FreeQ[dlogexpr,_sqrt],SmartPrintTemporary[" Rationalizing the denominator...",OptionValue["Verbose"]];
dlogexpr=SmartMonitor[Table[den2rat[dlogexpr[[i,j]]],{i,Length[alpexpr]},{j,Length[vars]}],ProgressIndicator[i*Length[vars]+j,{1,Length[alpexpr]*Length[vars]}],OptionValue["Verbose"]]]];sqrtlist1=SortBy[Union@Cases[alpexpr,_sqrt,\[Infinity]],ByteCount];];
If[!AllTrue[Denominator/@Flatten[dlogexpr],RationalExpressionQ[#,vars]&],Print["Error: unable to rationalize the denominator"];Return[$Failed]];
print["Square roots: ",sqrtlist1/.sqrt->Sqrt];print["Symbolic sqrt-reduced dlog vector generated. Time elapsed: ",AbsoluteTime[]-time];{dlogexpr,vars,sqrtlist1}]


Options[GenIntRelMat]={"Samples"->Automatic,"Tries"->100,"Threads"->0,"Verbose"->True};
GenIntRelMat[{dlogexpr_?MatrixQ,vars_?VectorQ,sqrtlist_?VectorQ},OptionsPattern[]]:=With[{n$samp=If[OptionValue["Samples"]===Automatic,10+Ceiling[Binomial[Length[dlogexpr],2]/Binomial[Length[vars],2]],OptionValue["Samples"]],
den$lcm=Times@@Union@Flatten[ToList[Times]/@Factor[Denominator/@Flatten[dlogexpr]]],genRand=(Thread[vars->RandomPrime[{2,3*Length[dlogexpr]},Length[vars]]]&)},
(* Generate numeric dlog \wedge dlog coefficient vectors *)
Module[{dlogsqexpr,mat,tocoeff,print,sq,cnt,time},time=AbsoluteTime[];print=If[TrueQ[OptionValue["Verbose"]],Print,List];print["Numeric sampling points: ",n$samp];dlogsqexpr=dlogexpr/.Dispatch[Thread[sqrtlist->Array[sq,Length[sqrtlist]]]];tocoeff=If[sqrtlist==={},Identity,Flatten[CoefficientArrays[#,Array[sq,Length[sqrtlist]]]]&];
mat=SmartMonitor[Table[With[{numrule=Dispatch[NestWhile[genRand,genRand[],(den$lcm/.#)===0&,1,OptionValue["Tries"]]]},If[(den$lcm/.numrule)===0,Print["Error: unable to avoid zero denominator within ",OptionValue["Tries"]," sampling tries"];Return[$Failed]];
With[{dlogexprnum=dlogsqexpr/.numrule,radiclistnum=sqrtlist[[All,1]]/.numrule},tocoeff/@Table[ExpandNumerator[Det[dlogexprnum[[idx$dlog,idx$dvar]]]]/.{sq[i_]^2:>radiclistnum[[i]]},{idx$dlog,Subsets[Range[Length[dlogexpr]],{2}]},{idx$dvar,Subsets[Range[Length[vars]],{2}]}]]],{cnt,n$samp}],ProgressIndicator[cnt,{1,n$samp}],OptionValue["Verbose"]];
With[{len=Max[Flatten[Map[Length,mat,{2}]]]},mat=CanonSparseArray[SparseArray[Flatten[Map[PadRight[#,len]&,mat,{2}],{{1,3},{2}}]]]];
(* Row reduce numeric dlog \wedge dlog coefficient vectors *)
SmartPrintTemporary[" Row reducing numeric system...",OptionValue["Verbose"]];mat=SparseRREF`SparseRREF[mat,"Method"->"Right","Threads"->OptionValue["Threads"]];If[mat===$Failed,Return[$Failed]];
print["Raw integrability relations generated. Time elapsed: ",AbsoluteTime[]-time];mat]]


Options[GenLettRelMat]={"Samples"->Automatic,"Tries"->100,"Threads"->0,"Verbose"->True};
GenLettRelMat[{dlogexprold_?MatrixQ,dlogexprnew_?MatrixQ,vars_?VectorQ,sqrtlist_?VectorQ},OptionsPattern[]]:=With[{n=Length[dlogexprold],n$samp=If[OptionValue["Samples"]===Automatic,10+2Ceiling[Length[dlogexprold]/Length[vars]],OptionValue["Samples"]],
den$lcm=Times@@Union@Flatten[ToList[Times]/@Factor[Denominator/@Flatten[Join[dlogexprold,dlogexprnew]]]],genRand=(Thread[vars->RandomPrime[{2,3 Max[1,Length[dlogexprold]]},Length[vars]]]&)},
(* Generate numeric old/new dlog coefficient vectors *)
Module[{dlogsqexpr,mat,tocoeff,print,sq,cnt,pivots,time},time=AbsoluteTime[];print=If[TrueQ[OptionValue["Verbose"]],Print,List];print["Numeric sampling points: ",n$samp];
dlogsqexpr=Join[dlogexprold,dlogexprnew]/.Dispatch[Thread[sqrtlist->Array[sq,Length[sqrtlist]]]];tocoeff=If[sqrtlist==={},Identity,Flatten[CoefficientArrays[#,Array[sq,Length[sqrtlist]]]]&];
mat=SmartMonitor[Table[With[{numrule=Dispatch[NestWhile[genRand,genRand[],(den$lcm/.#)===0&,1,OptionValue["Tries"]]]},If[(den$lcm/.numrule)===0,Print["Error: unable to avoid zero denominator within ",OptionValue["Tries"]," sampling tries"];Return[$Failed]];
With[{dlogexprnum=dlogsqexpr/.numrule,radiclistnum=sqrtlist[[All,1]]/.numrule},tocoeff/@Table[ExpandNumerator[dlogexprnum[[idx$dlog,idx$dvar]]]/.{sq[i_]^2:>radiclistnum[[i]]},{idx$dlog,2n},{idx$dvar,Length[vars]}]]],{cnt,n$samp}],ProgressIndicator[cnt,{1,n$samp}],OptionValue["Verbose"]];
With[{len=Max[Flatten[Map[Length,mat,{2}]]]},mat=CanonSparseArray[SparseArray[Flatten[Map[PadRight[#,len]&,mat,{2}],{{1,3},{2}}]]]];
(* Row reduce numeric old/new letter linear relation system *)
SmartPrintTemporary[" Row reducing numeric system...",OptionValue["Verbose"]];mat=SparseRREF`SparseRREF[mat,"OutputMode"->"RREF,Pivots","Method"->"Right","Threads"->OptionValue["Threads"],"Verbose"->OptionValue["Verbose"]];If[mat===$Failed,Return[$Failed]];{mat,pivots}=mat;pivots=SortBy[pivots,Last];
If[Length[pivots]<n,Print["Error: alphabet dlog vectors are not linearly independent; try increasing Samples"];Return[$Failed]];
If[Length[pivots]>n,Print["Error: mapped alphabet dlogs are not in the span of the original alphabet dlogs; try increasing Samples"];Return[$Failed]];
With[{oldblock=mat[[pivots[[All,1]],Range[n]]],newblock=mat[[pivots[[All,1]],n+1;;2n]]},
If[pivots[[All,2]]===Range[n]&&oldblock===IdentityMatrix[n,SparseArray],mat=SparseArray[Transpose[newblock]],
(* Recover the linear relation in the original old-letter basis. *)
SmartPrintTemporary[" Inverting old dlog block...",OptionValue["Verbose"]];mat=SparseRREF`SparseMatInv[oldblock,"Threads"->OptionValue["Threads"]];If[mat===$Failed,Print["Error: alphabet dlog vectors are not linearly independent; try increasing Samples"];Return[$Failed]];mat=SparseArray[Transpose[mat . newblock]]]];
print["Letter relation matrix generated. Total time elapsed: ",AbsoluteTime[]-time];mat]]


Options[GenDlogmatInt]=Options[GenIntRelMat];
GenDlogmatInt[alphabetExpr_?VectorQ,opts:OptionsPattern[]]:=Module[{n,result,upper,print,time},
time=AbsoluteTime[];print=If[TrueQ[OptionValue["Verbose"]],Print,List];n=Length[alphabetExpr];
upper=SparseArray[Flatten@Table[{i,j,(i-1)(2n-i)/2+j-i}->1,{i,n},{j,i+1,n}],{n,n,Binomial[n,2]}];
result=GenSqrtD[alphabetExpr,"Verbose"->OptionValue["Verbose"]];If[result===$Failed,Return[$Failed]];
result=GenIntRelMat[result,opts];If[result===$Failed,Return[$Failed]];
result=SparseArray[(upper-Transpose[upper,{2,1,3}]) . Transpose[CanonSparseArray[result]]];
print["Integrability condition tensor (dlogmat) generated. Total time elapsed: ",AbsoluteTime[]-time];result]


GenDlogmatCA[alphabet_?VectorQ,adjpairs_?MatrixQ]:=Module[{n,adjpairIds,S,W},n=Length[alphabet];adjpairIds=adjpairs/.Dispatch[Thread[alphabet->Range[n]]];
SmartArrayReshape[Transpose[SparseArray[NullSpace[CoefficientArrays[S@@@DeleteDuplicates[adjpairIds],Flatten[Array[S,{n,n}],1]][[2]]]]],{n,n,_},"Verbose"->False]]


GenDlogmatES[alphabet_?VectorQ,nonadjpairs_?MatrixQ]:=GenDlogmatCA[alphabet,Complement[Tuples[alphabet,{2}],Normal[nonadjpairs]]]


Options[GenLettTransMat]=Options[GenLettRelMat];
GenLettTransMat[alphabetExpr_?VectorQ,kineMap_,opts:OptionsPattern[]]:=Module[{n,result,print,time},
time=AbsoluteTime[];print=If[TrueQ[OptionValue["Verbose"]],Print,List];n=Length[alphabetExpr];
result=GenSqrtD[Join[alphabetExpr,alphabetExpr/.kineMap],"Verbose"->OptionValue["Verbose"]];If[result===$Failed,Return[$Failed]];
result=GenLettRelMat[{result[[1,;;n]],result[[1,n+1;;2n]],result[[2]],result[[3]]},opts];If[result===$Failed,Return[$Failed]];
print["Letter transformation matrix generated. Total time elapsed: ",AbsoluteTime[]-time];result]


GenFEC[alphabet_?VectorQ,firstentry_?VectorQ]:=With[{n=Length[alphabet]},With[{firstentryIds=firstentry/.Dispatch[Thread[alphabet->Range[n]]]},SparseArray[Table[{i,1,firstentryIds[[i]]}->1,{i,Length[firstentryIds]}],{Length[firstentryIds],1,n}]]]


GenLEC[alphabet_?VectorQ,lastentry_?VectorQ]:=With[{n=Length[alphabet]},With[{lastentryIds=lastentry/.Dispatch[Thread[alphabet->Range[n]]]},SparseArray[Table[{i,lastentryIds[[i]],1}->1,{i,Length[lastentryIds]}],{Length[lastentryIds],n,1}]]]


Options[GetIntegrabilityTensor]=Options[GenDlogmatInt];
GetIntegrabilityTensor[name_String,opts:OptionsPattern[]]:=Module[{result},
If[!RequireAlphabet[name],Return[$Failed]];
If[!ConditionSetQ[name,"Expression"],Print["Error: alphabet expression for \"",name,"\" has not been set"];Return[$Failed]];
If[ResultCachedQ[name,"Integrability"],Return[$result[name]["Integrability"]]];
result=GenDlogmatInt[$condition[name]["Expression"],opts];If[result===$Failed,Return[$Failed]];
SetResultValue[name,"Integrability",result];result]


Options[GetClusterAdjacencyTensor]={"Verbose"->True};
GetClusterAdjacencyTensor[name_String,OptionsPattern[]]:=Module[{result},
If[!RequireAlphabet[name],Return[$Failed]];
If[!ConditionSetQ[name,"Cluster Adjacency"],Print["Error: cluster adjacency condition for \"",name,"\" has not been set"];Return[$Failed]];
If[ResultCachedQ[name,"Cluster Adjacency"],Return[$result[name]["Cluster Adjacency"]]];
result=GenDlogmatCA[$condition[name]["Alphabet"],$condition[name]["Cluster Adjacency"]];If[result===$Failed,Return[$Failed]];
SetResultValue[name,"Cluster Adjacency",result];result]


Options[GetExtendedSteinmannTensor]={"Verbose"->True};
GetExtendedSteinmannTensor[name_String,OptionsPattern[]]:=Module[{result},
If[!RequireAlphabet[name],Return[$Failed]];
If[!ConditionSetQ[name,"Extended Steinmann"],Print["Error: extended Steinmann condition for \"",name,"\" has not been set"];Return[$Failed]];
If[ResultCachedQ[name,"Extended Steinmann"],Return[$result[name]["Extended Steinmann"]]];
result=GenDlogmatES[$condition[name]["Alphabet"],$condition[name]["Extended Steinmann"]];If[result===$Failed,Return[$Failed]];
SetResultValue[name,"Extended Steinmann",result];result]


Options[GetFirstEntryTensor]={"Verbose"->True};
GetFirstEntryTensor[name_String,OptionsPattern[]]:=Module[{result},
If[!RequireAlphabet[name],Return[$Failed]];
If[!ConditionSetQ[name,"First Entry"],Print["Error: first-entry condition for \"",name,"\" has not been set"];Return[$Failed]];
If[ResultCachedQ[name,"First Entry"],Return[$result[name]["First Entry"]]];
result=GenFEC[$condition[name]["Alphabet"],$condition[name]["First Entry"]];If[result===$Failed,Return[$Failed]];
SetResultValue[name,"First Entry",result];result]


Options[GetLastEntryTensor]={"Verbose"->True};
GetLastEntryTensor[name_String,OptionsPattern[]]:=Module[{result},
If[!RequireAlphabet[name],Return[$Failed]];
If[!ConditionSetQ[name,"Last Entry"],Print["Error: last-entry condition for \"",name,"\" has not been set"];Return[$Failed]];
If[ResultCachedQ[name,"Last Entry"],Return[$result[name]["Last Entry"]]];
result=GenLEC[$condition[name]["Alphabet"],$condition[name]["Last Entry"]];If[result===$Failed,Return[$Failed]];
SetResultValue[name,"Last Entry",result];result]


Options[GetLetterTransformationTensor]=Options[GenLettTransMat];
GetLetterTransformationTensor[name_String,transName_String,opts:OptionsPattern[]]:=Module[{result},
If[!RequireAlphabet[name],Return[$Failed]];
If[!ConditionSetQ[name,"Expression"],Print["Error: alphabet expression for \"",name,"\" has not been set"];Return[$Failed]];
If[!KeyExistsQ[$condition[name]["Letter Transformation"],transName],Print["Error: letter transformation \"",transName,"\" for alphabet \"",name,"\" has not been set"];Return[$Failed]];
If[KeyExistsQ[$result[name]["Letter Transformation"],transName],Return[$result[name]["Letter Transformation"][transName]]];
result=GenLettTransMat[$condition[name]["Expression"],$condition[name]["Letter Transformation"][transName],opts];If[result===$Failed,Return[$Failed]];
SetResultValue[name,"Letter Transformation",Join[$result[name]["Letter Transformation"],<|transName->result|>]];result]


Options[GetAlphabetConditionTensor]={"Verbose"->True};
GetVerboseOption[opts___Rule]:=("Verbose"/.{opts}/.Options[GetAlphabetConditionTensor])
$dlogmatConditionKeys={"Integrability","Cluster Adjacency","Extended Steinmann"};
GetAlphabetConditionTensor[name_String,"Integrability",opts___Rule]:=GetIntegrabilityTensor[name,"Verbose"->GetVerboseOption[opts]]
GetAlphabetConditionTensor[name_String,"Cluster Adjacency",opts___Rule]:=GetClusterAdjacencyTensor[name,"Verbose"->GetVerboseOption[opts]]
GetAlphabetConditionTensor[name_String,"Extended Steinmann",opts___Rule]:=GetExtendedSteinmannTensor[name,"Verbose"->GetVerboseOption[opts]]
GetAlphabetConditionTensor[name_String,"First Entry",opts___Rule]:=GetFirstEntryTensor[name,"Verbose"->GetVerboseOption[opts]]
GetAlphabetConditionTensor[name_String,"Last Entry",opts___Rule]:=GetLastEntryTensor[name,"Verbose"->GetVerboseOption[opts]]
GetAlphabetConditionTensor[name_String,{"Letter Transformation",transName_String},opts___Rule]:=GetLetterTransformationTensor[name,transName,"Verbose"->GetVerboseOption[opts]]
GetAlphabetConditionTensor[name_String,{},opts___Rule]:=(Print["Error: no alphabet condition tensors were requested"];$Failed)
GetAlphabetConditionTensor[name_String,conditionNames:{__String},opts___Rule]:=Module[{invalid,tensors,verbose},
invalid=Complement[conditionNames,$dlogmatConditionKeys];
If[invalid=!={},Print["Error: cannot combine non-dlogmat or unknown condition tensors ",invalid];Return[$Failed]];
verbose=GetVerboseOption[opts];
tensors=GetAlphabetConditionTensor[name,#,"Verbose"->verbose]&/@conditionNames;If[MemberQ[tensors,$Failed],Return[$Failed]];
CombineConditionTensor[Sequence@@tensors,"Verbose"->verbose]]
GetAlphabetConditionTensor[name_String,conditionName_String,opts___Rule]:=(Print["Error: unknown alphabet condition tensor \"",conditionName,"\""];$Failed)


End[]

EndPackage[]
