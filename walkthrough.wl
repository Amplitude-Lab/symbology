(* ::Package:: *)

(* ::Section:: *)
(*Symbol Bootstrap*)


(* ::Subsection:: *)
(*Load Package*)


SetDirectory[NotebookDirectory[]];
Get["SymbolBootstrap.wl"]


(* ::Subsection:: *)
(*Example: heptagon E6 letters a-like*)


(* ::Subsubsection:: *)
(*Momentum twistor helper*)


abrules[mat_?MatrixQ]:=With[{detrules=(ab[as___]:>Det[mat[[{as}]]])},{ab[a_,cap1[{b_,c_},{d_,e_},{f_,g_}]]:>(ab[a,b,d,e]ab[a,c,f,g]-ab[a,b,f,g]ab[a,c,d,e]/.detrules),detrules}]


(* ::Subsubsection:: *)
(*Declare alphabet with letters*)


alphabetE6={a11,a12,a13,a14,a15,a16,a17,a21,a22,a23,a24,a25,a26,a27,a31,a32,a33,a34,a35,a36,a37,a41,a42,a43,a44,a45,a46,a47,a51,a52,a53,a54,a55,a56,a57,a61,a62,a63,a64,a65,a66,a67};


DeclareAlphabet["E6",alphabetE6]


(* ::Subsubsection:: *)
(*Set alphabet parameterization*)


MatE6={{1,0,0,0},{(1+f5)/f5,1,0,0},{(1+f2+f1 f2+f1 f2 f4+f1 f2 f3 f4)/(f1 f2 f3 f4 f5),(1+f2+f1 f2+f6+f2 f6+f1 f2 f6+f1 f2 f4 f6+f1 f2 f3 f4 f6)/(f1 f2 f3 f4 f6),1,0},{(1+f2+f1 f2+f1 f2 f4)/(f1 f2 f3 f4 f5),((1+f2+f1 f2+f1 f2 f4) (1+f6))/(f1 f2 f3 f4 f6),1+f2+f1 f2+f1 f2 f4,1},{0,1/(f3 f6),f2 (1+f1+f1 f4),1},{0,0,(f2 (1+f2+f1 f2+f1 f2 f4))/(1+f2),1},{0,0,0,1}};


rewrite={a11->(ab[1,2,3,4] ab[1,5,6,7] ab[2,3,6,7])/(ab[1,2,3,7] ab[1,2,6,7] ab[3,4,5,6]),a21->(ab[1,2,3,4] ab[2,5,6,7])/(ab[1,2,6,7] ab[2,3,4,5]),a31->(ab[1,5,6,7] ab[2,3,4,7])/(ab[1,2,3,7] ab[4,5,6,7]),a41->(ab[2,4,5,7] ab[3,4,5,6])/(ab[2,3,4,5] ab[4,5,6,7]),a51->ab[1,cap1[{2,3},{4,5},{6,7}]]/(ab[1,2,3,4] ab[1,5,6,7]),a61->-(ab[1,cap1[{2,7},{3,4},{5,6}]]/(ab[1,2,3,4] ab[1,5,6,7])),a12->(ab[1,2,6,7] ab[1,3,4,7] ab[2,3,4,5])/(ab[1,2,3,4] ab[1,2,3,7] ab[4,5,6,7]),a22->(ab[1,3,6,7] ab[2,3,4,5])/(ab[1,2,3,7] ab[3,4,5,6]),a32->(ab[1,2,6,7] ab[1,3,4,5])/(ab[1,2,3,4] ab[1,5,6,7]),a42->(ab[1,3,5,6] ab[4,5,6,7])/(ab[1,5,6,7] ab[3,4,5,6]),a52->ab[2,cap1[{1,7},{3,4},{5,6}]]/(ab[1,2,6,7] ab[2,3,4,5]),a62->-(ab[2,cap1[{1,3},{4,5},{6,7}]]/(ab[1,2,6,7] ab[2,3,4,5])),a13->(ab[1,2,3,7] ab[1,2,4,5] ab[3,4,5,6])/(ab[1,2,3,4] ab[1,5,6,7] ab[2,3,4,5]),a23->(ab[1,2,4,7] ab[3,4,5,6])/(ab[1,2,3,4] ab[4,5,6,7]),a33->(ab[1,2,3,7] ab[2,4,5,6])/(ab[1,2,6,7] ab[2,3,4,5]),a43->(ab[1,5,6,7] ab[2,4,6,7])/(ab[1,2,6,7] ab[4,5,6,7]),a53->ab[3,cap1[{1,2},{4,5},{6,7}]]/(ab[1,2,3,7] ab[3,4,5,6]),a63->-(ab[3,cap1[{1,7},{2,4},{5,6}]]/(ab[1,2,3,7] ab[3,4,5,6])),a14->(ab[1,2,3,4] ab[2,3,5,6] ab[4,5,6,7])/(ab[1,2,6,7] ab[2,3,4,5] ab[3,4,5,6]),a24->(ab[1,2,3,5] ab[4,5,6,7])/(ab[1,5,6,7] ab[2,3,4,5]),a34->(ab[1,2,3,4] ab[3,5,6,7])/(ab[1,2,3,7] ab[3,4,5,6]),a44->(ab[1,2,6,7] ab[1,3,5,7])/(ab[1,2,3,7] ab[1,5,6,7]),a54->ab[4,cap1[{1,7},{2,3},{5,6}]]/(ab[1,2,3,4] ab[4,5,6,7]),a64->-(ab[4,cap1[{1,2},{3,5},{6,7}]]/(ab[1,2,3,4] ab[4,5,6,7])),a15->(ab[1,5,6,7] ab[2,3,4,5] ab[3,4,6,7])/(ab[1,2,3,7] ab[3,4,5,6] ab[4,5,6,7]),a25->(ab[1,5,6,7] ab[2,3,4,6])/(ab[1,2,6,7] ab[3,4,5,6]),a35->(ab[1,4,6,7] ab[2,3,4,5])/(ab[1,2,3,4] ab[4,5,6,7]),a45->(ab[1,2,3,7] ab[1,2,4,6])/(ab[1,2,3,4] ab[1,2,6,7]),a55->ab[5,cap1[{1,2},{3,4},{6,7}]]/(ab[1,5,6,7] ab[2,3,4,5]),a65->-(ab[5,cap1[{1,7},{2,3},{4,6}]]/(ab[1,5,6,7] ab[2,3,4,5])),a16->(ab[1,2,6,7] ab[1,4,5,7] ab[3,4,5,6])/(ab[1,2,3,4] ab[1,5,6,7] ab[4,5,6,7]),a26->(ab[1,2,6,7] ab[3,4,5,7])/(ab[1,2,3,7] ab[4,5,6,7]),a36->(ab[1,2,5,7] ab[3,4,5,6])/(ab[1,5,6,7] ab[2,3,4,5]),a46->(ab[1,2,3,4] ab[2,3,5,7])/(ab[1,2,3,7] ab[2,3,4,5]),a56->ab[6,cap1[{1,7},{2,3},{4,5}]]/(ab[1,2,6,7] ab[3,4,5,6]),a66->-(ab[6,cap1[{1,2},{3,4},{5,7}]]/(ab[1,2,6,7] ab[3,4,5,6])),a17->(ab[1,2,3,7] ab[1,2,5,6] ab[4,5,6,7])/(ab[1,2,6,7] ab[1,5,6,7] ab[2,3,4,5]),a27->(ab[1,2,3,7] ab[1,4,5,6])/(ab[1,2,3,4] ab[1,5,6,7]),a37->(ab[1,2,3,6] ab[4,5,6,7])/(ab[1,2,6,7] ab[3,4,5,6]),a47->(ab[1,3,4,6] ab[2,3,4,5])/(ab[1,2,3,4] ab[3,4,5,6]),a57->ab[7,cap1[{1,2},{3,4},{5,6}]]/(ab[1,2,3,7] ab[4,5,6,7]),a67->-(ab[7,cap1[{1,6},{2,3},{4,5}]]/(ab[1,2,3,7] ab[4,5,6,7]))};


(*alphabets in momentum-twistor parameterization*)
alphabetExprE6=Simplify[alphabetE6/.rewrite/.abrules[MatE6]];


SetAlphabetExpression["E6",alphabetExprE6]


(* ::Subsubsection:: *)
(*Get integrability condition tensor*)


dlogmatE6Int=GetIntegrabilityTensor["E6"]


(* Equivalently *)
dlogmatE6Int==GetAlphabetConditionTensor["E6","Integrability"]


(* ::Subsubsection:: *)
(*Set extended Steinmann*)


nonadjpairE6ES={{a11,a12},{a11,a13},{a12,a13},{a12,a14},{a13,a14},{a13,a15},{a14,a15},{a14,a16},
{a15,a16},{a15,a17},{a16,a17},{a16,a11},{a17,a11},{a17,a12}};


SetExtendedSteinmann["E6",nonadjpairE6ES]


(* ::Subsubsection:: *)
(*Get extended Steinmann condition tensor*)


dlogmatE6ES=GetExtendedSteinmannTensor["E6"]


(* Equivalently *)
dlogmatE6ES==GetAlphabetConditionTensor["E6","Extended Steinmann"]


(* ::Subsubsection:: *)
(*Set cluster adjacency*)


adjpairE6CA={{a11,a11},{a11,a14},{a11,a15},{a11,a21},{a11,a22},{a11,a24},{a11,a25},{a11,a26},{a11,a31},{a11,a33},{a11,a34},{a11,a35},{a11,a37},{a11,a41},{a11,a43},{a11,a46},{a11,a51},{a11,a53},{a11,a56},{a11,a62},{a11,a67},{a12,a12},{a12,a15},{a12,a16},{a12,a22},{a12,a23},{a12,a25},{a12,a26},{a12,a27},{a12,a31},{a12,a32},{a12,a34},{a12,a35},{a12,a36},{a12,a42},{a12,a44},{a12,a47},{a12,a52},{a12,a54},{a12,a57},{a12,a61},{a12,a63},{a13,a13},{a13,a16},{a13,a17},{a13,a21},{a13,a23},{a13,a24},{a13,a26},{a13,a27},{a13,a32},{a13,a33},{a13,a35},{a13,a36},{a13,a37},{a13,a41},{a13,a43},{a13,a45},{a13,a51},{a13,a53},{a13,a55},{a13,a62},{a13,a64},{a14,a11},{a14,a14},{a14,a17},{a14,a21},{a14,a22},{a14,a24},{a14,a25},{a14,a27},{a14,a31},{a14,a33},{a14,a34},{a14,a36},{a14,a37},{a14,a42},{a14,a44},{a14,a46},{a14,a52},{a14,a54},{a14,a56},{a14,a63},{a14,a65},{a15,a11},{a15,a12},{a15,a15},{a15,a21},{a15,a22},{a15,a23},{a15,a25},{a15,a26},{a15,a31},{a15,a32},{a15,a34},{a15,a35},{a15,a37},{a15,a43},{a15,a45},{a15,a47},{a15,a53},{a15,a55},{a15,a57},{a15,a64},{a15,a66},{a16,a12},{a16,a13},{a16,a16},{a16,a22},{a16,a23},{a16,a24},{a16,a26},{a16,a27},{a16,a31},{a16,a32},{a16,a33},{a16,a35},{a16,a36},{a16,a41},{a16,a44},{a16,a46},{a16,a51},{a16,a54},{a16,a56},{a16,a65},{a16,a67},{a17,a13},{a17,a14},{a17,a17},{a17,a21},{a17,a23},{a17,a24},{a17,a25},{a17,a27},{a17,a32},{a17,a33},{a17,a34},{a17,a36},{a17,a37},{a17,a42},{a17,a45},{a17,a47},{a17,a52},{a17,a55},{a17,a57},{a17,a61},{a17,a66},{a21,a11},{a21,a13},{a21,a14},{a21,a15},{a21,a17},{a21,a21},{a21,a23},{a21,a24},{a21,a25},{a21,a26},{a21,a31},{a21,a33},{a21,a34},{a21,a36},{a21,a37},{a21,a41},{a21,a43},{a21,a45},{a21,a46},{a21,a52},{a21,a53},{a21,a55},{a21,a57},{a21,a62},{a21,a64},{a21,a66},{a22,a11},{a22,a12},{a22,a14},{a22,a15},{a22,a16},{a22,a22},{a22,a24},{a22,a25},{a22,a26},{a22,a27},{a22,a31},{a22,a32},{a22,a34},{a22,a35},{a22,a37},{a22,a42},{a22,a44},{a22,a46},{a22,a47},{a22,a51},{a22,a53},{a22,a54},{a22,a56},{a22,a63},{a22,a65},{a22,a67},{a23,a12},{a23,a13},{a23,a15},{a23,a16},{a23,a17},{a23,a21},{a23,a23},{a23,a25},{a23,a26},{a23,a27},{a23,a31},{a23,a32},{a23,a33},{a23,a35},{a23,a36},{a23,a41},{a23,a43},{a23,a45},{a23,a47},{a23,a52},{a23,a54},{a23,a55},{a23,a57},{a23,a61},{a23,a64},{a23,a66},{a24,a11},{a24,a13},{a24,a14},{a24,a16},{a24,a17},{a24,a21},{a24,a22},{a24,a24},{a24,a26},{a24,a27},{a24,a32},{a24,a33},{a24,a34},{a24,a36},{a24,a37},{a24,a41},{a24,a42},{a24,a44},{a24,a46},{a24,a51},{a24,a53},{a24,a55},{a24,a56},{a24,a62},{a24,a65},{a24,a67},{a25,a11},{a25,a12},{a25,a14},{a25,a15},{a25,a17},{a25,a21},{a25,a22},{a25,a23},{a25,a25},{a25,a27},{a25,a31},{a25,a33},{a25,a34},{a25,a35},{a25,a37},{a25,a42},{a25,a43},{a25,a45},{a25,a47},{a25,a52},{a25,a54},{a25,a56},{a25,a57},{a25,a61},{a25,a63},{a25,a66},{a26,a11},{a26,a12},{a26,a13},{a26,a15},{a26,a16},{a26,a21},{a26,a22},{a26,a23},{a26,a24},{a26,a26},{a26,a31},{a26,a32},{a26,a34},{a26,a35},{a26,a36},{a26,a41},{a26,a43},{a26,a44},{a26,a46},{a26,a51},{a26,a53},{a26,a55},{a26,a57},{a26,a62},{a26,a64},{a26,a67},{a27,a12},{a27,a13},{a27,a14},{a27,a16},{a27,a17},{a27,a22},{a27,a23},{a27,a24},{a27,a25},{a27,a27},{a27,a32},{a27,a33},{a27,a35},{a27,a36},{a27,a37},{a27,a42},{a27,a44},{a27,a45},{a27,a47},{a27,a51},{a27,a52},{a27,a54},{a27,a56},{a27,a61},{a27,a63},{a27,a65},{a31,a11},{a31,a12},{a31,a14},{a31,a15},{a31,a16},{a31,a21},{a31,a22},{a31,a23},{a31,a25},{a31,a26},{a31,a31},{a31,a33},{a31,a34},{a31,a35},{a31,a36},{a31,a41},{a31,a43},{a31,a44},{a31,a46},{a31,a52},{a31,a54},{a31,a56},{a31,a57},{a31,a63},{a31,a65},{a31,a67},{a32,a12},{a32,a13},{a32,a15},{a32,a16},{a32,a17},{a32,a22},{a32,a23},{a32,a24},{a32,a26},{a32,a27},{a32,a32},{a32,a34},{a32,a35},{a32,a36},{a32,a37},{a32,a42},{a32,a44},{a32,a45},{a32,a47},{a32,a51},{a32,a53},{a32,a55},{a32,a57},{a32,a61},{a32,a64},{a32,a66},{a33,a11},{a33,a13},{a33,a14},{a33,a16},{a33,a17},{a33,a21},{a33,a23},{a33,a24},{a33,a25},{a33,a27},{a33,a31},{a33,a33},{a33,a35},{a33,a36},{a33,a37},{a33,a41},{a33,a43},{a33,a45},{a33,a46},{a33,a51},{a33,a52},{a33,a54},{a33,a56},{a33,a62},{a33,a65},{a33,a67},{a34,a11},{a34,a12},{a34,a14},{a34,a15},{a34,a17},{a34,a21},{a34,a22},{a34,a24},{a34,a25},{a34,a26},{a34,a31},{a34,a32},{a34,a34},{a34,a36},{a34,a37},{a34,a42},{a34,a44},{a34,a46},{a34,a47},{a34,a52},{a34,a53},{a34,a55},{a34,a57},{a34,a61},{a34,a63},{a34,a66},{a35,a11},{a35,a12},{a35,a13},{a35,a15},{a35,a16},{a35,a22},{a35,a23},{a35,a25},{a35,a26},{a35,a27},{a35,a31},{a35,a32},{a35,a33},{a35,a35},{a35,a37},{a35,a41},{a35,a43},{a35,a45},{a35,a47},{a35,a51},{a35,a53},{a35,a54},{a35,a56},{a35,a62},{a35,a64},{a35,a67},{a36,a12},{a36,a13},{a36,a14},{a36,a16},{a36,a17},{a36,a21},{a36,a23},{a36,a24},{a36,a26},{a36,a27},{a36,a31},{a36,a32},{a36,a33},{a36,a34},{a36,a36},{a36,a41},{a36,a42},{a36,a44},{a36,a46},{a36,a52},{a36,a54},{a36,a55},{a36,a57},{a36,a61},{a36,a63},{a36,a65},{a37,a11},{a37,a13},{a37,a14},{a37,a15},{a37,a17},{a37,a21},{a37,a22},{a37,a24},{a37,a25},{a37,a27},{a37,a32},{a37,a33},{a37,a34},{a37,a35},{a37,a37},{a37,a42},{a37,a43},{a37,a45},{a37,a47},{a37,a51},{a37,a53},{a37,a55},{a37,a56},{a37,a62},{a37,a64},{a37,a66},{a41,a11},{a41,a13},{a41,a16},{a41,a21},{a41,a23},{a41,a24},{a41,a26},{a41,a31},{a41,a33},{a41,a35},{a41,a36},{a41,a41},{a41,a43},{a41,a46},{a41,a51},{a41,a62},{a41,a67},{a42,a12},{a42,a14},{a42,a17},{a42,a22},{a42,a24},{a42,a25},{a42,a27},{a42,a32},{a42,a34},{a42,a36},{a42,a37},{a42,a42},{a42,a44},{a42,a47},{a42,a52},{a42,a61},{a42,a63},{a43,a11},{a43,a13},{a43,a15},{a43,a21},{a43,a23},{a43,a25},{a43,a26},{a43,a31},{a43,a33},{a43,a35},{a43,a37},{a43,a41},{a43,a43},{a43,a45},{a43,a53},{a43,a62},{a43,a64},{a44,a12},{a44,a14},{a44,a16},{a44,a22},{a44,a24},{a44,a26},{a44,a27},{a44,a31},{a44,a32},{a44,a34},{a44,a36},{a44,a42},{a44,a44},{a44,a46},{a44,a54},{a44,a63},{a44,a65},{a45,a13},{a45,a15},{a45,a17},{a45,a21},{a45,a23},{a45,a25},{a45,a27},{a45,a32},{a45,a33},{a45,a35},{a45,a37},{a45,a43},{a45,a45},{a45,a47},{a45,a55},{a45,a64},{a45,a66},{a46,a11},{a46,a14},{a46,a16},{a46,a21},{a46,a22},{a46,a24},{a46,a26},{a46,a31},{a46,a33},{a46,a34},{a46,a36},{a46,a41},{a46,a44},{a46,a46},{a46,a56},{a46,a65},{a46,a67},{a47,a12},{a47,a15},{a47,a17},{a47,a22},{a47,a23},{a47,a25},{a47,a27},{a47,a32},{a47,a34},{a47,a35},{a47,a37},{a47,a42},{a47,a45},{a47,a47},{a47,a57},{a47,a61},{a47,a66},{a51,a11},{a51,a13},{a51,a16},{a51,a22},{a51,a24},{a51,a26},{a51,a27},{a51,a32},{a51,a33},{a51,a35},{a51,a37},{a51,a41},{a51,a51},{a51,a53},{a51,a56},{a51,a62},{a51,a67},{a52,a12},{a52,a14},{a52,a17},{a52,a21},{a52,a23},{a52,a25},{a52,a27},{a52,a31},{a52,a33},{a52,a34},{a52,a36},{a52,a42},{a52,a52},{a52,a54},{a52,a57},{a52,a61},{a52,a63},{a53,a11},{a53,a13},{a53,a15},{a53,a21},{a53,a22},{a53,a24},{a53,a26},{a53,a32},{a53,a34},{a53,a35},{a53,a37},{a53,a43},{a53,a51},{a53,a53},{a53,a55},{a53,a62},{a53,a64},{a54,a12},{a54,a14},{a54,a16},{a54,a22},{a54,a23},{a54,a25},{a54,a27},{a54,a31},{a54,a33},{a54,a35},{a54,a36},{a54,a44},{a54,a52},{a54,a54},{a54,a56},{a54,a63},{a54,a65},{a55,a13},{a55,a15},{a55,a17},{a55,a21},{a55,a23},{a55,a24},{a55,a26},{a55,a32},{a55,a34},{a55,a36},{a55,a37},{a55,a45},{a55,a53},{a55,a55},{a55,a57},{a55,a64},{a55,a66},{a56,a11},{a56,a14},{a56,a16},{a56,a22},{a56,a24},{a56,a25},{a56,a27},{a56,a31},{a56,a33},{a56,a35},{a56,a37},{a56,a46},{a56,a51},{a56,a54},{a56,a56},{a56,a65},{a56,a67},{a57,a12},{a57,a15},{a57,a17},{a57,a21},{a57,a23},{a57,a25},{a57,a26},{a57,a31},{a57,a32},{a57,a34},{a57,a36},{a57,a47},{a57,a52},{a57,a55},{a57,a57},{a57,a61},{a57,a66},{a61,a12},{a61,a17},{a61,a23},{a61,a25},{a61,a27},{a61,a32},{a61,a34},{a61,a36},{a61,a42},{a61,a47},{a61,a52},{a61,a57},{a61,a61},{a62,a11},{a62,a13},{a62,a21},{a62,a24},{a62,a26},{a62,a33},{a62,a35},{a62,a37},{a62,a41},{a62,a43},{a62,a51},{a62,a53},{a62,a62},{a63,a12},{a63,a14},{a63,a22},{a63,a25},{a63,a27},{a63,a31},{a63,a34},{a63,a36},{a63,a42},{a63,a44},{a63,a52},{a63,a54},{a63,a63},{a64,a13},{a64,a15},{a64,a21},{a64,a23},{a64,a26},{a64,a32},{a64,a35},{a64,a37},{a64,a43},{a64,a45},{a64,a53},{a64,a55},{a64,a64},{a65,a14},{a65,a16},{a65,a22},{a65,a24},{a65,a27},{a65,a31},{a65,a33},{a65,a36},{a65,a44},{a65,a46},{a65,a54},{a65,a56},{a65,a65},{a66,a15},{a66,a17},{a66,a21},{a66,a23},{a66,a25},{a66,a32},{a66,a34},{a66,a37},{a66,a45},{a66,a47},{a66,a55},{a66,a57},{a66,a66},{a67,a11},{a67,a16},{a67,a22},{a67,a24},{a67,a26},{a67,a31},{a67,a33},{a67,a35},{a67,a41},{a67,a46},{a67,a51},{a67,a56},{a67,a67}};


SetClusterAdjacency["E6",adjpairE6CA]


(* ::Subsubsection:: *)
(*Get cluster adjacency condition tensor*)


dlogmatE6CA=GetClusterAdjacencyTensor["E6"]


(* Equivalently *)
dlogmatE6CA==GetAlphabetConditionTensor["E6","Cluster Adjacency"]


(* ::Subsubsection:: *)
(*Get condition tensor of all above*)


dlogmatE6=GetAlphabetConditionTensor["E6",{"Integrability","Extended Steinmann","Cluster Adjacency"}]


(* ::Subsubsection:: *)
(*Set first entry*)


SetFirstEntry["E6",{a11,a12,a13,a14,a15,a16,a17}]


(* ::Subsubsection:: *)
(*Get FEC_1*)


FEC[1]=GetFirstEntryTensor["E6"]


(* ::Subsubsection:: *)
(*Set last entry*)


SetLastEntry["E6",{a21,a22,a23,a24,a25,a26,a27,a31,a32,a33,a34,a35,a36,a37}]


(* ::Subsubsection:: *)
(*Get LEC_1*)


LEC[1]=GetLastEntryTensor["E6"]


(* ::Subsection:: *)
(*Example: pentagon integrability and letter transformations*)


(* ::Subsubsection:: *)
(*Declare pentagon alphabet with square-root letters*)


alphabetPentagon=Table[W[i],{i,31}];
pentagonDelta=s12^2 s15^2-2s12^2 s15 s23+s12^2 s23^2+2s12 s15 s23 s34-2s12 s23^2 s34+s23^2 s34^2-2s12 s15^2 s45+2s12 s15 s23 s45+2s12 s15 s34 s45+2s12 s23 s34 s45+2s15 s23 s34 s45-2s23 s34^2 s45+s15^2 s45^2-2s15 s34 s45^2+s34^2 s45^2;
alphabetExprPentagon={s12,s23,s34,s45,s15,s34+s45,s15+s45,s12+s15,s12+s23,s23+s34,s12-s45,-s15+s23,-s12+s34,-s23+s45,s15-s34,s12+s23-s45,-s15+s23+s34,-s12+s34+s45,s15-s23+s45,s12+s15-s34,-s12-s23+s34+s45,s15-s23-s34+s45,s12+s15-s34-s45,s12-s15+s23-s45,-s12-s15+s23+s34,
(-Sqrt[pentagonDelta]-s12 s15+s12 s23-s23 s34-s15 s45+s34 s45)/(Sqrt[pentagonDelta]-s12 s15+s12 s23-s23 s34-s15 s45+s34 s45),
(-Sqrt[pentagonDelta]-s12 s15-s12 s23+s23 s34+s15 s45-s34 s45)/(Sqrt[pentagonDelta]-s12 s15-s12 s23+s23 s34+s15 s45-s34 s45),
(-Sqrt[pentagonDelta]+s12 s15-s12 s23-s23 s34-s15 s45+s34 s45)/(Sqrt[pentagonDelta]+s12 s15-s12 s23-s23 s34-s15 s45+s34 s45),
(-Sqrt[pentagonDelta]-s12 s15+s12 s23-s23 s34+s15 s45-s34 s45)/(Sqrt[pentagonDelta]-s12 s15+s12 s23-s23 s34+s15 s45-s34 s45),
(-Sqrt[pentagonDelta]+s12 s15-s12 s23+s23 s34-s15 s45-s34 s45)/(Sqrt[pentagonDelta]+s12 s15-s12 s23+s23 s34-s15 s45-s34 s45),pentagonDelta};


DeclareAlphabet["Pentagon",alphabetPentagon]


SetAlphabetExpression["Pentagon",alphabetExprPentagon]


(* ::Subsubsection:: *)
(*Get integrability condition tensor*)


dlogmatPentagonInt=GetIntegrabilityTensor["Pentagon"]


(* Equivalently *)
dlogmatPentagonInt==GetAlphabetConditionTensor["Pentagon","Integrability"]


(* ::Subsubsection:: *)
(*Cyclic transformation*)


cyclicMapPentagon={s12->s23,s23->s34,s34->s45,s45->s15,s15->s12};


SetLetterTransformation["Pentagon","Cyclic",cyclicMapPentagon]


cycmatPentagon=GetLetterTransformationTensor["Pentagon","Cyclic"]


MatrixPower[cycmatPentagon,5]==IdentityMatrix[31,SparseArray]


(* Equivalently *)
cycmatPentagon==GetAlphabetConditionTensor["Pentagon",{"Letter Transformation","Cyclic"}]


(* ::Subsubsection:: *)
(*Flip transformation*)


flipMapPentagon={s12->s15,s23->s45,s34->s34,s45->s23,s15->s12};


SetLetterTransformation["Pentagon","Flip",flipMapPentagon]


flipmatPentagon=GetLetterTransformationTensor["Pentagon","Flip"]


MatrixPower[flipmatPentagon,2]==IdentityMatrix[31,SparseArray]


(* Equivalently *)
flipmatPentagon==GetAlphabetConditionTensor["Pentagon",{"Letter Transformation","Flip"}]


(* ::Subsection:: *)
(*Example: 4pt Form Factor integrability and letter transformations*)


(* ::Subsubsection:: *)
(*Declare 4pt form factor alphabet with square-root letters*)


alphabetFF4pt=Table[W[i],{i,93}];
alphabetExprFF4pt=BinaryDeserialize[ByteArray[CompressedData["
1:eJwBxQs69CFib1JiAQAAALgLAAA4Qzp4nO2dzY7cNhaF4arq+CfeZJONgbxD
VT9CL7KZRTDIdoBxgDRgwDMG0rHz8ok9lXHSJVHn3stL/VSp/O0aXZREkZff
OZeUqPt/Pez+8ebh14fn379999Prt/9+vz/9eTj9eXu/e9j98Pb9w92Tzu/3
m4ebH9/85+eHu0+P//2wl/89dAp0/ns6rTqqW5lbfVpVL1kD4wR/lnVqEP5+
uN9+/j28vDgmvGdRj+A6qqDup8da6Iufmjaq5SHuJHEl/97s0CoKBmGhK5yK
j2H/xtUcHFPRRje1/bt1u6Mz/szGvvnh3W8//9L5391GdYG4UH/wDqps3duw
Rvt+eATjZVjlw7HKz0/x2f/1/fGG5I18tNFSnOH2eP4JGnrUXYftX4eeuKKD
waJ67JCq2qkeRud99TcSTkU3+kZ0/7qtG7ZjvyeeeqCT9fNDuj/c1eUTBSxB
DVrwNmggv4FvO6cykDMqgvz+/T+Q8tHfSJ+nzijSGDF0MJCVhBroFu+NwBtP
vPudE3JsYEdkpLdVcMKIrO7e3uh/USVsuiODLquWhUPUclI3zKExu7C4XqSJ
x9XClAjMRm5Wj5yKlhTB5QtDmxpkXEeodoG0B8NpvPIPg6t+bDbdfhDb9UN7
iqE1SKp+L4LHk8KN3eaOxepkrl0DHIzjbtlOVOYs9EBbQzlM0DGh8YOfy07x
g65b9uVkdkF0a6s1yBHgscAfBgO2Tmv2S6qWO8h6hZ0gAsMI7G+KMRRUsxtV
G1n21NRRMtnplI9R2hgMtr8czeC/2xq7nZk3Mu+scPnW3E1tNG2MfnejyaBn
VNumIX1rlP3ajKdAGw6JTpEeMBo00iFVu8JQT8o8xTmrYUa3bu+2d4kXCR/2
/Uiwu6/a1M/GvIyyDGsQt44fFofqsNAlVVAcrNw8OKvhB7d/01wKqmLsXnKz
G4om5bfBeRM+pWpcmlnrPqhrT8WDumYsQSwX1elDVKso9lPeoMGgRa5Mx0jj
zJ02Um2akOjEetGwpKBR84PucqbUAjin5q3HVjMjSK+8xYZOVMWrEf4AEv/9
k6KPVdX5V7FUW2dADxp4Rg0SZleOgKjicVltgT/qVuiWvakf8+Elqmx0XQ4f
eq7apKF/gUjbuoEf6mAxpdvQjp2r7WJ4V18idZg1qNVhRjRIPIrmK82pCiKj
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
P4d/p4ceSuAu
"]]]


DeclareAlphabet["FF4pt",alphabetFF4pt]


SetAlphabetExpression["FF4pt",alphabetExprFF4pt]


(* ::Subsubsection:: *)
(*Get integrability condition tensor*)


dlogmatFF4ptInt=GetIntegrabilityTensor["FF4pt"]


(* Equivalently *)
dlogmatFF4ptInt==GetAlphabetConditionTensor["FF4pt","Integrability"]


(* ::Subsubsection:: *)
(*Cyclic transformation*)


cyclicMapFF4pt={u1->u2,u2->u3,u3->1+u2-v1-v2,v1->v2,v2->1+u3-u1-v2};


SetLetterTransformation["FF4pt","Cyclic",cyclicMapFF4pt]


cycmatFF4pt=GetLetterTransformationTensor["FF4pt","Cyclic"]


MatrixPower[cycmatFF4pt,4]==IdentityMatrix[93,SparseArray]


(* Equivalently *)
cycmatFF4pt==GetAlphabetConditionTensor["FF4pt",{"Letter Transformation","Cyclic"}]


(* ::Subsubsection:: *)
(*Flip transformation*)


flipMapFF4pt={u1->1+u2-v1-v2,u2->u3,u3->u2,v1->1+u3-u1-v2,v2->v2};


SetLetterTransformation["FF4pt","Flip",flipMapFF4pt]


flipmatFF4pt=GetLetterTransformationTensor["FF4pt","Flip"]


MatrixPower[flipmatFF4pt,2]==IdentityMatrix[93,SparseArray]


(* Equivalently *)
flipmatFF4pt==GetAlphabetConditionTensor["FF4pt",{"Letter Transformation","Flip"}]
