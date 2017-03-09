(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Mar 8, 2017 *)

BeginPackage["CoupledModes`"]
(* Exported symbols added here with SymbolName::usage *) 

SymmeterizeModes::usage =
	"Adds in coupling terms corresponding to the negative frequencies"

CoupledModeSystem::usage =
	"Gives a langevin matrix"

ShowGraph::usage = "Displays the full coupling graph for the langevin matrix"

ShowGraphSimplified::usage = "Displays the full coupling graph for the langevin matrix, but with edges collapsed"

ConnectedModes::usage = "List of modes which are connected in the langevin graph"

ExtractModes::usage = "Reduce a set of modes and conversion/gain processes to ones given in a list"

Cofactor::usage = "Cofactor[a,i,j] of a matrix a C_ij"

Begin["`Private`"]
(* Implementation of the package *)

SymmeterizeModes[np_Integer, nq_Integer, modes__, conversion__, gain__] := Module[
  {m, c, g},
  m = Join[modes, modes];
  c = Join[{#[[1]], #[[2]] + np, #[[3]] + np} & /@ conversion, conversion];
  g = Join[{#[[1]], #[[3]], #[[2]] + np} & /@ gain, {#[[1]], #[[2]], #[[3]] + np} & /@ gain];
  {np, nq, m, c, g}
  ]

CoupledModeSystem[np_Integer,nq_Integer,modes__,conversion__,gain__] := Module[
	{elt, bsector,cprocs,asector,dsector,diag,sym},

	bsector = Sort[{#[[2]], #[[3]]}] -> #[[1]] & /@ gain;
	cprocs = GroupBy[conversion, ((# > np &) /@ # &)@*Rest];
	elt = Lookup[cprocs, {{False, False}, {True, True}}, {}];
	asector = Rest[#] -> #[[1]] & /@ elt[[1]];
	dsector = Rest[#] -> -Conjugate[#[[1]]] & /@ elt[[2]];
	diag = Table[{i, i} -> 
		If[i > np, -Conjugate[modes[[i]]], modes[[i]]], {i, 1, np + nq}];
	sym[s_, l_] := 
	  Join[Reverse[#[[1]]] -> s*Conjugate[#[[2]]] & /@ l, l];
	Join[sym[1, asector], sym[-1, bsector], sym[1, dsector], diag] // 
	  SparseArray // Normal
]

ShowGraph[lm__] := Module[
	{lms,graphEdges,graphLabels},
	lms = SparseArray[lm];
	graphEdges = #[[1]] -> #[[2]] & /@ lms["NonzeroPositions"];
	graphLabels = TraditionalForm /@ lms["NonzeroValues"];
	GraphPlot[Transpose[{graphEdges, graphLabels}], DirectedEdges -> True, PackingMethod -> "LayeredLeft"]
]

ShowGraphSimplified[lm__] := Module[{map, f, lms, p},
  lms = SparseArray[lm];
  
  f[{x_}] := x;
  f[{x_, y_} /; x == Conjugate[y]] := x;
  f[{x_, y_} /; x == -Conjugate[y]] := Framed[x];
  
  map = GroupBy[
    Transpose@{Sort /@ lms["NonzeroPositions"], lms["NonzeroValues"]},
     First -> Last];
  p = f /@ map;
  GraphPlot[Transpose@{Rule @@@ Keys[p], Values[p]}]]

ConnectedModes[lm__] := Sort /@ WeaklyConnectedComponents[Rule @@@ SparseArray[lm]["NonzeroPositions"]]

ExtractModes[ci__,em__] := Module[
	{np,nq,modes,conversion,gain,npn,npq,nmodes,map,rw,nconv,ngain},
	
	{np,nq,modes,conversion,gain} = ci;

	{npn, npq} = {Count[em, x_ /; x <= np], Count[em, x_ /; x > np]};
	nmodes = modes[[em]];
	map = Table[em[[i]] -> i, {i, 1, Length[em]}];
	rw[{a_, b_, c_}] := {a, b /. map, c /. map};
	nconv = rw /@ Select[conversion, MemberQ[em, #[[2]]] && MemberQ[em, #[[3]]] &];
	ngain = rw /@ Select[gain, MemberQ[em, #[[2]]] && MemberQ[em, #[[3]]] &];
	{npn, npq, nmodes, nconv, ngain}
]

Cofactor[a_, i_, j_] := Det[Drop[a, {i}, {j}]]*(2*Mod[i*j, 2] - 1)

End[]

EndPackage[]

