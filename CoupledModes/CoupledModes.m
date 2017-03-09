(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Mar 8, 2017 *)

BeginPackage["CoupledModes`"]
(* Exported symbols added here with SymbolName::usage *) 

CoupledModeSystem::usage =
	"Gives a langevin matrix"

Begin["`Private`"]
(* Implementation of the package *)

CoupledModeSystem[np_Integer,nq_Integer,modes__,conversion__,gain__] := Module[
	{coordify, bsector,cprocs,asector,dsector,diag,sym},

	coordify[i_, j_] := {Min[i, j], Max[i, j]};
	bsector = coordify[#[[2]], #[[3]]] -> #[[1]] & /@ gain;
	
	cprocs = GroupBy[conversion, ((# > np &) /@ # &)@*Rest];
	asector = {#[[2]], #[[3]]} -> #[[1]] & /@ cprocs[{False, False}];
	dsector = {#[[2]], #[[3]]} -> -Conjugate[#[[1]]] & /@ cprocs[{True, True}];
	diag = Table[{i, i} -> If[i > np, -Conjugate[modes[[i]]], modes[[i]]], {i, 1, np + nq}];
	
	sym[s_, l_] := Join[Reverse[#[[1]]] -> s*Conjugate[#[[2]]] & /@ l, l];
	
	Join[sym[1, asector], sym[-1, bsector], sym[1, dsector], diag] // SparseArray // Normal
]

End[]

EndPackage[]

