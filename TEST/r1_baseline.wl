(* R1 regression helper: generate the MC .cpp for representative specs and
   save them under TEST/INTERFILES/r1_<tag>_<phase>.cpp so a before/after diff
   proves the default (MC) codegen is byte-identical.  Usage:
     wolframscript -file TEST/r1_baseline.wl before
     wolframscript -file TEST/r1_baseline.wl after
*)
SetDirectory[DirectoryName[$InputFileName]];
phase = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "before"];
Get[FileNameJoin[{ParentDirectory[], "tropical_eval.wl"}]];
outDir = FileNameJoin[{Directory[], "INTERFILES"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

genMC[tag_, spec_] := Module[{verts, fan, dv, sl, conv, f},
  verts = PolytopeVertices[(Times @@ spec["Polynomials"])^(-1), spec["Variables"]];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dv, sl} = fan;
  conv = Select[
    Table[ProcessSector[spec, dv, sl[[s]], s, "Verbose" -> False], {s, Length[sl]}],
    (AssociationQ[#] && !#["IsDivergent"]) &];
  f = FileNameJoin[{outDir, "r1_" <> tag <> "_" <> phase <> ".cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, f];   (* default Integrator -> "MC" *)
  Print["wrote ", f, "  (", Length[conv], " sectors)"];
];

(* param-free, real exponent (n=2 simplex) *)
genMC["simplex2",
  <|"Polynomials" -> {1 + x[1] + x[2]}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-4}, "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>];

(* complex exponent, kinematic symbol (Example 17) *)
genMC["ex17",
  <|"Polynomials" -> {1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-(2 + I/2)},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {lam},
    "RegulatorSymbol" -> None|>];

Print["R1 ", phase, " snapshots done."];
