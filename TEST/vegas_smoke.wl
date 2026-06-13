(* R2 VEGAS smoke + R3 CUBA-absent test, driven through EvaluateTropicalMC. *)
SetDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{ParentDirectory[], "tropical_eval.wl"}]];
interfiles = FileNameJoin[{Directory[], "INTERFILES"}];
If[!DirectoryQ[interfiles], CreateDirectory[interfiles]];

(* n=2 simplex: P = 1 + x1 + x2, A=0, B=-4, exact = 1/(n+1)! = 1/6 *)
spec = <|"Polynomials" -> {1 + x[1] + x[2]}, "MonomialExponents" -> {0, 0},
  "PolynomialExponents" -> {-4}, "Variables" -> {x[1], x[2]},
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
verts = PolytopeVertices[(Times @@ spec["Polynomials"])^(-1), spec["Variables"]];
fan = ComputeDecomposition[verts, "ShowProgress" -> False];
exact = 1/6.;

Print["===== R2: VEGAS smoke (n=2 simplex, exact = ", exact, ") ====="];
rVE = EvaluateTropicalMC[spec, fan, {{}},
  "Integrator" -> "VEGAS", "NSamples" -> 200000,
  "RunChecks" -> False, "Verbose" -> True,
  "WorkingDirectory" -> interfiles];
If[AssociationQ[rVE],
  Module[{r = rVE["Results"][[1]], val, finite},
    val = r["Re"] + I r["Im"];
    finite = AllTrue[{r["Re"], r["Im"], r["ReErr"], r["ImErr"]}, NumericQ[#] && Abs[#] < 10^30 &];
    Print["  VEGAS = ", val, " +/- (", r["ReErr"], ", ", r["ImErr"], ")"];
    Print["  |rel dev| from exact = ", Abs[(r["Re"] - exact)/exact]];
    Print["  4 finite numbers: ", finite];
    Print["  R2 ", If[finite && Abs[(r["Re"] - exact)/exact] < 1*^-2, "PASS", "FAIL"]];
  ],
  Print["  R2 FAIL: EvaluateTropicalMC[VEGAS] returned ", rVE]];

Print[];
Print["===== R3: CUBA-absent => nocuba + $Failed (no silent MC fallback) ====="];
(* Override the private detector to simulate no CUBA. *)
Module[{old = DownValues[TropicalEval`Private`findCubaPrefix], got, msgs},
  TropicalEval`Private`findCubaPrefix[] := $Failed;
  {got, msgs} = Reap[
    Quiet[
      Check[
        EvaluateTropicalMC[spec, fan, {{}}, "Integrator" -> "VEGAS",
          "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> interfiles],
        "MSG_FIRED", TropicalEval::nocuba],
      TropicalEval::nocuba]];
  DownValues[TropicalEval`Private`findCubaPrefix] = old;
  Print["  returned: ", got, "   (expect MSG_FIRED, i.e. nocuba fired and $Failed path taken)"];
  Print["  R3 ", If[got === "MSG_FIRED", "PASS", "FAIL"]];
];

Print[];
Print["===== bad integrator => badintegrator + $Failed ====="];
Module[{got},
  got = Quiet[Check[
    EvaluateTropicalMC[spec, fan, {{}}, "Integrator" -> "NOPE",
      "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> interfiles],
    "MSG_FIRED", TropicalEval::badintegrator],
    TropicalEval::badintegrator];
  Print["  returned: ", got];
  Print["  badintegrator ", If[got === "MSG_FIRED", "PASS", "FAIL"]];
];
