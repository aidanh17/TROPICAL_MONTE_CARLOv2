(* ============================================================================
   bench_lift_variance.wl  --  §8.3 Benchmark: lifted vs unlifted variance

   Test 6 Cases A, B, C at fixed 10^6 MC samples.
   Unlifted: EvaluateTropicalMC with unit-coefficient proxy fan.
   Lifted:   EvaluateTropicalMCLifted with explicit k-rules.
   Output:   per-case tables to stdout + SANDBOX/benchmark_results.md.

   Case A: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2}
           Lift: explicit k=2 rule on 10^6 x1^2.  Automatic fan.
   Case B: P = 10^-4 + 10^4 x1^2 + 10^-4 x2^2 + 10^4 x1 x2^2 + x1^2 x2,  B={-2}
           (i) two-rule k=1: primary {2,0} + secondary {1,2} with residual.
           (ii) single-rule k=2 on primary {2,0}.
   Case C: P = 1 + 10^8 x1^3 x2 + x2^3, B={-3}
           Degenerate lifted polytope (3 monomials in 3D lift).
           Explicit fan provided for both k=2 (z0=10^4) and k=4 (z0=10^2).
           The automatic fan detection fires liftdegenerate for both k values.

   Run:
     cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file SANDBOX/bench_lift_variance.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Print["tropical_eval.wl loaded"];
Print[];

NSamples = 1000000;


(* ============================================================================
   Helper: run unlifted EvaluateTropicalMC
   Returns result Association.  Variable names use camelCase (no underscores
   to avoid Mathematica pattern-matching issues).
   ============================================================================ *)

runUnlifted[spec_Association, proxyPoly_, caseLabel_String] :=
Module[
  {vars, verts, fan, timing, res, direct,
   value, err, sigma, allSD, cfmAll},

  vars  = spec["Variables"];
  verts = PolytopeVertices[proxyPoly^(-1), vars];
  fan   = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["  ", caseLabel, " unlifted: ", Length[fan[[1]]], " rays, ",
        Length[fan[[2]]], " sectors"];

  direct = Quiet @ NIntegrate[
    (Times @@ MapThread[Power, {spec["Polynomials"], spec["PolynomialExponents"]}]) *
    (Times @@ MapThread[Power, {vars, spec["MonomialExponents"]}]),
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 20, PrecisionGoal -> 5
  ];
  Print["  NIntegrate reference = ", direct];

  {timing, res} = AbsoluteTiming[
    Quiet @ EvaluateTropicalMC[spec, fan, {{}},
      "NSamples"   -> NSamples,
      "RunChecks"  -> False,
      "Verbose"    -> False]
  ];

  If[!AssociationQ[res] || !KeyExistsQ[res, "Results"],
    Print["  ", caseLabel, " unlifted FAILED: ", res];
    Return[<|"caseLabel" -> caseLabel, "mode" -> "unlifted", "failed" -> True,
             "direct" -> direct, "fan" -> fan|>]
  ];

  value  = res["Results"][[1]]["Re"];
  err    = res["Results"][[1]]["ReErr"];
  sigma  = err * Sqrt[NSamples];

  (* CFM on first few sectors *)
  allSD = Table[
    ProcessSector[spec, fan[[1]], fan[[2, s]], s],
    {s, Length[fan[[2]]]}
  ];
  cfmAll = Table[
    If[!TrueQ[sd["IsDivergent"]],
      CheckFlatteningMagnitude[sd, 20, {}],
      <|"Mean" -> 0, "Max" -> 0, "Min" -> 0|>],
    {sd, allSD}
  ];

  Print["  Unlifted result: ", value, " +/- ", err,
        "  sigma=", sigma, "  relErr=", N[Abs[(value-direct)/direct], 3]];
  Print["  Wall-clock: ", timing, " s"];
  Print["  CFM min=", Min[#["Min"] & /@ cfmAll], "  max=", Max[#["Max"] & /@ cfmAll]];

  <|"caseLabel" -> caseLabel, "mode" -> "unlifted", "failed" -> False,
    "value" -> value, "err" -> err, "sigma" -> sigma,
    "direct" -> direct, "relErr" -> Abs[(value-direct)/direct],
    "wallClock" -> timing, "nSectors" -> Length[fan[[2]]],
    "fan" -> fan,
    "cfmMin" -> Min[#["Min"] & /@ cfmAll],
    "cfmMax" -> Max[#["Max"] & /@ cfmAll],
    "cfmMean" -> Mean[#["Mean"] & /@ cfmAll]
  |>
];


(* ============================================================================
   Helper: run lifted EvaluateTropicalMCLifted
   ============================================================================ *)

runLifted[spec_Association, liftRules_, fanDataOpt_, caseLabel_String, liftLabel_String,
          direct_] :=
Module[
  {timing, res, liftRes, liftedSpec, liftData, liftedFan,
   value, err, sigma, nSects, nDropped,
   allLiftedSD, cfmAll, hasConstList},

  Print["  ", caseLabel, " ", liftLabel, " (", Length[liftRules],
        " rule(s)): running ..."];

  {timing, res} = AbsoluteTiming[
    Quiet @ EvaluateTropicalMCLifted[spec, {{}},
      "LiftRules"  -> liftRules,
      "FanData"    -> fanDataOpt,
      "NSamples"   -> NSamples,
      "RunChecks"  -> False,
      "Verbose"    -> True]
  ];

  If[!AssociationQ[res] || !KeyExistsQ[res, "Results"],
    Print["  ", caseLabel, " ", liftLabel, " FAILED: ", res];
    Return[<|"caseLabel" -> caseLabel, "liftLabel" -> liftLabel,
             "mode" -> "lifted", "failed" -> True|>]
  ];

  value  = res["Results"][[1]]["Re"];
  err    = res["Results"][[1]]["ReErr"];
  sigma  = err * Sqrt[NSamples];
  nSects = Lookup[res, "ConvergentSectors", 0];

  (* Re-process sectors to collect HasConstantTerm and CFM *)
  liftRes = LiftCoefficients[spec, liftRules];
  If[AssociationQ[liftRes],
    liftedSpec = liftRes["LiftedSpec"];
    liftData   = liftRes["LiftData"];
    liftedFan  = If[fanDataOpt =!= Automatic, fanDataOpt,
      Module[{vs = Quiet @ PolytopeVertices[
        (Times @@ liftedSpec["Polynomials"])^(-1), liftedSpec["Variables"]]},
        Quiet @ ComputeDecomposition[vs, "ShowProgress" -> False]
      ]
    ];
    allLiftedSD = If[ListQ[liftedFan] && Length[liftedFan] >= 2,
      Table[
        Quiet @ ProcessSectorLifted[liftedSpec, liftedFan[[1]],
                                    liftedFan[[2, s]], s, liftData],
        {s, Length[liftedFan[[2]]]}
      ],
      {}
    ];
    cfmAll = Table[
      If[AssociationQ[sd] && !TrueQ[Lookup[sd, "EmptyDomain", False]],
        Quiet @ CheckFlatteningMagnitude[sd, 20, {}],
        <|"Mean" -> 0, "Max" -> 0, "Min" -> 0|>],
      {sd, allLiftedSD}
    ];
    hasConstList = Table[
      If[AssociationQ[sd] && !TrueQ[Lookup[sd, "EmptyDomain", False]],
        Lookup[sd, "HasConstantTerm", Missing[]],
        "dropped"],
      {sd, allLiftedSD}
    ];
    nDropped = Count[allLiftedSD,
      a_ /; AssociationQ[a] && TrueQ[Lookup[a, "EmptyDomain", False]]];
    ,
    cfmAll = {}; hasConstList = {}; nDropped = 0
  ];

  Print["  Lifted result: ", value, " +/- ", err,
        "  sigma=", sigma, "  relErr=", N[Abs[(value-direct)/direct], 3]];
  Print["  Wall-clock: ", timing, " s"];
  Print["  ConvergentSectors=", nSects, "  EmptyDomainDrops=", nDropped];
  If[Length[cfmAll] > 0,
    Print["  CFM min=", Min[Cases[#["Min"] & /@ cfmAll, _?NumberQ]],
          "  max=", Max[Cases[#["Max"] & /@ cfmAll, _?NumberQ]]]
  ];
  If[MemberQ[hasConstList, False],
    Print["  WARNING: HasConstantTerm=False in sector(s) — heavy-tail risk"]
  ];

  <|"caseLabel"  -> caseLabel,
    "liftLabel"  -> liftLabel,
    "mode"       -> "lifted",
    "failed"     -> False,
    "value"      -> value,
    "err"        -> err,
    "sigma"      -> sigma,
    "direct"     -> direct,
    "relErr"     -> Abs[(value-direct)/direct],
    "wallClock"  -> timing,
    "nSectors"   -> nSects,
    "nDropped"   -> nDropped,
    "cfmMin"     -> If[Length[cfmAll]>0, Min[Cases[#["Min"] & /@ cfmAll, _?NumberQ]], 0],
    "cfmMax"     -> If[Length[cfmAll]>0, Max[Cases[#["Max"] & /@ cfmAll, _?NumberQ]], 0],
    "cfmMean"    -> If[Length[cfmAll]>0,
                     With[{v=Cases[#["Mean"] & /@ cfmAll, _?NumberQ]},
                          If[Length[v]>0, Mean[v], 0]], 0],
    "hasConstTermAnyFalse" -> MemberQ[hasConstList, False],
    "hasConstList"         -> hasConstList
  |>
];

vrFactor[sigU_, sigL_] := If[NumericQ[sigU] && NumericQ[sigL] && sigL > 0,
  (sigU/sigL)^2, "N/A"];

fmtN[x_] := If[NumericQ[x], ToString[N[x, 4]], ToString[x]];


(* ============================================================================
   CASE A
   ============================================================================ *)

Print[""];
Print["================================================================="];
Print["CASE A: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2, B={-2}"];
Print["================================================================="];

specA = <|
  "Polynomials"         -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

resAU = runUnlifted[specA, 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2, "CaseA"];
Print["--- Case A unlifted done ---"];

resALk2 = runLifted[specA,
  {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>},
  Automatic, "CaseA", "lifted-k2", resAU["direct"]];
Print["--- Case A lifted k=2 done ---"];

vrA = vrFactor[resAU["sigma"], resALk2["sigma"]];
Print["Case A variance reduction (sigma_U/sigma_L)^2 = ", N[vrA, 4]];
Print["Case A >= 10x PASS/FAIL: ", If[NumericQ[vrA] && vrA >= 10, "PASS", "FAIL"]];
Print[];


(* ============================================================================
   CASE B
   ============================================================================ *)

Print[""];
Print["================================================================="];
Print["CASE B: P = 10^-4 + 10^4 x1^2 + 10^-4 x2^2 + 10^4 x1 x2^2 + x1^2 x2, B={-2}"];
Print["================================================================="];

specB = <|
  "Polynomials"         -> {10^-4 + 10^4 x[1]^2 + 10^-4 x[2]^2 +
                            10^4 x[1] x[2]^2 + x[1]^2 x[2]},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

directB = Quiet @ NIntegrate[
  1 / (10^-4 + 10^4 t1^2 + 10^-4 t2^2 + 10^4 t1 t2^2 + t1^2 t2)^2,
  {t1, 0, Infinity}, {t2, 0, Infinity},
  MaxRecursion -> 20, PrecisionGoal -> 5
];
Print["  NIntegrate reference = ", directB];

resBU = runUnlifted[specB,
  10^-4 + 10^4 x[1]^2 + 10^-4 x[2]^2 + 10^4 x[1] x[2]^2 + x[1]^2 x[2],
  "CaseB"];
Print["--- Case B unlifted done ---"];

(* Two-rule k=1: primary {2,0} + secondary {1,2} with residual coefficient *)
Print["Trying two-rule k=1 lift ..."];
resBL2rule = runLifted[specB,
  {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->1|>,
   <|"PolyIndex"->1, "ExponentVector"->{1,2}, "k"->1|>},
  Automatic, "CaseB", "2rule-k1", directB];
Print["--- Case B 2-rule k=1 done ---"];

(* Single-rule k=2 on primary *)
Print["Trying single-rule k=2 lift ..."];
resBL1k2 = runLifted[specB,
  {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>},
  Automatic, "CaseB", "1rule-k2", directB];
Print["--- Case B 1-rule k=2 done ---"];

vrB2rule = vrFactor[resBU["sigma"],
  If[!TrueQ[resBL2rule["failed"]], resBL2rule["sigma"], Infinity]];
vrB1k2   = vrFactor[resBU["sigma"],
  If[!TrueQ[resBL1k2["failed"]], resBL1k2["sigma"], Infinity]];
Print["Case B variance reductions: 2-rule-k1=", N[vrB2rule,4],
      "  1-rule-k2=", N[vrB1k2,4]];
Print["Case B is report-only (no >= 10x gate)"];
Print[];


(* ============================================================================
   CASE C
   ============================================================================ *)

Print[""];
Print["================================================================="];
Print["CASE C: P = 1 + 10^8 x1^3 x2 + x2^3, B={-3}"];
Print["================================================================="];

specC = <|
  "Polynomials"         -> {1 + 10^8 x[1]^3 x[2] + x[2]^3},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-3},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

resCU = runUnlifted[specC, 1 + x[1]^3 x[2] + x[2]^3, "CaseC"];
Print["--- Case C unlifted done ---"];

(* Both k=2 and k=4 have degenerate lifted polytopes (3 monomials in 3D).
   The automatic fan fires liftdegenerate.  Provide explicit fans.
   Fan structure for Case C: 5 rays + 6 simplices (same as Toy-1 fan form
   adapted for the lineality direction of each k).
   k=2: lineality direction (2,0,-3); k=4: lineality direction (4,0,-3).
   Validated via ValidateLiftedDecomposition (RelErr < 1% confirmed). *)
dvCk2 = {{1,0,0},{0,1,0},{-1,-3,0},{2,0,-3},{-2,0,3}};
slC   = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}};
dvCk4 = {{1,0,0},{0,1,0},{-1,-3,0},{4,0,-3},{-4,0,3}};

fanCk2 = {dvCk2, slC};
fanCk4 = {dvCk4, slC};

Print["Trying k=2 (z0=10^4) with explicit fan ..."];
resCLk2 = runLifted[specC,
  {<|"PolyIndex"->1, "ExponentVector"->{3,1}, "k"->2|>},
  fanCk2, "CaseC", "k2-explicit", resCU["direct"]];
Print["--- Case C lifted k=2 done ---"];

Print["Trying k=4 (z0=10^2) with explicit fan ..."];
resCLk4 = runLifted[specC,
  {<|"PolyIndex"->1, "ExponentVector"->{3,1}, "k"->4|>},
  fanCk4, "CaseC", "k4-explicit", resCU["direct"]];
Print["--- Case C lifted k=4 done ---"];

vrCk2 = vrFactor[resCU["sigma"],
  If[!TrueQ[resCLk2["failed"]], resCLk2["sigma"], Infinity]];
vrCk4 = vrFactor[resCU["sigma"],
  If[!TrueQ[resCLk4["failed"]], resCLk4["sigma"], Infinity]];
vrCbest = Max[{vrCk2, vrCk4} /. {"N/A" -> 0}];

Print["Case C variance reductions: k=2=", N[vrCk2,4], "  k=4=", N[vrCk4,4]];
Print["Case C best variance reduction = ", N[vrCbest,4]];
Print["Case C >= 10x PASS/FAIL: ", If[NumericQ[vrCbest] && vrCbest >= 10, "PASS", "FAIL"]];
Print[];


(* ============================================================================
   WRITE benchmark_results.md
   ============================================================================ *)

mdFile = FileNameJoin[{DirectoryName[$InputFileName], "benchmark_results.md"}];

Module[{lines, sigUA, sigLA, sigUB, sigLBbest, vrBbest,
        sigUC, sigLCbest, hcA, hcB, hcC},

  sigUA = resAU["sigma"];
  sigLA = If[!TrueQ[resALk2["failed"]], resALk2["sigma"], Infinity];

  sigUB = resBU["sigma"];
  sigLBbest = Min[{
    If[!TrueQ[resBL2rule["failed"]], resBL2rule["sigma"], Infinity],
    If[!TrueQ[resBL1k2["failed"]], resBL1k2["sigma"], Infinity]
  }];
  vrBbest = vrFactor[sigUB, sigLBbest];

  sigUC = resCU["sigma"];
  sigLCbest = Min[{
    If[!TrueQ[resCLk2["failed"]], resCLk2["sigma"], Infinity],
    If[!TrueQ[resCLk4["failed"]], resCLk4["sigma"], Infinity]
  }];
  vrCbest2 = vrFactor[sigUC, sigLCbest];

  hcA = TrueQ[Lookup[resALk2, "hasConstTermAnyFalse", False]];
  hcB = TrueQ[Lookup[resBL2rule, "hasConstTermAnyFalse", False]] ||
        TrueQ[Lookup[resBL1k2,   "hasConstTermAnyFalse", False]];
  hcC = TrueQ[Lookup[resCLk2, "hasConstTermAnyFalse", False]] ||
        TrueQ[Lookup[resCLk4, "hasConstTermAnyFalse", False]];

  lines = {
    "# TROPICAL_MONTE_CARLO v2: Lift-Variance Benchmark",
    "",
    "**Date:** 2026-06-12  ",
    "**Samples per run:** 1,000,000  ",
    "**Caveat:** Wall-clock includes C++ compilation and optional NIntegrate; MC-loop-only time not separately isolated by the driver.  Timings are total AbsoluteTiming of the driver call.",
    "",
    "---",
    "",
    "## Case A: `P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2`, B={-2}",
    "",
    "Lift rule: explicit `{PolyIndex=1, ExponentVector={2,0}, k=2}`.  Automatic fan (Polymake, full-dimensional 3D Newton polytope).  Also ran Automatic k=1 in test_lifted.wl (informational; see §23B).",
    "",
    "| Mode | MC Value | StdErr | Per-sample σ | NIntegrate | Rel.Dev. | Sectors | Dropped | Wall(s) |",
    "|------|----------|--------|--------------|------------|----------|---------|---------|---------|",
    "| Unlifted | " <> fmtN[resAU["value"]] <>
      " | " <> fmtN[resAU["err"]] <>
      " | " <> fmtN[sigUA] <>
      " | " <> fmtN[resAU["direct"]] <>
      " | " <> fmtN[resAU["relErr"]] <>
      " | " <> ToString[resAU["nSectors"]] <>
      " | 0 | " <> fmtN[resAU["wallClock"]] <> " |",
    "| Lifted k=2 | " <>
      If[!TrueQ[resALk2["failed"]],
        fmtN[resALk2["value"]] <> " | " <> fmtN[resALk2["err"]] <>
        " | " <> fmtN[sigLA] <>
        " | " <> fmtN[resALk2["direct"]] <>
        " | " <> fmtN[resALk2["relErr"]] <>
        " | " <> ToString[resALk2["nSectors"]] <>
        " | " <> ToString[resALk2["nDropped"]] <>
        " | " <> fmtN[resALk2["wallClock"]],
        "FAILED | - | - | - | - | - | - | -"] <> " |",
    "",
    "**Variance-reduction factor = " <> fmtN[vrA] <> "**  -->  " <>
      If[NumericQ[vrA] && vrA >= 10, "**PASS** (>= 10x)", "**FAIL** (< 10x)"],
    "CFM (20 samples/sector): Unlifted min/max = " <>
      fmtN[resAU["cfmMin"]] <> " / " <> fmtN[resAU["cfmMax"]] <>
      "; Lifted k=2 min/max = " <>
      If[!TrueQ[resALk2["failed"]], fmtN[resALk2["cfmMin"]] <> " / " <> fmtN[resALk2["cfmMax"]], "N/A"] <> ".",
    If[hcA, "**HasConstantTerm=False** in at least one lifted sector.", "HasConstantTerm=True in all lifted sectors."],
    "",
    "---",
    "",
    "## Case B: `P = 10^-4 + 10^4 x1^2 + 10^-4 x2^2 + 10^4 x1 x2^2 + x1^2 x2`, B={-2}",
    "",
    "Two lift strategies tried: **(i)** two-rule k=1 (primary {2,0} + secondary {1,2} with residual); **(ii)** single-rule k=2 on primary {2,0}.  Case B is report-only.",
    "",
    "| Mode | MC Value | Per-sample σ | NIntegrate | Rel.Dev. | Sectors | Dropped | Var.Red. | Wall(s) |",
    "|------|----------|--------------|------------|----------|---------|---------|----------|---------|",
    "| Unlifted | " <> fmtN[resBU["value"]] <>
      " | " <> fmtN[sigUB] <>
      " | " <> fmtN[directB] <>
      " | " <> fmtN[resBU["relErr"]] <>
      " | " <> ToString[resBU["nSectors"]] <>
      " | 0 | - | " <> fmtN[resBU["wallClock"]] <> " |",
    "| Lifted 2-rule k=1 | " <>
      If[!TrueQ[resBL2rule["failed"]],
        fmtN[resBL2rule["value"]] <> " | " <>
        fmtN[resBL2rule["sigma"]] <> " | " <> fmtN[directB] <>
        " | " <> fmtN[resBL2rule["relErr"]] <>
        " | " <> ToString[resBL2rule["nSectors"]] <>
        " | " <> ToString[resBL2rule["nDropped"]] <>
        " | " <> fmtN[vrFactor[sigUB, resBL2rule["sigma"]]] <>
        " | " <> fmtN[resBL2rule["wallClock"]],
        "FAILED | - | " <> fmtN[directB] <> " | - | - | - | - | -"] <> " |",
    "| Lifted 1-rule k=2 | " <>
      If[!TrueQ[resBL1k2["failed"]],
        fmtN[resBL1k2["value"]] <> " | " <>
        fmtN[resBL1k2["sigma"]] <> " | " <> fmtN[directB] <>
        " | " <> fmtN[resBL1k2["relErr"]] <>
        " | " <> ToString[resBL1k2["nSectors"]] <>
        " | " <> ToString[resBL1k2["nDropped"]] <>
        " | " <> fmtN[vrFactor[sigUB, resBL1k2["sigma"]]] <>
        " | " <> fmtN[resBL1k2["wallClock"]],
        "FAILED | - | " <> fmtN[directB] <> " | - | - | - | - | -"] <> " |",
    "",
    If[hcB, "**HasConstantTerm=False** in at least one Case-B lifted sector.", "HasConstantTerm=True in all Case-B lifted sectors."],
    "Case B is report-only (no >= 10x success criterion).",
    "",
    "---",
    "",
    "## Case C: `P = 1 + 10^8 x1^3*x2 + x2^3`, B={-3}",
    "",
    "The lifted Newton polytope (3 monomials in 3D) is lower-dimensional for both k=2 and k=4.  Automatic fan detection fires `liftdegenerate`.  **Explicit fans supplied** (5 rays: {1,0,0},{0,1,0},{-1,-3,0} plus the lineality pair {2,0,-3}/{-2,0,3} for k=2 and {4,0,-3}/{-4,0,3} for k=4; 6 simplices of the Toy-1 pattern).",
    "",
    "| Mode | MC Value | Per-sample σ | NIntegrate | Rel.Dev. | Sectors | Dropped | Var.Red. | Wall(s) |",
    "|------|----------|--------------|------------|----------|---------|---------|----------|---------|",
    "| Unlifted | " <> fmtN[resCU["value"]] <>
      " | " <> fmtN[sigUC] <>
      " | " <> fmtN[resCU["direct"]] <>
      " | " <> fmtN[resCU["relErr"]] <>
      " | " <> ToString[resCU["nSectors"]] <>
      " | 0 | - | " <> fmtN[resCU["wallClock"]] <> " |",
    "| Lifted k=2 | " <>
      If[!TrueQ[resCLk2["failed"]],
        fmtN[resCLk2["value"]] <> " | " <>
        fmtN[resCLk2["sigma"]] <> " | " <> fmtN[resCU["direct"]] <>
        " | " <> fmtN[resCLk2["relErr"]] <>
        " | " <> ToString[resCLk2["nSectors"]] <>
        " | " <> ToString[resCLk2["nDropped"]] <>
        " | " <> fmtN[vrCk2] <>
        " | " <> fmtN[resCLk2["wallClock"]],
        "FAILED | - | " <> fmtN[resCU["direct"]] <> " | - | - | - | - | -"] <> " |",
    "| Lifted k=4 | " <>
      If[!TrueQ[resCLk4["failed"]],
        fmtN[resCLk4["value"]] <> " | " <>
        fmtN[resCLk4["sigma"]] <> " | " <> fmtN[resCU["direct"]] <>
        " | " <> fmtN[resCLk4["relErr"]] <>
        " | " <> ToString[resCLk4["nSectors"]] <>
        " | " <> ToString[resCLk4["nDropped"]] <>
        " | " <> fmtN[vrCk4] <>
        " | " <> fmtN[resCLk4["wallClock"]],
        "FAILED | - | " <> fmtN[resCU["direct"]] <> " | - | - | - | - | -"] <> " |",
    "",
    "**Best variance-reduction = " <> fmtN[vrCbest2] <> "**  -->  " <>
      If[NumericQ[vrCbest2] && vrCbest2 >= 10, "**PASS** (>= 10x)", "**FAIL** (< 10x)"],
    If[hcC, "**HasConstantTerm=False** in at least one Case-C lifted sector.", "HasConstantTerm=True in all Case-C lifted sectors."],
    "",
    "---",
    "",
    "## Summary Table",
    "",
    "| Case | sigma_unlifted | sigma_lifted (best) | Var.reduction | >= 10x? |",
    "|------|---------------|---------------------|---------------|---------|",
    "| A    | " <> fmtN[sigUA] <> " | " <> fmtN[sigLA] <>
      " | " <> fmtN[vrA] <>
      " | " <> If[NumericQ[vrA] && vrA >= 10, "PASS", "FAIL"] <> " |",
    "| B    | " <> fmtN[sigUB] <> " | " <> fmtN[sigLBbest] <>
      " | " <> fmtN[vrBbest] <>
      " | report-only |",
    "| C    | " <> fmtN[sigUC] <> " | " <> fmtN[sigLCbest] <>
      " | " <> fmtN[vrCbest2] <>
      " | " <> If[NumericQ[vrCbest2] && vrCbest2 >= 10, "PASS", "FAIL"] <> " |",
    "",
    "---",
    "",
    "## Important Caveats on Variance Estimation",
    "",
    "**Sample-based sigma can severely underestimate true sigma for heavy-tailed integrands.**",
    "See `SANDBOX/phase1_results.txt` for a full true-sigma analysis: the unlifted Toy-1 MC",
    "missed the spike mass entirely (sample means ~927 vs true sum 5e5); true-sigma analysis",
    "found infinite-variance configurations in BOTH lifted and unlifted Toy-1 sectors.",
    "For Toy 2, one lifted sector had divergent I2 (0.2% contribution); excluding it, the",
    "remaining 5 sectors showed 4.2x true-sigma improvement.",
    "",
    "**When HasConstantTerm = False** for a lifted sector (re-clearing loses the constant term),",
    "the sector's polynomial Qtilde vanishes at a cube corner, making Qtilde^B non-integrable",
    "in L^2 for B << 0.  The sample-based sigma (from finite MC) is then an optimistic",
    "(possibly dramatically underestimated) proxy for the true variance.  Per-case status:",
    "",
    "- Case A lifted sectors: " <>
      If[hcA, "SOME have HasConstantTerm=False. Sample sigma may be OPTIMISTIC.", "All HasConstantTerm=True."],
    "- Case B lifted sectors: " <>
      If[hcB, "SOME have HasConstantTerm=False. Sample sigma may be OPTIMISTIC.", "All HasConstantTerm=True."],
    "- Case C lifted sectors: " <>
      If[hcC, "SOME have HasConstantTerm=False. Sample sigma may be OPTIMISTIC.", "All HasConstantTerm=True."]
  };

  Export[mdFile, StringRiffle[Select[lines, StringQ], "\n"], "Text"];
  Print["Wrote ", mdFile];

  Print[""];
  Print["================================================================="];
  Print["  BENCHMARK SUMMARY"];
  Print["================================================================="];
  Print["  Case A: sigma_unlifted=", fmtN[sigUA], "  sigma_lifted(k=2)=", fmtN[sigLA]];
  Print["          Var.reduction=", fmtN[vrA],
        "  [", If[NumericQ[vrA] && vrA >= 10, "PASS >= 10x", "FAIL < 10x"], "]"];
  Print["  Case B (report-only): sigma_unlifted=", fmtN[sigUB],
        "  best=", fmtN[sigLBbest], "  Var.red.=", fmtN[vrBbest]];
  Print["  Case C: sigma_unlifted=", fmtN[sigUC],
        "  k2=", fmtN[If[!TrueQ[resCLk2["failed"]],resCLk2["sigma"],Infinity]],
        "  k4=", fmtN[If[!TrueQ[resCLk4["failed"]],resCLk4["sigma"],Infinity]]];
  Print["          Best Var.reduction=", fmtN[vrCbest2],
        "  [", If[NumericQ[vrCbest2] && vrCbest2 >= 10, "PASS >= 10x", "FAIL < 10x"], "]"];
  Print["  HasConstantTerm=False warnings: A=", hcA, "  B=", hcB, "  C=", hcC];
  Print["================================================================="];
];
