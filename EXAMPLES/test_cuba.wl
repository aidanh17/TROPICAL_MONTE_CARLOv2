(* ============================================================================
   test_cuba.wl — CUBA cross-check of ALL test-suite integrals

   Integrates every convergent integral appearing in the validation suite
   (tropical_eval_examples.wl Tests 1, 2, 5, 6, 7) and in the lifted-
   pipeline tests (test_lifted.wl 23A/23D) with the CUBA library, applied
   DIRECTLY to the original integrand on [0,inf)^n (compactification
   x = t/(1-t); no tropical decomposition).  For each integral it compares:

       CUBA Cuhre        deterministic adaptive cubature   (reference)
       sector sum        tropical decomposition, per-sector NIntegrate
       tropical MC       the C++ Monte Carlo pipeline
       exact value       where a closed form exists (23A, 23D)

   Test 3v2 is a pair of negative tests (error guards on non-convergent
   specs) — nothing to integrate.  Test 23B's integral is identical to
   Test 6 Case A, which is covered (including its lifted MC).

   PASS gates:
     - sector sum vs Cuhre:  case-dependent tolerance (1%-5%, matching
       the precision goals the suite itself uses)
     - tropical MC vs Cuhre: within 5 sigma of the MC error bar; cases
       whose cleared polynomial has a SMALL constant term (Test 6B) are
       heavy-tailed for plain uniform MC and are reported without a gate
       (see Example 20 / SUMMARY §8)
     - exact value vs Cuhre and vs sector sum where available

   Requires: tropical_fan.wl, tropical_eval.wl, Polymake, g++, CUBA
   (`brew install cuba`).  Runtime: roughly 5-15 minutes (the 4D case
   and the 10^7-eval Cuhre runs dominate).

   Run:
     wolframscript -file EXAMPLES/test_cuba.wl
   ============================================================================ *)

(* --- Load package and CUBA helpers --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Get[FileNameJoin[{Directory[], "EXAMPLES", "cuba_common.wl"}]];

Print["tropical_eval.wl + cuba_common.wl loaded"];
Print[];

interfilesDir = FileNameJoin[{Directory[], "INTERFILES"}];
$cubaRows = {};

(* --------------------------------------------------------------------------
   runCubaCase: one standard (unlifted) test integral.
   opts (rules): "PrecisionGoal" (sector validation), "MCSamples"
   (0 = skip MC), "MCGate" (True -> 5-sigma gate, False -> report only),
   "SectorTol", "MaxEval"/"CuhreEpsRel" (passed to CUBA), "Exact".
   -------------------------------------------------------------------------- *)

Options[runCubaCase] = {
  "PrecisionGoal" -> 3, "MCSamples" -> 1000000, "MCGate" -> True,
  "SectorTol" -> 0.01, "MaxEval" -> 2000000, "CuhreEpsRel" -> 10^-6,
  "Exact" -> None};

runCubaCase[label_, tag_, spec_, OptionsPattern[]] := Module[
  {pg, ns, mcGate, sTol, exact, verts, fanData, vr, sectorSum,
   cuba, cuhre, mcRes, mcVal, mcSig,
   sectorDev, mcDevSigma, casePass = True, notes = {}},

  pg     = OptionValue["PrecisionGoal"];
  ns     = OptionValue["MCSamples"];
  mcGate = OptionValue["MCGate"];
  sTol   = OptionValue["SectorTol"];
  exact  = OptionValue["Exact"];

  Print["--- ", label, " ---"];

  (* tropical decomposition + per-sector NIntegrate sum *)
  verts   = PolytopeVertices[(Times @@ spec["Polynomials"])^(-1),
                             spec["Variables"]];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, pg];
  sectorSum = If[AssociationQ[vr], vr["SectorSum"], $Failed];

  (* CUBA Cuhre on the direct integrand *)
  cuba = runCubaCheck[spec, tag, "RunVegas" -> False,
    "MaxEval" -> OptionValue["MaxEval"],
    "CuhreEpsRel" -> OptionValue["CuhreEpsRel"]];
  cuhre = If[AssociationQ[cuba] && KeyExistsQ[cuba, "CUHRE"],
    cuba["CUHRE"], $Failed];

  If[cuhre === $Failed,
    Print["  CUBA unavailable or failed — cannot cross-check; FAIL"];
    AppendTo[$cubaRows, <|"Label" -> label, "Pass" -> False,
      "Notes" -> "CUBA failed"|>];
    Return[False]
  ];

  Print["  CUBA Cuhre:   ", cuhre["Value"], "  +/- ",
        ToString[cuhre["Error"], InputForm],
        "  (", cuhre["NEval"], " evals, fail=", cuhre["Fail"], ")"];

  If[exact =!= None,
    Print["  Exact:        ", N[exact],
          "   Cuhre |rel dev| = ",
          Abs[(cuhre["Value"] - exact)/exact]]];

  (* sector sum vs Cuhre *)
  If[sectorSum =!= $Failed,
    sectorDev = Abs[(sectorSum - cuhre["Value"])/cuhre["Value"]];
    Print["  Sector sum:   ", sectorSum, "   |rel dev| = ", sectorDev,
          "  (gate ", sTol, ")"];
    If[!TrueQ[sectorDev < sTol],
      casePass = False; AppendTo[notes, "sector-sum gate exceeded"];
    ],
    casePass = False; AppendTo[notes, "ValidateDecomposition failed"];
    sectorDev = Infinity;
  ];

  (* tropical C++ MC *)
  mcDevSigma = None;
  If[ns > 0,
    mcRes = EvaluateTropicalMC[spec, fanData, {{}},
      "NSamples" -> ns, "RunChecks" -> False, "Verbose" -> False,
      "WorkingDirectory" -> interfilesDir];
    If[AssociationQ[mcRes],
      Module[{r = mcRes["Results"][[1]]},
        mcVal = r["Re"] + I r["Im"];
        mcSig = Sqrt[r["ReErr"]^2 + r["ImErr"]^2];
        mcDevSigma = Abs[mcVal - cuhre["Value"]]/mcSig;
        Print["  Tropical MC:  ", mcVal, "  +/- (",
              r["ReErr"], ", ", r["ImErr"], ")  [", ns, " samples]"];
        Print["                |dev| from Cuhre = ",
              Abs[mcVal - cuhre["Value"]], "  = ", mcDevSigma, " sigma",
              If[!mcGate, "   [report only: heavy-tailed, no gate]", ""]];
        If[mcGate && !TrueQ[mcDevSigma < 5],
          casePass = False; AppendTo[notes, "MC > 5 sigma from Cuhre"];
        ];
      ],
      casePass = False; AppendTo[notes, "EvaluateTropicalMC failed"];
    ];
  ];

  Print["  ", If[casePass, "PASS", "FAIL"],
        If[notes =!= {}, "  (" <> StringRiffle[notes, "; "] <> ")", ""]];
  Print[];

  AppendTo[$cubaRows, <|"Label" -> label, "Pass" -> casePass,
    "SectorDev" -> sectorDev, "MCSigma" -> mcDevSigma,
    "Notes" -> StringRiffle[notes, "; "]|>];
  casePass
];

(* --------------------------------------------------------------------------
   Test 1: P = 1 + 2x1^2 + x2^2 + x1*x2^2 + 3x1^2*x2, B = -2 and -3
   -------------------------------------------------------------------------- *)

Module[{poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
  Do[
    runCubaCase["Test 1 (A=" <> ToString[A] <> ")", "t1a" <> ToString[A],
      <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
        "PolynomialExponents" -> {-A}, "Variables" -> {x[1], x[2]},
        "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>],
    {A, {2, 3}}
  ]
];

(* --------------------------------------------------------------------------
   Test 2: same polynomial, complex exponent B = -(2 + 0.5 I)
   -------------------------------------------------------------------------- *)

Module[{poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
  runCubaCase["Test 2 (A=2+0.5I)", "t2",
    <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
      "PolynomialExponents" -> {-(2 + 0.5 I)}, "Variables" -> {x[1], x[2]},
      "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>]
];

(* --------------------------------------------------------------------------
   Test 5: P = 1 + lam*x1^2 + x2^2 + x1*x2^2, B = -(2 + 0.5 I),
   at the five representative lam values the suite itself checks.
   -------------------------------------------------------------------------- *)

Module[{lamValues = {0.1, 2.5, 5.0, 7.5, 10.0}},
  Do[
    Module[{poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2},
      runCubaCase["Test 5 (lam=" <> ToString[lam] <> ")",
        "t5l" <> StringReplace[ToString[lam], "." -> "p"],
        <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
          "PolynomialExponents" -> {-(2 + 0.5 I)},
          "Variables" -> {x[1], x[2]},
          "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>]
    ],
    {lam, lamValues}
  ]
];

(* --------------------------------------------------------------------------
   Test 6: extreme coefficients.
   Case A: 10^6 coefficient (bounded integrand -> MC gated; also re-checked
           with the LIFTED MC, cf. Test 23B).
   Case B: constant term 10^-4 -> cleared polynomial has small constant
           term, plain uniform MC is heavy-tailed: MC reported, not gated.
   Case C: 10^8 coefficient, exponent -3 (bounded; PG2 in the suite).
   -------------------------------------------------------------------------- *)

runCubaCase["Test 6A (10^6 coeff)", "t6a",
  <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>];

(* Test 6A again through the LIFTED pipeline (= Test 23B's integral) *)
Module[{spec, resLifted, cuhreVal, dev, sig, pass},
  Print["--- Test 6A / 23B: lifted MC (k=2) vs Cuhre ---"];
  spec = <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  cuhreVal = Module[{c = runCubaCheck[spec, "t6a", "RunVegas" -> False]},
    If[AssociationQ[c] && KeyExistsQ[c, "CUHRE"], Re[c["CUHRE"]["Value"]],
      $Failed]];
  resLifted = Quiet@EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules" -> {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
    "NSamples" -> 1000000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> interfilesDir];
  pass = False;
  If[AssociationQ[resLifted] && NumericQ[cuhreVal],
    Module[{r = resLifted["Results"][[1]]},
      dev = r["Re"] - cuhreVal; sig = r["ReErr"];
      Print["  Lifted MC:    ", r["Re"], "  +/- ", r["ReErr"],
            "   dev from Cuhre = ", dev, "  = ", Abs[dev]/sig, " sigma"];
      pass = TrueQ[Abs[dev]/sig < 5];
    ],
    Print["  lifted run or CUBA failed"];
  ];
  Print["  ", If[pass, "PASS", "FAIL"]];
  Print[];
  AppendTo[$cubaRows, <|"Label" -> "Test 6A/23B lifted MC",
    "Pass" -> pass, "SectorDev" -> None,
    "MCSigma" -> If[NumericQ[dev] && NumericQ[sig], Abs[dev]/sig, None],
    "Notes" -> ""|>];
];

runCubaCase["Test 6B (coeffs 10^-4..10^4)", "t6b",
  <|"Polynomials" -> {10^-4 + 10^4 x[1]^2 + 10^-4 x[2]^2 +
                      10^4 x[1] x[2]^2 + x[1]^2 x[2]},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  "MCGate" -> False];

runCubaCase["Test 6C (10^8 coeff, B=-3)", "t6c",
  <|"Polynomials" -> {1 + 10^8 x[1]^3 x[2] + x[2]^3},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-3},
    "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  "PrecisionGoal" -> 2, "SectorTol" -> 0.05, "MaxEval" -> 10000000];

(* --------------------------------------------------------------------------
   Test 7: 3D and 4D
   -------------------------------------------------------------------------- *)

runCubaCase["Test 7A (3D)", "t7a",
  <|"Polynomials" -> {1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3]},
    "MonomialExponents" -> {0, 0, 0}, "PolynomialExponents" -> {-3},
    "Variables" -> {x[1], x[2], x[3]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>];

runCubaCase["Test 7B (4D)", "t7b",
  <|"Polynomials" -> {1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 +
                      x[1] x[2] + x[3] x[4]},
    "MonomialExponents" -> {0, 0, 0, 0}, "PolynomialExponents" -> {-4},
    "Variables" -> {x[1], x[2], x[3], x[4]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  "PrecisionGoal" -> 2, "SectorTol" -> 0.02, "MCSamples" -> 500000,
  "MaxEval" -> 5000000, "CuhreEpsRel" -> 10^-5];

(* --------------------------------------------------------------------------
   Test 23A: P = 1 + x1 + 10^-6 x2, B = -3.  Exact = 500000.
   The tropical side uses the LIFTED decomposition with the explicit fan
   (degenerate lifted polytope), exactly as in RunTest23A.
   -------------------------------------------------------------------------- *)

Module[{spec, exact, liftRules, lcRes, liftedSpec, liftData, explicitFan,
        vl, sectorSum, cuba, cuhre, sectorDev, cuhreDev, pass = True,
        notes = {}},
  Print["--- Test 23A (1 + x1 + 10^-6 x2)^-3, exact 500000 ---"];
  exact = 500000;
  spec = <|"Polynomials" -> {1 + x[1] + 10^-6 x[2]},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-3},
    "Variables" -> {x[1], x[2]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 1}, "k" -> 2|>};
  lcRes = LiftCoefficients[spec, liftRules];
  liftedSpec = lcRes["LiftedSpec"]; liftData = lcRes["LiftData"];
  explicitFan = {{{1,0,0},{0,1,0},{-1,-1,0},{0,2,-1},{0,-2,1}},
                 {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}}};
  vl = Quiet@ValidateLiftedDecomposition[
    spec, liftedSpec, explicitFan, liftData, {}, 3];
  sectorSum = If[AssociationQ[vl], vl["SectorSum"], $Failed];

  cuba = runCubaCheck[spec, "t23a", "RunVegas" -> False,
    "MaxEval" -> 10000000];
  cuhre = If[AssociationQ[cuba] && KeyExistsQ[cuba, "CUHRE"],
    cuba["CUHRE"], $Failed];

  If[cuhre =!= $Failed,
    cuhreDev = Abs[(cuhre["Value"] - exact)/exact];
    Print["  CUBA Cuhre:        ", cuhre["Value"], "  +/- ",
          ToString[cuhre["Error"], InputForm],
          "  (", cuhre["NEval"], " evals, fail=", cuhre["Fail"], ")"];
    Print["  Exact:             ", N[exact],
          "   Cuhre |rel dev| = ", cuhreDev];
    If[!TrueQ[cuhreDev < 0.01],
      pass = False; AppendTo[notes, "Cuhre vs exact gate exceeded"]],
    pass = False; AppendTo[notes, "CUBA failed"];
  ];

  If[sectorSum =!= $Failed,
    sectorDev = Abs[(sectorSum - exact)/exact];
    Print["  Lifted sector sum: ", sectorSum,
          "   |rel dev| from exact = ", sectorDev, "  (gate 0.01)"];
    If[!TrueQ[sectorDev < 0.01],
      pass = False; AppendTo[notes, "lifted sector sum gate exceeded"]],
    pass = False; AppendTo[notes, "ValidateLiftedDecomposition failed"];
  ];

  Print["  ", If[pass, "PASS", "FAIL"],
        If[notes =!= {}, "  (" <> StringRiffle[notes, "; "] <> ")", ""]];
  Print[];
  AppendTo[$cubaRows, <|"Label" -> "Test 23A (lifted, exact 5e5)",
    "Pass" -> pass, "SectorDev" -> sectorDev, "MCSigma" -> None,
    "Notes" -> StringRiffle[notes, "; "]|>];
];

(* --------------------------------------------------------------------------
   Test 23D: P = 1 + 10^6 x1 (1D), B = -2.  Exact = 10^-6.
   Tropical side: lifted decomposition with the explicit k=1 fan, as in
   RunTest23D.  CUBA side: 1D integrand padded to ndim=2 for Cuhre.
   NOTE: the entire integral lives in a sliver x < ~10^-6 (t < 10^-6
   after compactification).  This row primarily gates the SECTOR SUM
   against the exact value; the Cuhre result is reported and gated only
   loosely (5%) — finding an invisible 10^-6 boundary sliver is exactly
   the kind of problem direct quadrature can fail at silently, which is
   the reason the tropical decomposition exists.
   -------------------------------------------------------------------------- *)

Module[{spec, exact, liftRules, lcRes, liftedSpec, liftData, raysK1, sectsK1,
        explicitFan, vl, sectorSum, cuba, cuhre, sectorDev, cuhreDev,
        pass = True, notes = {}},
  Print["--- Test 23D (1 + 10^6 x)^-2 (1D), exact 10^-6 ---"];
  exact = 10^-6;
  spec = <|"Polynomials" -> {1 + 10^6 x[1]},
    "MonomialExponents" -> {0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1]},
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 1|>};
  lcRes = LiftCoefficients[spec, liftRules];
  liftedSpec = lcRes["LiftedSpec"]; liftData = lcRes["LiftData"];
  raysK1  = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
  sectsK1 = Table[{i, If[i < Length[raysK1], i+1, 1]}, {i, Length[raysK1]}];
  explicitFan = {raysK1, sectsK1};
  vl = Quiet@ValidateLiftedDecomposition[
    spec, liftedSpec, explicitFan, liftData, {}, 4];
  sectorSum = If[AssociationQ[vl], vl["SectorSum"], $Failed];

  cuba = runCubaCheck[spec, "t23d", "RunVegas" -> False,
    "MaxEval" -> 10000000];
  cuhre = If[AssociationQ[cuba] && KeyExistsQ[cuba, "CUHRE"],
    cuba["CUHRE"], $Failed];

  If[sectorSum =!= $Failed,
    sectorDev = Abs[(sectorSum - exact)/exact];
    Print["  Lifted sector sum: ", sectorSum,
          "   |rel dev| from exact = ", sectorDev, "  (gate 0.005)"];
    If[!TrueQ[sectorDev < 0.005],
      pass = False; AppendTo[notes, "lifted sector sum gate exceeded"]],
    pass = False; AppendTo[notes, "ValidateLiftedDecomposition failed"];
  ];

  If[cuhre =!= $Failed,
    cuhreDev = Abs[(cuhre["Value"] - exact)/exact];
    Print["  CUBA Cuhre (1D padded): ", cuhre["Value"], "  +/- ",
          ToString[cuhre["Error"], InputForm],
          "  (", cuhre["NEval"], " evals, fail=", cuhre["Fail"], ")"];
    Print["    |rel dev| from exact = ", cuhreDev, "  (gate 0.05)"];
    If[!TrueQ[cuhreDev < 0.05],
      pass = False;
      AppendTo[notes, "Cuhre missed the 10^-6 sliver (see header note)"]],
    pass = False; AppendTo[notes, "CUBA failed"];
  ];

  Print["  ", If[pass, "PASS", "FAIL"],
        If[notes =!= {}, "  (" <> StringRiffle[notes, "; "] <> ")", ""]];
  Print[];
  AppendTo[$cubaRows, <|"Label" -> "Test 23D (lifted 1D, exact 1e-6)",
    "Pass" -> pass, "SectorDev" -> sectorDev, "MCSigma" -> None,
    "Notes" -> StringRiffle[notes, "; "]|>];
];

(* --------------------------------------------------------------------------
   Summary
   -------------------------------------------------------------------------- *)

Module[{nPass, nFail},
  nPass = Count[$cubaRows, r_ /; TrueQ[r["Pass"]]];
  nFail = Length[$cubaRows] - nPass;
  Print["================================================================"];
  Print["  CUBA cross-check summary (", Length[$cubaRows], " cases)"];
  Print["================================================================"];
  Module[{fmt},
    (* single-line scientific notation (avoids OutputForm superscripts) *)
    fmt[v_] := ToString[CForm[SetPrecision[v, 3]]];
    Do[
      Print["  ", StringPadRight[row["Label"], 34],
        If[TrueQ[row["Pass"]], "PASS", "FAIL"],
        If[NumericQ[row["SectorDev"]],
          "   sector dev " <> fmt[row["SectorDev"]], ""],
        If[NumericQ[row["MCSigma"]],
          "   MC " <> fmt[row["MCSigma"]] <> " sigma", ""],
        If[row["Notes"] =!= "", "   [" <> row["Notes"] <> "]", ""]],
      {row, $cubaRows}
    ];
  ];
  Print["----------------------------------------------------------------"];
  Print["  ", nPass, " PASSED, ", nFail, " FAILED"];
  Print["================================================================"];
];
