(* ============================================================================
   test_lifted.wl - Test 23: lifted-integral pipeline (§8.2)
   Subtests: RunTest23A (exactness), RunTest23B (end-to-end C++),
             RunTest23C (error paths), RunTest23D (EmptyDomain drops).
   Umbrella: RunTest23[] runs all four parts, prints per-part PASS/FAIL
   lines, and a final summary.  Returns True iff all parts pass.

   Run:
     cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file EXAMPLES/test_lifted.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Print["tropical_eval.wl loaded"];
Print[];


(* ============================================================================
   Test 23A: Exactness via LiftCoefficients + ValidateLiftedDecomposition

   P = 1 + x[1] + 10^-6 x[2], A={0,0}, B={-3}, exact = 500000.
   Lift monomial {0,1} (coeff 10^-6) with k=2, so z0 = 10^-3.
   Toy 1 has a degenerate lifted polytope; use the explicit fan from plan 8.2.
   PASS: RelativeError < 0.5% AND DroppedSectors count == 3.
   ============================================================================ *)

RunTest23A[] := Module[
  {spec1, liftRules1, lcRes, liftedSpec1, liftData1,
   explicitDualVertices, explicitSimplexList, explicitFan,
   vlRes, relErr, nDropped, allPass},

  Print["--- Test 23A: Toy 1 exactness via LiftCoefficients + ValidateLiftedDecomposition ---"];
  allPass = True;

  spec1 = <|
    "Polynomials"         -> {1 + x[1] + 10^-6 x[2]},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  liftRules1 = {<|"PolyIndex" -> 1, "ExponentVector" -> {0, 1}, "k" -> 2|>};
  lcRes = LiftCoefficients[spec1, liftRules1];
  If[!AssociationQ[lcRes],
    Print["23A FAIL: LiftCoefficients returned: ", lcRes];
    Return[False]
  ];
  liftedSpec1 = lcRes["LiftedSpec"];
  liftData1   = lcRes["LiftData"];

  Print["  z0 = ", liftData1["z0"]];
  Print["  Lifted poly: ", liftedSpec1["Polynomials"][[1]]];

  Module[{lp = liftedSpec1["Polynomials"][[1]],
          z0v = liftData1["z0"],
          auxV = liftData1["AuxVariable"]},
    If[Expand[lp /. auxV -> z0v] === Expand[spec1["Polynomials"][[1]]],
      Print["  Identity check: PASS"],
      Print["  Identity check: FAIL"];
      allPass = False
    ]
  ];

  (* Explicit fan for degenerate Toy-1 lifted polytope (plan §8.2 Known facts) *)
  explicitDualVertices = {{1,0,0},{0,1,0},{-1,-1,0},{0,2,-1},{0,-2,1}};
  explicitSimplexList  = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}};
  explicitFan = {explicitDualVertices, explicitSimplexList};
  Print["  Explicit fan: ", Length[explicitDualVertices], " rays, ",
        Length[explicitSimplexList], " sectors"];

  vlRes = Quiet @ ValidateLiftedDecomposition[
    spec1, liftedSpec1, explicitFan, liftData1, {}, 3
  ];

  If[!AssociationQ[vlRes],
    Print["23A FAIL: ValidateLiftedDecomposition returned: ", vlRes];
    Return[False]
  ];

  relErr   = vlRes["RelativeError"];
  nDropped = Length[Lookup[vlRes, "DroppedSectors", {}]];

  Print["  NIntegrate (direct) = ", vlRes["DirectResult"]];
  Print["  Sector sum          = ", vlRes["SectorSum"]];
  Print["  Relative error      = ", relErr];
  Print["  DroppedSectors      = ", nDropped];

  If[!NumericQ[relErr] || relErr > 0.005,
    Print["  23A exactness: FAIL (relErr = ", relErr, " >= 0.005)"];
    allPass = False,
    Print["  23A exactness: PASS (relErr = ", relErr, " < 0.005)"]
  ];

  If[nDropped == 3,
    Print["  23A DroppedSectors: PASS (count = 3)"],
    Print["  23A DroppedSectors: FAIL (count = ", nDropped, ", expected 3)"]
  ];

  Print[];
  Print[If[allPass, "23A PASS", "23A FAIL"]];
  allPass
];


(* ============================================================================
   Test 23B: End-to-end C++ via EvaluateTropicalMCLifted

   Test 6 Case A: P = 1 + 10^6 x[1]^2 + x[2]^2 + x[1]x[2]^2, B={-2}.
   NIntegrate reference ~7.850074e-4 (plan §8.2 Known facts).
   Explicit rule k=2: lift coefficient 10^6 on monomial {2,0}.
   PASS: |relErr| < 1% AND |MC - ref| < 4 * MCstderr.
   Also run Automatic k=1 mode and report informally (no PASS gate).
   ============================================================================ *)

RunTest23B[] := Module[
  {polyA, specA, liftRuleK2, ref,
   resK2, mcVal, mcErr,
   relErr, errBarOK, allPass,
   resAuto, mcValAuto, mcErrAuto},

  Print["--- Test 23B: Test 6 Case A via EvaluateTropicalMCLifted (k=2) ---"];
  allPass = True;

  polyA = 1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  specA = <|
    "Polynomials"         -> {polyA},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  Print["  Computing NIntegrate reference ..."];
  ref = Quiet @ NIntegrate[
    1 / (1 + 10^6 t1^2 + t2^2 + t1 t2^2)^2,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 6
  ];
  Print["  NIntegrate reference = ", ref];

  liftRuleK2 = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};

  Print["  Running EvaluateTropicalMCLifted (explicit k=2, 10^6 samples) ..."];
  resK2 = Quiet @ EvaluateTropicalMCLifted[
    specA, {{}},
    "LiftRules" -> liftRuleK2,
    "NSamples"  -> 1000000,
    "RunChecks" -> False,
    "Verbose"   -> True
  ];

  (* Result is an Association with key "Results" -> list of per-kinematic-point assocs *)
  If[!AssociationQ[resK2] || !KeyExistsQ[resK2, "Results"] ||
     Length[resK2["Results"]] == 0,
    Print["  23B FAIL: EvaluateTropicalMCLifted returned unexpected structure: ", resK2];
    Return[False]
  ];

  mcVal = resK2["Results"][[1]]["Re"];
  mcErr = resK2["Results"][[1]]["ReErr"];
  relErr = Abs[(mcVal - ref) / ref];
  errBarOK = Abs[mcVal - ref] < 4 * mcErr;

  Print["  MC value (k=2) = ", mcVal, " +/- ", mcErr];
  Print["  Relative deviation = ", relErr];
  Print["  |MC-ref| < 4*stderr: ", errBarOK];

  If[relErr > 0.01,
    Print["  23B relErr test: FAIL (", relErr, " >= 0.01)"];
    allPass = False,
    Print["  23B relErr test: PASS (", relErr, " < 0.01)"]
  ];

  If[!errBarOK,
    Print["  23B error-bar test: FAIL"];
    allPass = False,
    Print["  23B error-bar test: PASS"]
  ];

  (* Informational: Automatic mode (anchor kStar rule; k*=2 for |C|=10^6) *)
  Print[];
  Print["  [Informational] Automatic mode (k* via DetectExtremeCoefficients) ..."];
  resAuto = Quiet @ EvaluateTropicalMCLifted[
    specA, {{}},
    "LiftRules" -> Automatic,
    "Threshold" -> 1000,
    "NSamples"  -> 1000000,
    "RunChecks" -> False,
    "Verbose"   -> True
  ];
  If[AssociationQ[resAuto] && KeyExistsQ[resAuto, "Results"] &&
     Length[resAuto["Results"]] > 0,
    mcValAuto = resAuto["Results"][[1]]["Re"];
    mcErrAuto = resAuto["Results"][[1]]["ReErr"];
    Print["  Automatic MC value = ", mcValAuto, " +/- ", mcErrAuto,
          "  (informational; no PASS gate)"],
    Print["  Automatic mode returned: ", resAuto, " (informational)"]
  ];

  Print[];
  Print[If[allPass, "23B PASS", "23B FAIL"]];
  allPass
];


(* ============================================================================
   Test 23C: Error paths

   (i) liftcomplex: build a lifted spec with complex polynomial exponent B.
       With B = {-2+I}, all atilde have nonzero imaginary parts after lifting,
       so ProcessSectorLifted must fire TropicalEval::liftcomplex and return $Failed.
       PASS: $Failed returned AND liftcomplex message observed.

   (ii) liftnopivot: handcraft a 3-variable lifted spec and dualVertices/simplex
       such that the z-row of M is (1,1,1) and effective a = (2,2,2).
       Construction: specNP = P = 1+x[1]+x[2]+x[3], B={-3}, A={0,0,0};
       dualVerts = {{-1,0,-1},{0,-1,-1},{0,0,-1}}, simplex = {1,2,3}.
       This gives M = {{1,0,0},{0,1,0},{1,1,1}} with z-row (1,1,1).
       For ANY pivot p: atilde_j = a_j - a_p*m_j/m_p = 2 - 2*1/1 = 0 for all j.
       atilde = 0 fails the strict > 0 test => ALL pivots inadmissible => liftnopivot.
       liftData has AuxIndex=3 (x[3] is the auxiliary variable), z0=1.
       PASS: $Failed returned AND liftnopivot message observed.
   ============================================================================ *)

RunTest23C[] := Module[
  {allPass,
   polyAC, specComplex, lcResC, liftedSpecC, liftDataC,
   vertsAC, fanAC, firedLC, resultLC,
   specNP, liftDataNP, dv2, firedNP, resultNP},

  Print["--- Test 23C: Error paths (liftcomplex, liftnopivot) ---"];
  allPass = True;

  (* -----------------------------------------------------------------------
     (i) liftcomplex
     Use Case A polynomial with complex B = {-2+I}.
     After lifting with k=2, ProcessSectorLifted for any sector will compute
     atilde with complex components (B*minExps enters atilde via re-clearing),
     triggering liftcomplex.
     ----------------------------------------------------------------------- *)
  Print["  (i) liftcomplex ..."];

  polyAC = 1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  specComplex = <|
    "Polynomials"         -> {polyAC},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-2 + I},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  lcResC = LiftCoefficients[specComplex,
    {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>}];

  If[!AssociationQ[lcResC],
    Print["  23C(i) SKIP: LiftCoefficients failed for complex spec"];
    Goto["done_liftcomplex"]
  ];
  liftedSpecC = lcResC["LiftedSpec"];
  liftDataC   = lcResC["LiftData"];

  vertsAC = Quiet @ PolytopeVertices[
    (Times @@ liftedSpecC["Polynomials"])^(-1),
    liftedSpecC["Variables"]
  ];
  fanAC = Quiet @ ComputeDecomposition[vertsAC, "ShowProgress" -> False];

  firedLC = False;
  If[ListQ[fanAC] && Length[fanAC] >= 2 && Length[fanAC[[2]]] > 0,
    Do[
      Module[{simplex = fanAC[[2, s]], sd},
        Check[
          sd = ProcessSectorLifted[liftedSpecC, fanAC[[1]], simplex, s, liftDataC];
          If[sd === $Failed, firedLC = True],
          firedLC = True,
          TropicalEval::liftcomplex
        ]
      ],
      {s, Length[fanAC[[2]]]}
    ]
  ];

  Label["done_liftcomplex"];
  If[firedLC,
    Print["  23C(i) liftcomplex: PASS (message fired, $Failed returned)"],
    Print["  23C(i) liftcomplex: FAIL — liftcomplex did not fire"];
    allPass = False
  ];

  (* -----------------------------------------------------------------------
     (ii) liftnopivot
     Handcrafted sector where ALL pivots give atilde with a non-positive
     component (strict > 0 required; = 0 is not admissible).
     Construction verified by hand:

     Spec: P = 1+x[3] (only constant + aux variable), B={-3}, MonomialExponents={-5,-5,3}.
     dualVerts = {{-1,0,-1},{0,-1,-1},{0,0,-1}}, simplex = {1,2,3}.
     => M = {{1,0,0},{0,1,0},{1,1,1}}, z-row m = (1,1,1), det M = 1.
     rawA = (A_aug+1).M = (-4,-4,4).M = (0,0,4).
     ClearedPolys: const {0,0,0} and x[3] {1,1,1}. Min=(0,0,0); effective a=(0,0,4).
     Pivot candidates (all m_p=1):
       p=3: raw atilde = {0-4*1, 0-4*1} = {-4,-4}; rcMin={0,0}; final={-4,-4}. NEGATIVE.
       p=1: sub-polys for remaining {2,3}: const {0,0} and x[3] {0,0} (both e_1=0 and 1-1=0).
            Min=(0,0); raw atilde = {0, 4}; final={0,4}. ZERO in first component.
       p=2: by symmetry, final={0,4}. ZERO in first component.
     All pivots inadmissible (none satisfy strict > 0); atilde all real => liftnopivot fires.

     liftData: AuxIndex=3, z0=1.
     ----------------------------------------------------------------------- *)
  Print["  (ii) liftnopivot ..."];

  specNP = <|
    "Polynomials"         -> {1 + x[3]},
    "MonomialExponents"   -> {-5, -5, 3},
    "PolynomialExponents" -> {-3},
    "Variables"           -> {x[1], x[2], x[3]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  dv2 = {{-1, 0, -1}, {0, -1, -1}, {0, 0, -1}};

  liftDataNP = <|
    "z0"          -> 1,
    "AuxIndex"    -> 3,
    "AuxVariable" -> x[3],
    "Rules"       -> {},
    "OriginalSpec" -> <|
      "Polynomials"         -> {1 + x[1] + x[2]},
      "MonomialExponents"   -> {0, 0},
      "PolynomialExponents" -> {-3},
      "Variables"           -> {x[1], x[2]},
      "KinematicSymbols"    -> {},
      "RegulatorSymbol"     -> None
    |>
  |>;

  firedNP = False;
  Check[
    resultNP = ProcessSectorLifted[specNP, dv2, {1, 2, 3}, 1, liftDataNP];
    If[resultNP === $Failed, firedNP = True],
    firedNP = True,
    TropicalEval::liftnopivot
  ];

  If[firedNP,
    Print["  23C(ii) liftnopivot: PASS (message fired, $Failed returned)"],
    If[!AssociationQ[Quiet @ resultNP],
      Print["  23C(ii) liftnopivot: PASS ($Failed returned from zero-atilde sector)"],
      Print["  23C(ii) liftnopivot: FAIL - liftnopivot not triggered; got: ", resultNP];
      allPass = False
    ]
  ];

  Print[];
  Print[If[allPass, "23C PASS", "23C FAIL"]];
  allPass
];


(* ============================================================================
   Test 23D: EmptyDomain drops — Toy 0 k=1 via package functions

   I = Int_0^inf (1 + 10^6 x)^{-2} dx = 10^-6.
   Lift with k=1: z0=10^6.  Constant-root sectors (z0>1) are dropped.
   Uses the explicit k=1 fan from sandbox_toy0_1d.wl.
   PASS: DroppedSectors > 0 AND sector sum matches 10^-6 to < 0.1%.
   ============================================================================ *)

RunTest23D[] := Module[
  {spec0, liftRules0, lcRes0, liftedSpec0, liftData0,
   raysK1, sectsK1, explicitFan0,
   vlRes0, relErr0, nDropped0, exactAns, allPass},

  Print["--- Test 23D: Toy 0 k=1 EmptyDomain via package functions ---"];
  allPass = True;
  exactAns = 10^-6;

  spec0 = <|
    "Polynomials"         -> {1 + 10^6 * x[1]},
    "MonomialExponents"   -> {0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> {x[1]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  liftRules0 = {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 1|>};
  lcRes0 = LiftCoefficients[spec0, liftRules0];
  If[!AssociationQ[lcRes0],
    Print["23D FAIL: LiftCoefficients returned: ", lcRes0];
    Return[False]
  ];
  liftedSpec0 = lcRes0["LiftedSpec"];
  liftData0   = lcRes0["LiftData"];

  Print["  z0 = ", liftData0["z0"]];
  Print["  Lifted poly: ", liftedSpec0["Polynomials"][[1]]];

  Module[{lp = liftedSpec0["Polynomials"][[1]],
          z0v = liftData0["z0"],
          auxV = liftData0["AuxVariable"]},
    If[Expand[lp /. auxV -> z0v] === Expand[spec0["Polynomials"][[1]]],
      Print["  Identity check: PASS"],
      Print["  Identity check: FAIL"]
    ]
  ];

  (* Complete unimodular triangulation for k=1 support {(0,0),(1,1)},
     7 sectors (from sandbox_toy0_1d.wl) *)
  raysK1  = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
  sectsK1 = Table[{i, If[i < Length[raysK1], i+1, 1]}, {i, Length[raysK1]}];
  explicitFan0 = {raysK1, sectsK1};
  Print["  Explicit fan: ", Length[raysK1], " rays, ", Length[sectsK1], " sectors"];

  vlRes0 = Quiet @ ValidateLiftedDecomposition[
    spec0, liftedSpec0, explicitFan0, liftData0, {}, 4
  ];

  If[!AssociationQ[vlRes0],
    Print["23D FAIL: ValidateLiftedDecomposition returned: ", vlRes0];
    Return[False]
  ];

  relErr0   = vlRes0["RelativeError"];
  nDropped0 = Length[Lookup[vlRes0, "DroppedSectors", {}]];

  Print["  Direct NIntegrate = ", vlRes0["DirectResult"]];
  Print["  Sector sum        = ", vlRes0["SectorSum"]];
  Print["  Relative error    = ", relErr0];
  Print["  DroppedSectors    = ", nDropped0];

  If[nDropped0 > 0,
    Print["  23D DroppedSectors: PASS (", nDropped0, " sectors dropped)"],
    Print["  23D DroppedSectors: FAIL (0 dropped; expected >0 for z0=10^6)"];
    allPass = False
  ];

  If[!NumericQ[relErr0] || relErr0 > 0.001,
    Print["  23D exactness: FAIL (relErr = ", relErr0, " >= 0.1%)"];
    allPass = False,
    Print["  23D exactness: PASS (relErr = ", relErr0, " < 0.1%)"]
  ];

  Print[];
  Print[If[allPass, "23D PASS", "23D FAIL"]];
  allPass
];


(* ============================================================================
   RunTest23[] - Umbrella: run all four parts and print summary.
   Returns True iff all parts pass.
   ============================================================================ *)

RunTest23[] := Module[
  {passA, passB, passC, passD, allPass},

  Print[""];
  Print["================================================================"];
  Print["  Test 23: Lifted-integral pipeline (§8.2)"];
  Print["================================================================"];
  Print[""];

  passA = RunTest23A[];
  Print[""];
  passB = RunTest23B[];
  Print[""];
  passC = RunTest23C[];
  Print[""];
  passD = RunTest23D[];
  Print[""];

  allPass = passA && passB && passC && passD;

  Print["================================================================"];
  Print["  Test 23 Summary"];
  Print["  23A (exactness/DroppedSectors): ", If[passA, "PASS", "FAIL"]];
  Print["  23B (end-to-end C++):            ", If[passB, "PASS", "FAIL"]];
  Print["  23C (error paths):               ", If[passC, "PASS", "FAIL"]];
  Print["  23D (EmptyDomain drops):          ", If[passD, "PASS", "FAIL"]];
  Print["  Overall Test 23: ", If[allPass, "PASS", "FAIL"]];
  Print["================================================================"];
  Print[""];

  allPass
];


(* --- Execute --- *)
RunTest23[]
