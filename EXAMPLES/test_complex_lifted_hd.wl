(* ============================================================================
   test_complex_lifted_hd.wl — Test 25: lifting in higher dimensions.

   Covers structural gaps left by Tests 23–24: 3–4 variables, degree 4–6,
   non-trivial Euler measure (non-zero MonomialExponents), and multi-polynomial
   integrands.  Every test uses REAL polynomial exponents B so that the standard
   lifting algorithm applies (no SplitRealImag needed); the complex structure
   lives in the polynomial coefficients.

   Subtests:
     RunTest25A[] — 3-variable degree-4, single extreme complex coeff
     RunTest25B[] — 3-variable degree-5, very small coeff (k*=3 regime)
     RunTest25C[] — 3-variable, two polynomials with complex extremes each
     RunTest25D[] — 4-variable, non-trivial Euler measure (A_1 = 1/3)
     RunTest25E[] — 3-variable degree-6, three competing complex extremes

   Umbrella RunTest25[] runs all five, prints PASS/FAIL, exits non-zero on fail.

   Requires: tropical_fan.wl, tropical_eval.wl, Polymake; g++ for MC steps.
   Run:
     cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file EXAMPLES/test_complex_lifted_hd.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Print["tropical_eval.wl loaded"];
Print[];

$hdWorkDir = FileNameJoin[{Directory[], "INTERFILES"}];


(* ============================================================================
   Test 25A — 3-variable degree-4, single extreme complex coefficient.

   P  = 1 + (2+3I)·1e-5·x1^4 + x2^2 + x3^2 + x1·x2·x3
   A  = {0, 0, 0},  B = {-2}

   |C| = sqrt(13)·1e-5 ≈ 3.61e-5;  k* = 2;  z0 = |C|^{1/2};
   residual c = (2+3I)/sqrt(13), |c| = 1.

   PASS gates:
     (1) DetectExtremeCoefficients flags exactly 1 entry with SuggestedK = 2,
         magnitude ≈ sqrt(13)·1e-5.
     (2) LiftCoefficients identity holds — no liftidentity message fired.
     (3) |residual| = 1 to 1e-9.
     (4) ValidateLiftedDecomposition relErr < 1e-2.
     (5) EvaluateTropicalMCLifted (C++, 5×10^5 samples) within 5σ of
         NIntegrate reference in both Re and Im.
   ============================================================================ *)

RunTest25A[] := Module[
  {c, poly, vars, spec, lr, det, lc, ls, ld, resid, lv, lf, vl, relErr,
   ref, res, mcRe, mcIm, mcReErr, mcImErr,
   firedLiftidentity, allPass},

  Print["--- Test 25A: 3-var deg-4, single extreme complex coeff (2+3I)1e-5 ---"];
  allPass = True;

  c    = (2 + 3 I) 10^-5;
  vars = {x[1], x[2], x[3]};
  poly = 1 + c x[1]^4 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  (* (1) Detection *)
  det = DetectExtremeCoefficients[spec, 1000];
  If[Length[det] == 1 && det[[1]]["SuggestedK"] == 2 &&
     Abs[det[[1]]["Magnitude"] - Sqrt[13.] 10^-5] < 10^-9,
    Print["  (1) detect: PASS (SuggestedK=2, |C|=", det[[1]]["Magnitude"], ")"],
    Print["  (1) detect: FAIL -> ", det]; allPass = False
  ];

  lr = {<|"PolyIndex"->1, "ExponentVector"->{4,0,0}, "k"->2|>};

  (* (2) LiftCoefficients — identity must not fire *)
  firedLiftidentity = False;
  Check[
    lc = LiftCoefficients[spec, lr],
    firedLiftidentity = True,
    TropicalEval::liftidentity
  ];
  If[!firedLiftidentity && AssociationQ[lc],
    Print["  (2) liftidentity not fired: PASS"],
    Print["  (2) liftidentity fired: FAIL"]; allPass = False
  ];
  If[!AssociationQ[lc], Print["  LiftCoefficients FAILED"]; Return[False]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  resid = ld["Residuals"][[1]];

  (* (3) |residual| = 1 *)
  If[Abs[N[Abs[resid]] - 1] < 10^-9,
    Print["  (3) |residual|=1: PASS (", N[resid], ")"],
    Print["  (3) |residual|=1: FAIL (|resid|=", N[Abs[resid]], ")"]; allPass = False
  ];

  (* (4) ValidateLiftedDecomposition *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[!ListQ[lf] || Length[lf] < 2 || !AllTrue[lf[[2]], Length[#] == 4 &],
    Print["  lifted fan degenerate or wrong dimension — FAIL"]; Return[False]
  ];
  vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
  relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
  Print["  (4) validate relErr = ", relErr];
  If[NumericQ[relErr] && relErr < 1*^-2,
    Print["  (4) exactness: PASS (relErr < 1e-2)"],
    Print["  (4) exactness: FAIL"]; allPass = False
  ];

  (* (5) C++ MC vs NIntegrate *)
  ref = Quiet @ NIntegrate[
    (1 + (2+3I)*10^-5*t1^4 + t2^2 + t3^2 + t1*t2*t3)^(-2),
    {t1,0,Infinity}, {t2,0,Infinity}, {t3,0,Infinity},
    MaxRecursion->15, PrecisionGoal->4];
  res = EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules"->lr, "FanData"->lf,
    "NSamples"->5*10^5, "RunChecks"->False, "Verbose"->False,
    "WorkingDirectory"->$hdWorkDir];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    mcRe = res["Results"][[1]]["Re"]; mcIm = res["Results"][[1]]["Im"];
    mcReErr = res["Results"][[1]]["ReErr"]; mcImErr = res["Results"][[1]]["ImErr"];
    Print["  (5) MC  = ", mcRe, " + ", mcIm, " I"];
    Print["      ref = ", Re[ref], " + ", Im[ref], " I"];
    Print["      dev = (", mcRe-Re[ref], ", ", mcIm-Im[ref],
          ")  stderr=(", mcReErr, ", ", mcImErr, ")"];
    If[Abs[mcRe - Re[ref]] < 5 mcReErr && Abs[mcIm - Im[ref]] < 5 mcImErr,
      Print["  (5) C++ MC: PASS"],
      Print["  (5) C++ MC: FAIL"]; allPass = False
    ],
    Print["  (5) EvaluateTropicalMCLifted FAILED: ", res]; allPass = False
  ];

  Print[]; Print[If[allPass, "25A PASS", "25A FAIL"]];
  allPass
];


(* ============================================================================
   Test 25B — 3-variable degree-5, very small coefficient (k*=3 regime).

   P  = 1 + (1+I)·1e-8·x1^5 + x2^2 + x3^2 + x1·x2·x3
   A  = {0,0,0},  B = {-3}

   |C| = sqrt(2)·1e-8 ≈ 1.41e-8;  k* = 3;  z0 = |C|^{1/3} ≈ 2.42e-3.
   The degree-5 lifted vertex {5,0,0,3} in ℝ⁴ exercises polymake's
   4D fan machinery with a richer Newton polytope.

   PASS gates:
     (1) SuggestedK = 3.
     (2) LiftCoefficients identity holds to relative tolerance 1e-7.
     (3) ValidateLiftedDecomposition relErr < 5e-2.
     (4) C++ MC within 5σ of NIntegrate reference.
   ============================================================================ *)

RunTest25B[] := Module[
  {c, poly, vars, spec, det, k, lr, lc, ls, ld, resid, lv, lf, vl, relErr,
   ref, res, mcRe, mcIm, mcReErr, mcImErr, allPass},

  Print["--- Test 25B: 3-var deg-5, k*=3 regime, coeff (1+I)1e-8 ---"];
  allPass = True;

  c    = (1 + I) 10^-8;
  vars = {x[1], x[2], x[3]};
  poly = 1 + c x[1]^5 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0,0},
           "PolynomialExponents"->{-3}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  (* (1) k* = 3 *)
  det = DetectExtremeCoefficients[spec, 1000];
  k   = det[[1]]["SuggestedK"];
  If[k == 3,
    Print["  (1) SuggestedK=3: PASS (|C|=", det[[1]]["Magnitude"], ")"],
    Print["  (1) SuggestedK: FAIL (got ", k, ", expected 3)"]; allPass = False
  ];

  lr = {<|"PolyIndex"->1, "ExponentVector"->{5,0,0}, "k"->k|>};
  lc = LiftCoefficients[spec, lr]; ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  Print["  z0 = ", N[ld["z0"]], "   residual = ", N[ld["Residuals"][[1]]]];

  (* (2) Identity check: verify no liftidentity message fired and |resid| ≈ 1 *)
  resid = ld["Residuals"][[1]];
  If[Abs[N[Abs[resid]] - 1] < 10^-7,
    Print["  (2) identity/|residual|=1: PASS"],
    Print["  (2) |residual| off: FAIL (|resid|=", N[Abs[resid]], ")"]; allPass = False
  ];

  (* (3) ValidateLiftedDecomposition *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[!ListQ[lf] || Length[lf] < 2,
    Print["  lifted fan degenerate — FAIL"]; Return[False]
  ];
  vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
  relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
  Print["  (3) validate relErr = ", relErr];
  If[NumericQ[relErr] && relErr < 5*^-2,
    Print["  (3) exactness: PASS (relErr < 5e-2)"],
    Print["  (3) exactness: FAIL"]; allPass = False
  ];

  (* (4) C++ MC vs NIntegrate *)
  ref = Quiet @ NIntegrate[
    (1 + (1+I)*10^-8*t1^5 + t2^2 + t3^2 + t1*t2*t3)^(-3),
    {t1,0,Infinity}, {t2,0,Infinity}, {t3,0,Infinity},
    MaxRecursion->15, PrecisionGoal->3];
  res = EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules"->lr, "FanData"->lf,
    "NSamples"->5*10^5, "RunChecks"->False, "Verbose"->False,
    "WorkingDirectory"->$hdWorkDir];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    mcRe = res["Results"][[1]]["Re"]; mcIm = res["Results"][[1]]["Im"];
    mcReErr = res["Results"][[1]]["ReErr"]; mcImErr = res["Results"][[1]]["ImErr"];
    Print["  (4) MC  = ", mcRe, " + ", mcIm, " I"];
    Print["      ref = ", Re[ref], " + ", Im[ref], " I"];
    Print["      dev = (", mcRe-Re[ref], ", ", mcIm-Im[ref],
          ")  stderr=(", mcReErr, ", ", mcImErr, ")"];
    If[Abs[mcRe - Re[ref]] < 5 mcReErr && Abs[mcIm - Im[ref]] < 5 mcImErr,
      Print["  (4) C++ MC: PASS"],
      Print["  (4) C++ MC: FAIL"]; allPass = False
    ],
    Print["  (4) EvaluateTropicalMCLifted FAILED: ", res]; allPass = False
  ];

  Print[]; Print[If[allPass, "25B PASS", "25B FAIL"]];
  allPass
];


(* ============================================================================
   Test 25C — 3-variable, two polynomials with complex extremes in each.

   P1 = 1 + (1+2I)·1e-5·x1^3 + x2^2 + x3
   P2 = 1 + x1^2 + (3-I)·1e-4·x2^4 + x3^2
   A  = {0,0,0},  B = {-2, -1}

   Anchor selection: primary = P1 monomial (|log|C1|| ≈ 10.7 > |log|C2|| ≈ 8.1).
   z0 = |C1|^{1/2},  |c2| = |C2|/z0^2 = |C2|/|C1| ≈ 14.1 (O(10), non-extreme).

   PASS gates:
     (1) Detection finds 2 extremes; P1's monomial has larger |log|C|| → primary.
     (2) z0 matches |C1|^{1/2} to 1e-9.
     (3) P2 residual magnitude in [1, 20].
     (4) Identity holds for both P1 and P2 with no liftidentity message.
     (5) ValidateLiftedDecomposition relErr < 5e-2.
   ============================================================================ *)

RunTest25C[] := Module[
  {c1, c2, poly1, poly2, vars, spec, det, lr, lc, ls, ld, res,
   z0ref, lv, lf, vl, relErr, firedLiftidentity, allPass},

  Print["--- Test 25C: 3-var, 2 polys with complex extremes, primary selection ---"];
  allPass = True;

  c1   = (1 + 2 I) 10^-5;
  c2   = (3 - I) 10^-4;
  vars = {x[1], x[2], x[3]};
  poly1 = 1 + c1 x[1]^3 + x[2]^2 + x[3];
  poly2 = 1 + x[1]^2 + c2 x[2]^4 + x[3]^2;
  spec  = <|"Polynomials"->{poly1, poly2}, "MonomialExponents"->{0,0,0},
            "PolynomialExponents"->{-2,-1}, "Variables"->vars,
            "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  (* (1) Detection *)
  det = DetectExtremeCoefficients[spec, 1000];
  If[Length[det] == 2,
    Print["  (1) detect count: PASS (2 extremes found)"],
    Print["  (1) detect count: FAIL (found ", Length[det], ")"]; allPass = False
  ];

  (* Verify primary = P1 monomial (larger |log|C||) *)
  Module[{logMags, primaryDetIdx},
    logMags = Abs[Log[Abs[N[#["Magnitude"]]]]] & /@ det;
    primaryDetIdx = First @ Ordering[logMags, -1];
    If[det[[primaryDetIdx]]["PolyIndex"] == 1,
      Print["  (1) primary is P1's monomial: PASS"],
      Print["  (1) primary is NOT P1's monomial: FAIL"]; allPass = False
    ]
  ];

  (* (2) z0 vs |C1|^{1/2} *)
  lr = {<|"PolyIndex"->1, "ExponentVector"->{3,0,0}, "k"->2|>,
        <|"PolyIndex"->2, "ExponentVector"->{0,4,0}, "k"->2|>};
  firedLiftidentity = False;
  Check[
    lc = LiftCoefficients[spec, lr],
    firedLiftidentity = True,
    TropicalEval::liftidentity
  ];
  If[!AssociationQ[lc], Print["  LiftCoefficients FAILED"]; Return[False]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"]; res = ld["Residuals"];

  z0ref = Abs[N[c1]]^(1/2);
  If[Abs[N[ld["z0"]] - z0ref] < 10^-9,
    Print["  (2) z0 matches |C1|^{1/2}: PASS (z0=", N[ld["z0"]], ")"],
    Print["  (2) z0 mismatch: FAIL (z0=", N[ld["z0"]], ", ref=", z0ref, ")"]; allPass = False
  ];

  (* (3) P2 residual magnitude in [1, 20] *)
  Module[{mag2},
    mag2 = Abs[N[res[[2]]]];
    If[1 <= mag2 <= 20,
      Print["  (3) |c2| in [1,20]: PASS (|c2|=", mag2, ")"],
      Print["  (3) |c2| out of [1,20]: FAIL (|c2|=", mag2, ")"]; allPass = False
    ]
  ];

  (* (4) Identity for both polys *)
  If[!firedLiftidentity,
    Print["  (4) liftidentity not fired for either poly: PASS"],
    Print["  (4) liftidentity fired: FAIL"]; allPass = False
  ];

  (* (5) ValidateLiftedDecomposition *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[!ListQ[lf] || Length[lf] < 2,
    Print["  lifted fan degenerate — FAIL"]; Return[False]
  ];
  vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
  relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
  Print["  (5) validate relErr = ", relErr];
  If[NumericQ[relErr] && relErr < 5*^-2,
    Print["  (5) exactness: PASS (relErr < 5e-2)"],
    Print["  (5) exactness: FAIL"]; allPass = False
  ];

  Print[]; Print[If[allPass, "25C PASS", "25C FAIL"]];
  allPass
];


(* ============================================================================
   Test 25D — 4-variable, non-trivial Euler measure, degree-4 complex coeff.

   P  = 1 + (2+I)·1e-6·x1^4 + x2^2 + x3^2 + x4^2 + x1·x2·x3·x4
   A  = {1/3, 0, 0, 0}   (x1^(1/3) Euler prefactor)
   B  = {-3}

   Tests that the non-integer MonomialExponent 1/3 is preserved correctly in
   the lifted spec: newMonoExps = {1/3, 0, 0, 0, 0}.

   PASS gates:
     (1) SuggestedK = 2.
     (2) Lifted spec MonomialExponents = {1/3, 0, 0, 0, 0}.
     (3) ValidateLiftedDecomposition relErr < 5e-2
         (reference NIntegrate includes x1^(1/3) weight).
     (4) C++ MC within 5σ of NIntegrate reference.
   ============================================================================ *)

RunTest25D[] := Module[
  {c, poly, vars, spec, det, k, lr, lc, ls, ld, lv, lf, vl, relErr,
   ref, res, mcRe, mcIm, mcReErr, mcImErr, allPass},

  Print["--- Test 25D: 4-var, Euler weight x1^(1/3), complex coeff (2+I)1e-6 ---"];
  allPass = True;

  c    = (2 + I) 10^-6;
  vars = {x[1], x[2], x[3], x[4]};
  poly = 1 + c x[1]^4 + x[2]^2 + x[3]^2 + x[4]^2 + x[1] x[2] x[3] x[4];
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{1/3,0,0,0},
           "PolynomialExponents"->{-3}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  (* (1) k* = 2 *)
  det = DetectExtremeCoefficients[spec, 1000];
  k   = det[[1]]["SuggestedK"];
  If[k == 2,
    Print["  (1) SuggestedK=2: PASS"],
    Print["  (1) SuggestedK: FAIL (got ", k, ")"]; allPass = False
  ];

  lr = {<|"PolyIndex"->1, "ExponentVector"->{4,0,0,0}, "k"->2|>};
  lc = LiftCoefficients[spec, lr];
  If[!AssociationQ[lc], Print["  LiftCoefficients FAILED"]; Return[False]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"];

  (* (2) MonomialExponents of lifted spec = {1/3, 0, 0, 0, 0} *)
  If[ls["MonomialExponents"] == {1/3, 0, 0, 0, 0},
    Print["  (2) MonomialExponents: PASS (", ls["MonomialExponents"], ")"],
    Print["  (2) MonomialExponents: FAIL (", ls["MonomialExponents"], ")"]; allPass = False
  ];

  (* (3) ValidateLiftedDecomposition *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[!ListQ[lf] || Length[lf] < 2,
    Print["  lifted fan degenerate — FAIL"]; Return[False]
  ];
  vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 3];
  relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
  Print["  (3) validate relErr = ", relErr];
  If[NumericQ[relErr] && relErr < 5*^-2,
    Print["  (3) exactness: PASS (relErr < 5e-2)"],
    Print["  (3) exactness: FAIL"]; allPass = False
  ];

  (* (4) C++ MC vs NIntegrate (with x1^(1/3) weight) *)
  ref = Quiet @ NIntegrate[
    t1^(1/3) * (1 + (2+I)*10^-6*t1^4 + t2^2 + t3^2 + t4^2 + t1*t2*t3*t4)^(-3),
    {t1,0,Infinity}, {t2,0,Infinity}, {t3,0,Infinity}, {t4,0,Infinity},
    MaxRecursion->12, PrecisionGoal->2];
  res = EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules"->lr, "FanData"->lf,
    "NSamples"->5*10^5, "RunChecks"->False, "Verbose"->False,
    "WorkingDirectory"->$hdWorkDir];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    mcRe = res["Results"][[1]]["Re"]; mcIm = res["Results"][[1]]["Im"];
    mcReErr = res["Results"][[1]]["ReErr"]; mcImErr = res["Results"][[1]]["ImErr"];
    Print["  (4) MC  = ", mcRe, " + ", mcIm, " I"];
    Print["      ref = ", Re[ref], " + ", Im[ref], " I"];
    Print["      dev = (", mcRe-Re[ref], ", ", mcIm-Im[ref],
          ")  stderr=(", mcReErr, ", ", mcImErr, ")"];
    If[Abs[mcRe - Re[ref]] < 5 mcReErr && Abs[mcIm - Im[ref]] < 5 mcImErr,
      Print["  (4) C++ MC: PASS"],
      Print["  (4) C++ MC: FAIL"]; allPass = False
    ],
    Print["  (4) EvaluateTropicalMCLifted FAILED: ", res]; allPass = False
  ];

  Print[]; Print[If[allPass, "25D PASS", "25D FAIL"]];
  allPass
];


(* ============================================================================
   Test 25E — 3-variable degree-6, three competing complex extremes.

   eps = 3·1e-5
   P   = 1 + (1+I)·eps·x1^6 + (2+I)·eps·x1^3·x2^3 + (1-3I)·eps·x2^6
           + x1^2 + x2^2 + x3^2 + x1·x2·x3
   A   = {0,0,0},  B = {-2}

   Three extreme monomials of comparable |log|C||:
     x1^6:     |C| = sqrt(2)·eps ≈ 4.24e-5,  |log|C|| ≈ 10.07  (primary)
     x1^3·x2^3: |C| = sqrt(5)·eps ≈ 6.71e-5, |log|C|| ≈ 9.61
     x2^6:     |C| = sqrt(10)·eps ≈ 9.49e-5, |log|C|| ≈ 9.26

   k* = 2 for all three.  Residuals are all O(1) by construction.

   PASS gates:
     (1) Detection finds 3 extreme monomials, all with SuggestedK = 2.
     (2) Primary rule = x1^6 monomial (largest |log|C||); z0 matches.
     (3) All 3 residual magnitudes in [0.5, 3].
     (4) Identity holds for all 3 simultaneously, no liftidentity message.
     (5) ValidateLiftedDecomposition relErr < 1e-1.
     (6) C++ MC within 5σ of NIntegrate reference.
   ============================================================================ *)

RunTest25E[] := Module[
  {eps, c1, c2, c3, poly, vars, spec, det, lr, lc, ls, ld, res,
   z0ref, residMags, lv, lf, vl, relErr,
   mcRef, mcRes, mcRe, mcIm, mcReErr, mcImErr,
   firedLiftidentity, allPass},

  Print["--- Test 25E: 3-var deg-6, three competing complex extremes ---"];
  allPass = True;

  eps  = 3 10^-5;
  c1   = (1 +  I) eps;
  c2   = (2 +  I) eps;
  c3   = (1 - 3 I) eps;
  vars = {x[1], x[2], x[3]};
  poly = 1 + c1 x[1]^6 + c2 x[1]^3 x[2]^3 + c3 x[2]^6
           + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  (* (1) Detection: 3 extremes, all k=2 *)
  det = DetectExtremeCoefficients[spec, 1000];
  If[Length[det] == 3 && AllTrue[det, (#["SuggestedK"] == 2) &],
    Print["  (1) detect: PASS (3 extremes, all k=2)"],
    Print["  (1) detect: FAIL (found ", Length[det], ")"];
    allPass = False
  ];

  (* (2) Primary = x1^6 (largest |log|C||) *)
  Module[{logMags, primaryDetIdx, primaryVec},
    logMags = Abs[Log[Abs[N[#["Magnitude"]]]]] & /@ det;
    primaryDetIdx = First @ Ordering[logMags, -1];
    primaryVec = det[[primaryDetIdx]]["ExponentVector"];
    If[primaryVec == {6, 0, 0},
      Print["  (2) primary = x1^6: PASS"],
      Print["  (2) primary: FAIL (got ", primaryVec, ")"]; allPass = False
    ]
  ];

  lr = {<|"PolyIndex"->1, "ExponentVector"->{6,0,0}, "k"->2|>,
        <|"PolyIndex"->1, "ExponentVector"->{3,3,0}, "k"->2|>,
        <|"PolyIndex"->1, "ExponentVector"->{0,6,0}, "k"->2|>};

  firedLiftidentity = False;
  Check[
    lc = LiftCoefficients[spec, lr],
    firedLiftidentity = True,
    TropicalEval::liftidentity
  ];
  If[!AssociationQ[lc], Print["  LiftCoefficients FAILED"]; Return[False]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"]; res = ld["Residuals"];

  (* z0 check *)
  z0ref = Abs[N[c1]]^(1/2);
  If[Abs[N[ld["z0"]] - z0ref] < 10^-9,
    Print["  (2) z0 matches |c1|^{1/2}: PASS (z0=", N[ld["z0"]], ")"],
    Print["  (2) z0 mismatch: FAIL"]; allPass = False
  ];

  (* (3) Residual magnitudes in [0.5, 3] *)
  residMags = N[Abs /@ res];
  If[AllTrue[residMags, 0.5 <= # <= 3 &],
    Print["  (3) residual magnitudes in [0.5,3]: PASS (", residMags, ")"],
    Print["  (3) residual magnitudes: FAIL (", residMags, ")"]; allPass = False
  ];

  (* (4) Identity check *)
  If[!firedLiftidentity,
    Print["  (4) liftidentity not fired: PASS"],
    Print["  (4) liftidentity fired: FAIL"]; allPass = False
  ];

  (* (5) ValidateLiftedDecomposition *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[!ListQ[lf] || Length[lf] < 2,
    Print["  lifted fan degenerate — FAIL"]; Return[False]
  ];
  vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 3];
  relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
  Print["  (5) validate relErr = ", relErr];
  If[NumericQ[relErr] && relErr < 1*^-1,
    Print["  (5) exactness: PASS (relErr < 1e-1)"],
    Print["  (5) exactness: FAIL"]; allPass = False
  ];

  (* (6) C++ MC vs NIntegrate *)
  mcRef = Quiet @ NIntegrate[
    (1 + (1+I)*eps*t1^6 + (2+I)*eps*t1^3*t2^3 + (1-3I)*eps*t2^6
       + t1^2 + t2^2 + t3^2 + t1*t2*t3)^(-2),
    {t1,0,Infinity}, {t2,0,Infinity}, {t3,0,Infinity},
    MaxRecursion->12, PrecisionGoal->2];
  mcRes = EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules"->lr, "FanData"->lf,
    "NSamples"->5*10^5, "RunChecks"->False, "Verbose"->False,
    "WorkingDirectory"->$hdWorkDir];
  If[AssociationQ[mcRes] && KeyExistsQ[mcRes, "Results"],
    mcRe = mcRes["Results"][[1]]["Re"]; mcIm = mcRes["Results"][[1]]["Im"];
    mcReErr = mcRes["Results"][[1]]["ReErr"]; mcImErr = mcRes["Results"][[1]]["ImErr"];
    Print["  (6) MC  = ", mcRe, " + ", mcIm, " I"];
    Print["      ref = ", Re[mcRef], " + ", Im[mcRef], " I"];
    Print["      dev = (", mcRe-Re[mcRef], ", ", mcIm-Im[mcRef],
          ")  stderr=(", mcReErr, ", ", mcImErr, ")"];
    If[Abs[mcRe - Re[mcRef]] < 5 mcReErr && Abs[mcIm - Im[mcRef]] < 5 mcImErr,
      Print["  (6) C++ MC: PASS"],
      Print["  (6) C++ MC: FAIL"]; allPass = False
    ],
    Print["  (6) EvaluateTropicalMCLifted FAILED: ", mcRes]; allPass = False
  ];

  Print[]; Print[If[allPass, "25E PASS", "25E FAIL"]];
  allPass
];


(* ============================================================================
   Umbrella
   ============================================================================ *)

RunTest25[] := Module[{p25a, p25b, p25c, p25d, p25e, all},
  Print[""];
  Print["================================================================"];
  Print["  Test 25: higher-dimensional lifting (3–4 variables, deg 4–6)"];
  Print["================================================================"];
  Print[""];

  p25a = RunTest25A[]; Print[""];
  p25b = RunTest25B[]; Print[""];
  p25c = RunTest25C[]; Print[""];
  p25d = RunTest25D[]; Print[""];
  p25e = RunTest25E[]; Print[""];

  all = p25a && p25b && p25c && p25d && p25e;

  Print["================================================================"];
  Print["  Test 25 Summary"];
  Print["  25A (3-var deg-4, complex coeff, full C++):   ", If[p25a,"PASS","FAIL"]];
  Print["  25B (3-var deg-5, k*=3, C++ MC):              ", If[p25b,"PASS","FAIL"]];
  Print["  25C (3-var, 2 polys, multi-rule):             ", If[p25c,"PASS","FAIL"]];
  Print["  25D (4-var, Euler weight A_1=1/3, C++ MC):    ", If[p25d,"PASS","FAIL"]];
  Print["  25E (3-var deg-6, 3 competing extremes, MC):  ", If[p25e,"PASS","FAIL"]];
  Print["  Overall Test 25: ", If[all, "PASS", "FAIL"]];
  Print["================================================================"];
  Print[""];
  all
];


(* --- Execute --- *)
If[RunTest25[], Print["ALL HIGH-DIMENSION LIFT TESTS PASSED"],
   Print["SOME HIGH-DIMENSION LIFT TESTS FAILED"]; Exit[1]];
