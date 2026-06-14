(* ============================================================================
   test_complex_lifted.wl — Test 24: auxiliary-variable lifting with COMPLEX,
   small-magnitude polynomial coefficients.

   The lifting module (DetectExtremeCoefficients / LiftCoefficients /
   ProcessSectorLifted / ValidateLiftedDecomposition / EvaluateTropicalMCLifted)
   was originally exercised only on real coefficients.  These tests confirm it
   is correct when a polynomial carries a COMPLEX coefficient whose magnitude is
   very small — the magnitude is moved into the real anchor z0 = |C|^{1/k} and
   the PHASE is carried exactly by the residual c = C/z0^k (an O(1) complex
   number), so the lifted decomposition stays exact and the C++ Monte Carlo
   returns the correct complex value.

   Subtests (each returns True/False):
     RunTestCxA[] — single small complex coeff (1+I)1e-4, k=2: exactness
                    (Re & Im), end-to-end C++ MC, and clean complex codegen.
     RunTestCxB[] — very small (1+I)1e-7 with the AUTOMATIC anchor (k*=3):
                    the kStar rule scales with |C|; lifted decomposition exact.
     RunTestCxC[] — TWO small complex coeffs lifted against one anchor
                    (multi-rule): residual phases are O(1); decomposition exact.
     RunTestCxD[] — boundary/regression: a complex polynomial EXPONENT (not a
                    coefficient) must still fire TropicalEval::liftcomplex.

   Umbrella RunTestCx[] runs all four, prints PASS/FAIL, and Exit[]s non-zero
   if any part fails.

   Requires: tropical_fan.wl, tropical_eval.wl, Polymake; g++ for part A.
   Run:
     cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file EXAMPLES/test_complex_lifted.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Print["tropical_eval.wl loaded"];
Print[];

$cxWorkDir = FileNameJoin[{Directory[], "INTERFILES"}];


(* ============================================================================
   Test 24A: single small complex coefficient, full pipeline.

   I = Int_[0,inf)^2 dx1 dx2 / (1 + (1+I) 1e-4 x1^2 + x2^2 + x1 x2^2)^2
   (same Newton-polytope structure as the real Case A / Example 20, but the
   x1^2 coefficient is now complex of magnitude |C| = Sqrt[2] 1e-4 ~ 1.41e-4).

   Lift C x1^2 -> z^2 x1^2, k=2  =>  z0 = |C|^{1/2} ~ 0.0119 (real),
   residual c = C/z0^2 = (1+I)/Sqrt[2]  (unit modulus phase, O(1)).

   PASS gates:
     (1) DetectExtremeCoefficients flags it with SuggestedK == 2.
     (2) LiftCoefficients identity (lifted poly at z->z0 equals original) holds,
         and |residual| == 1.
     (3) ValidateLiftedDecomposition relErr < 1e-3 AND the sector sum matches
         direct NIntegrate in BOTH Re and Im.
     (4) EvaluateTropicalMCLifted (C++, 1e6 samples) lands within 5 sigma of the
         NIntegrate reference in BOTH Re and Im.
     (5) the generated C++ contains a cx(...) complex literal and no broken
         CForm token (Complex(...)).
   ============================================================================ *)

RunTestCxA[] := Module[
  {c, poly, vars, spec, lr, det, lc, ls, ld, resid,
   lv, lf, vl, relErr, ref, res, mcRe, mcIm, mcReErr, mcImErr,
   src, idOK, allPass},

  Print["--- Test 24A: single small complex coeff (1+I)1e-4, k=2 ---"];
  allPass = True;

  c    = (1 + I) 10^-4;
  vars = {x[1], x[2]};
  poly = 1 + c x[1]^2 + x[2]^2 + x[1] x[2]^2;
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;
  lr   = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

  (* (1) detection *)
  det = DetectExtremeCoefficients[spec, 1000];
  If[Length[det] == 1 && det[[1]]["SuggestedK"] == 2 &&
     Abs[det[[1]]["Magnitude"] - Sqrt[2.] 10^-4] < 10^-9,
    Print["  (1) detect: PASS (SuggestedK=2, |C|=", det[[1]]["Magnitude"], ")"],
    Print["  (1) detect: FAIL -> ", det]; allPass = False
  ];

  (* (2) lift + identity + unit residual *)
  lc = LiftCoefficients[spec, lr];
  If[!AssociationQ[lc], Print["  (2) LiftCoefficients FAILED"]; Return[False]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"]; resid = ld["Residuals"][[1]];
  idOK = Simplify[(ls["Polynomials"][[1]] /. ld["AuxVariable"]->ld["z0"]) - poly] === 0;
  If[idOK && Abs[N[Abs[resid]] - 1] < 10^-9,
    Print["  (2) lift: PASS (identity holds; z0=", N[ld["z0"]],
          ", residual=", N[resid], ", |residual|=1)"],
    Print["  (2) lift: FAIL (identity=", idOK, ", residual=", N[resid], ")"];
    allPass = False
  ];

  (* (3) exactness vs direct NIntegrate (complex) *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = ComputeDecomposition[lv, "ShowProgress"->False];
  vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
  If[AssociationQ[vl],
    relErr = vl["RelativeError"];
    Print["  (3) validate: direct = ", vl["DirectResult"]];
    Print["               sum    = ", vl["SectorSum"]];
    Print["               relErr = ", relErr];
    If[NumericQ[relErr] && relErr < 10^-3 &&
       Abs[Re[vl["SectorSum"] - vl["DirectResult"]]] < 10^-3 Abs[vl["DirectResult"]] &&
       Abs[Im[vl["SectorSum"] - vl["DirectResult"]]] < 10^-3 Abs[vl["DirectResult"]],
      Print["  (3) exactness: PASS (relErr < 1e-3, Re & Im both match)"],
      Print["  (3) exactness: FAIL (relErr=", relErr, ")"]; allPass = False
    ],
    Print["  (3) ValidateLiftedDecomposition FAILED"]; allPass = False
  ];

  (* (4) end-to-end C++ MC, Re and Im within 5 sigma *)
  ref = Quiet @ NIntegrate[1/(1 + c t1^2 + t2^2 + t1 t2^2)^2,
    {t1,0,Infinity}, {t2,0,Infinity}, MaxRecursion->25, PrecisionGoal->8];
  res = EvaluateTropicalMCLifted[spec, {{}}, "LiftRules"->lr,
    "NSamples"->10^6, "RunChecks"->False, "Verbose"->False,
    "WorkingDirectory"->$cxWorkDir];
  If[AssociationQ[res] && KeyExistsQ[res, "Results"],
    mcRe = res["Results"][[1]]["Re"];     mcIm = res["Results"][[1]]["Im"];
    mcReErr = res["Results"][[1]]["ReErr"]; mcImErr = res["Results"][[1]]["ImErr"];
    Print["  (4) MC = ", mcRe, " + ", mcIm, " I"];
    Print["      ref= ", Re[ref], " + ", Im[ref], " I"];
    Print["      dev= (", mcRe - Re[ref], ", ", mcIm - Im[ref],
          ")   stderr=(", mcReErr, ", ", mcImErr, ")"];
    If[Abs[mcRe - Re[ref]] < 5 mcReErr && Abs[mcIm - Im[ref]] < 5 mcImErr,
      Print["  (4) C++ MC: PASS (Re & Im within 5 sigma)"],
      Print["  (4) C++ MC: FAIL (outside 5 sigma)"]; allPass = False
    ],
    Print["  (4) EvaluateTropicalMCLifted FAILED (is g++ installed?): ", res];
    allPass = False
  ];

  (* (5) generated C++ uses a complex literal and has no broken CForm token *)
  src = Quiet @ Import[FileNameJoin[{$cxWorkDir,"tropical_mc_generated.cpp"}], "Text"];
  If[StringQ[src] && StringContainsQ[src, "cx("] && !StringContainsQ[src, "Complex("],
    Print["  (5) codegen: PASS (cx(...) present, no Complex(...) token)"],
    Print["  (5) codegen: FAIL"]; allPass = False
  ];

  Print[]; Print[If[allPass, "24A PASS", "24A FAIL"]];
  allPass
];


(* ============================================================================
   Test 24B: very small magnitude with the AUTOMATIC anchor.

   c = (1+I) 1e-7  =>  |C| ~ 1.41e-7,  k* = Ceiling[|log10|C||/3] = 3,
   z0 = |C|^{1/3} ~ 5.2e-3.  Confirms the kStar rule scales with the (complex)
   magnitude and that the lifted decomposition stays exact at k=3.

   PASS gates: SuggestedK == 3; ValidateLiftedDecomposition relErr < 1e-3.
   (If the k=3 lifted polytope were degenerate the automatic fan would need an
   explicit FanData — reported, not silently passed.)
   ============================================================================ *)

RunTestCxB[] := Module[
  {c, poly, vars, spec, det, k, lr, lc, ls, ld, lv, lf, vl, relErr, allPass},

  Print["--- Test 24B: very small (1+I)1e-7, automatic anchor (expect k*=3) ---"];
  allPass = True;

  c    = (1 + I) 10^-7;
  vars = {x[1], x[2]};
  poly = 1 + c x[1]^2 + x[2]^2 + x[1] x[2]^2;
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  det = DetectExtremeCoefficients[spec, 1000];
  k   = det[[1]]["SuggestedK"];
  If[k == 3,
    Print["  SuggestedK: PASS (k* = 3 for |C| ~ ", det[[1]]["Magnitude"], ")"],
    Print["  SuggestedK: FAIL (got ", k, ", expected 3)"]; allPass = False
  ];

  lr = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->k|>};
  lc = LiftCoefficients[spec, lr]; ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  Print["  z0 = ", N[ld["z0"]], "   residual = ", N[ld["Residuals"][[1]]]];

  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[ListQ[lf] && AllTrue[lf[[2]], Length[#] == 3 &],
    vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
    relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
    Print["  validate relErr = ", relErr];
    If[NumericQ[relErr] && relErr < 10^-3,
      Print["  exactness: PASS (relErr < 1e-3)"],
      Print["  exactness: FAIL"]; allPass = False
    ],
    Print["  lifted polytope degenerate at k=", k,
          " -> needs explicit FanData (liftdegenerate path); marking FAIL"];
    allPass = False
  ];

  Print[]; Print[If[allPass, "24B PASS", "24B FAIL"]];
  allPass
];


(* ============================================================================
   Test 24C: multi-rule lift — two small complex coefficients, one anchor.

   P = 1 + (1+I)1e-4 x1^2 + x2^2 + (1-2I)1e-4 x1 x2^2,  B = -2.
   Lift BOTH extreme monomials against the same auxiliary variable (k=2 each).
   The shared anchor z0 = |C_primary|^{1/2} fixes the magnitude; each residual
   c_i = C_i/z0^2 keeps its own O(1) phase, so the round-trip stays exact for
   several lifted monomials at once.

   PASS gates: identity holds; both |residual| are O(1) (in [0.1, 10]);
   ValidateLiftedDecomposition relErr < 1e-2.
   ============================================================================ *)

RunTestCxC[] := Module[
  {c1, c2, poly, vars, spec, lr, lc, ls, ld, res, lv, lf, vl, relErr,
   idOK, magsOK, allPass},

  Print["--- Test 24C: multi-rule, two small complex coeffs ---"];
  allPass = True;

  c1 = (1 + I) 10^-4;  c2 = (1 - 2 I) 10^-4;
  vars = {x[1], x[2]};
  poly = 1 + c1 x[1]^2 + x[2]^2 + c2 x[1] x[2]^2;
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;
  lr = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>,
        <|"PolyIndex"->1, "ExponentVector"->{1,2}, "k"->2|>};

  lc = LiftCoefficients[spec, lr];
  If[!AssociationQ[lc], Print["  LiftCoefficients FAILED"]; Return[False]];
  ls = lc["LiftedSpec"]; ld = lc["LiftData"]; res = ld["Residuals"];
  Print["  z0        = ", N[ld["z0"]]];
  Print["  residuals = ", N[res]];

  idOK   = Simplify[(ls["Polynomials"][[1]] /. ld["AuxVariable"]->ld["z0"]) - poly] === 0;
  magsOK = AllTrue[N[Abs /@ res], 0.1 <= # <= 10 &];
  If[idOK, Print["  identity: PASS"], Print["  identity: FAIL"]; allPass = False];
  If[magsOK, Print["  residual magnitudes O(1): PASS"],
             Print["  residual magnitudes O(1): FAIL -> ", N[Abs /@ res]]; allPass = False];

  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[ListQ[lf] && AllTrue[lf[[2]], Length[#] == 3 &],
    vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
    relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
    Print["  validate direct = ", If[AssociationQ[vl], vl["DirectResult"], "?"]];
    Print["  validate relErr = ", relErr];
    If[NumericQ[relErr] && relErr < 10^-2,
      Print["  exactness: PASS (relErr < 1e-2)"],
      Print["  exactness: FAIL"]; allPass = False
    ],
    Print["  multi-rule lifted polytope degenerate -> needs explicit FanData; FAIL"];
    allPass = False
  ];

  Print[]; Print[If[allPass, "24C PASS", "24C FAIL"]];
  allPass
];


(* ============================================================================
   Test 24D: boundary/regression — complex polynomial EXPONENT.

   Complex *coefficients* are supported; complex *exponents* B are not (a
   complex effective exponent admits no real-valued domain indicator).  With a
   real small coefficient but B = -(2 + I/2), ProcessSectorLifted must fire
   TropicalEval::liftcomplex and return $Failed.  This pins the supported
   envelope and guards against the two cases being conflated.

   PASS gate: liftcomplex fires (some sector returns $Failed).
   ============================================================================ *)

RunTestCxD[] := Module[
  {poly, vars, spec, lr, lc, ls, ld, lv, lf, sd, fired, allPass},

  Print["--- Test 24D: complex EXPONENT must fire liftcomplex (regression) ---"];
  allPass = True;

  vars = {x[1], x[2]};
  poly = 1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2;     (* real coeff: lift triggers *)
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-(2 + I/2)},        (* complex EXPONENT *)
           "Variables"->vars, "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;
  lr = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

  lc = LiftCoefficients[spec, lr]; ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = ComputeDecomposition[lv, "ShowProgress"->False];

  fired = False;
  Check[
    Do[ sd = ProcessSectorLifted[ls, lf[[1]], lf[[2,s]], s, ld];
        If[sd === $Failed, fired = True],
      {s, Length[lf[[2]]]}],
    fired = True,
    TropicalEval::liftcomplex
  ];

  If[fired,
    Print["  liftcomplex: PASS (fired for complex exponent, $Failed returned)"],
    Print["  liftcomplex: FAIL (did not fire)"]; allPass = False
  ];

  Print[]; Print[If[allPass, "24D PASS", "24D FAIL"]];
  allPass
];


(* ============================================================================
   Umbrella
   ============================================================================ *)

RunTestCx[] := Module[{pa, pb, pc, pd, all},
  Print[""];
  Print["================================================================"];
  Print["  Test 24: lifting with COMPLEX small-magnitude coefficients"];
  Print["================================================================"];
  Print[""];

  pa = RunTestCxA[]; Print[""];
  pb = RunTestCxB[]; Print[""];
  pc = RunTestCxC[]; Print[""];
  pd = RunTestCxD[]; Print[""];

  all = pa && pb && pc && pd;

  Print["================================================================"];
  Print["  Test 24 Summary"];
  Print["  24A (single complex coeff, full C++):  ", If[pa, "PASS", "FAIL"]];
  Print["  24B (very small, automatic k*=3):      ", If[pb, "PASS", "FAIL"]];
  Print["  24C (multi-rule complex residuals):    ", If[pc, "PASS", "FAIL"]];
  Print["  24D (complex exponent -> liftcomplex): ", If[pd, "PASS", "FAIL"]];
  Print["  Overall Test 24: ", If[all, "PASS", "FAIL"]];
  Print["================================================================"];
  Print[""];
  all
];


(* --- Execute --- *)
If[RunTestCx[], Print["ALL COMPLEX-LIFT TESTS PASSED"],
   Print["SOME COMPLEX-LIFT TESTS FAILED"]; Exit[1]];
