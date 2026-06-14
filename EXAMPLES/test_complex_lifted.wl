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
     RunTestCxE[] — Fix 1 regression: real floating-point extreme coeff must NOT
                    trip a false-positive liftidentity.
     RunTestCxF[] — Fix 2+3: complex EXPONENT via SplitRealImag on a real-positive
                    base (arg P = 0); pre-flight + value vs reference.
     RunTestCxG[] — BUG 1 (TMCv2_BUG_LOG.md): SplitRealImag on a COMPLEX base
                    (arg P != 0).  The phase must use the COMPLEX log std::log(P),
                    not the modulus std::log|P|, else the magnitude factor
                    exp(-Im(B) arg P) (~40% of the integral here) is dropped.
     RunTestCxH[] — BUG 2 (TMCv2_BUG_LOG.md): lifting on RAW-SYMBOL base vars must
                    not build Symbol[k] (the Symbol::string footgun).

   Umbrella RunTestCx[] runs all eight, prints PASS/FAIL, and Exit[]s non-zero
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
   Test 24E: Fix 1 validation — real floating-point extreme coefficient.

   P = 3.7 + 2.3e-6 x[1]^2 + x[1] + x[2]^2,  A = {0,0},  B = {-2}.
   The constant term 3.7 and the tiny coefficient 2.3e-6 (a machine float)
   make z0 = sqrt(2.3e-6) a machine float.  Before Fix 1, Simplify cannot
   reduce the floating-point roundtrip error to 0 and the liftidentity message
   fires as a false positive.  (The linear x[1] term gives the lifted Newton
   polytope its full dimension -- x[1] then appears in two monomials, x[1] and
   x[1]^2 -- so ValidateLiftedDecomposition's automatic fan is non-degenerate;
   without it the 3-monomial polytope would be lower-dimensional.)

   PASS gates:
     (1) DetectExtremeCoefficients finds exactly 1 extreme entry (x[1]^2),
         SuggestedK == 2.
     (2) liftidentity does NOT fire.
     (3) ValidateLiftedDecomposition relErr < 1e-2.
   ============================================================================ *)

RunTestCxE[] := Module[
  {poly, vars, spec, lr, det, lc, ls, ld, lv, lf, vl, relErr,
   firedLiftidentity, allPass},

  Print["--- Test 24E: Fix 1 — real float extreme coeff (liftidentity false-positive) ---"];
  allPass = True;

  vars = {x[1], x[2]};
  poly = 3.7 + 2.3*^-6 * x[1]^2 + x[1] + x[2]^2;
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;

  (* (1) Detection: exactly one extreme monomial (the 2.3e-6 x[1]^2); the
     constant 3.7 and the unit coefficients of x[1], x[2]^2 are non-extreme *)
  det = DetectExtremeCoefficients[spec, 1000];
  If[Length[det] == 1 &&
     AllTrue[det, (#["ExponentVector"] == {2,0} && #["SuggestedK"] == 2) &],
    Print["  (1) detect: PASS (1 extreme, x[1]^2, k=2)"],
    Print["  (1) detect: FAIL -> ", det]; allPass = False
  ];

  lr = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

  (* (2) LiftCoefficients must NOT fire liftidentity *)
  firedLiftidentity = False;
  Check[
    lc = LiftCoefficients[spec, lr],
    firedLiftidentity = True,
    TropicalEval::liftidentity
  ];
  If[!firedLiftidentity && AssociationQ[lc],
    Print["  (2) liftidentity not fired: PASS"],
    Print["  (2) liftidentity fired (false positive): FAIL"]; allPass = False
  ];
  If[!AssociationQ[lc], Print["  LiftCoefficients FAILED"]; Return[False]];

  ls = lc["LiftedSpec"]; ld = lc["LiftData"];
  Print["  z0 = ", N[ld["z0"]], "  residual = ", N[ld["Residuals"][[1]]]];

  (* (3) ValidateLiftedDecomposition *)
  lv = PolytopeVertices[(Times @@ ls["Polynomials"])^(-1), ls["Variables"]];
  lf = Quiet @ ComputeDecomposition[lv, "ShowProgress"->False];
  If[ListQ[lf] && AllTrue[lf[[2]], Length[#] == 3 &],
    vl = Quiet @ ValidateLiftedDecomposition[spec, ls, lf, ld, {}, 4];
    relErr = If[AssociationQ[vl], vl["RelativeError"], $Failed];
    Print["  validate relErr = ", relErr];
    If[NumericQ[relErr] && relErr < 1*^-2,
      Print["  (3) exactness: PASS (relErr < 1e-2)"],
      Print["  (3) exactness: FAIL"]; allPass = False
    ],
    Print["  lifted polytope degenerate -> FAIL"]; allPass = False
  ];

  Print[]; Print[If[allPass, "24E PASS", "24E FAIL"]];
  allPass
];


(* ============================================================================
   Test 24F: Fix 2 + Fix 3 — complex polynomial exponent with SplitRealImag.

   P = 1 + 1e-6 x[1]^2 + x[2]^2 + x[1] x[2]^2,  A = {0,0},
   B = -(2+I)  (complex exponent: gamma_re = -2, gamma_im = -1).

   With "ComplexExponentMode" -> "Reject" (default), the pre-flight check fires
   TropicalEval::liftcomplexexponents and returns $Failed.
   With "ComplexExponentMode" -> "SplitRealImag", the integrand is split into
   real-exponent weight P^{-2} and oscillatory phase exp(-I log|P|); lifting
   proceeds on the real exponent and the phase is multiplied per MC sample.

   PASS gates:
     (1) With Reject mode, liftcomplexexponents fires and $Failed is returned.
     (2) With SplitRealImag mode, liftcomplexexponents does NOT fire and the
         call returns an Association (not $Failed).
     (3) MC result within 5σ of NIntegrate reference (complex value).
   ============================================================================ *)

RunTestCxF[] := Module[
  {poly, vars, spec, lr, ref, res, mcRe, mcIm, mcReErr, mcImErr,
   firedPreflight, allPass, cubaHere, resV, vRe, vIm, relErrV},

  Print["--- Test 24F: Fix 2+3 — complex exponent, SplitRealImag mode ---"];
  allPass = True;

  vars = {x[1], x[2]};
  poly = 1 + 10^-6 * x[1]^2 + x[2]^2 + x[1] * x[2]^2;
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-(2 + I)}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;
  lr = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

  (* (1) Reject mode (default) should fire liftcomplexexponents and return $Failed *)
  firedPreflight = False;
  Quiet[
    Check[
      EvaluateTropicalMCLifted[spec, {{}},
        "LiftRules"->lr, "NSamples"->100,
        "RunChecks"->False, "Verbose"->False,
        "WorkingDirectory"->$cxWorkDir],
      firedPreflight = True,
      TropicalEval::liftcomplexexponents
    ],
    TropicalEval::liftcomplexexponents
  ];
  If[firedPreflight,
    Print["  (1) Reject mode fires liftcomplexexponents: PASS"],
    Print["  (1) Reject mode does not fire (should have): FAIL"]; allPass = False
  ];

  (* (2) SplitRealImag mode — pre-flight must NOT fire *)
  firedPreflight = False;
  res = Quiet[
    Check[
      EvaluateTropicalMCLifted[spec, {{}},
        "LiftRules"->lr,
        "ComplexExponentMode"->"SplitRealImag",
        "NSamples"->10^6, "RunChecks"->False, "Verbose"->False,
        "WorkingDirectory"->$cxWorkDir],
      (firedPreflight = True; $Failed),
      TropicalEval::liftcomplexexponents
    ],
    TropicalEval::liftcomplexexponents
  ];
  If[!firedPreflight,
    Print["  (2) SplitRealImag does not fire liftcomplexexponents: PASS"],
    Print["  (2) SplitRealImag fires liftcomplexexponents: FAIL"]; allPass = False
  ];

  (* (3) value check: SplitRealImag must reproduce the (complex) reference.
     The split puts the imaginary exponent into an oscillatory phase
     exp(i*Im(B)*log|P|); log|P| includes the tropical monomial factor cleared
     out of P, so the per-sector "MonoFactorLog" term is essential (omitting it
     -- the pre-fix behavior -- gave a wrong phase, ~50% off).  This lifted +
     OSCILLATORY integrand is heavy-tailed, so plain-MC error bars are optimistic
     (Manual sec:vegas); we therefore check the value with VEGAS at a larger
     budget against a high-precision reference, within the tolerance documented
     for the lifted+oscillatory regime.  (Needs CUBA; skipped if absent.) *)
  ref = Quiet @ NIntegrate[
    (1 + 10^-6 t1^2 + t2^2 + t1 t2^2)^(-(2+I)),
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion->40, PrecisionGoal->8, WorkingPrecision->30];
  cubaHere = AnyTrue[{"/opt/homebrew", "/usr/local", "/usr"},
    FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
    (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
     FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
     FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &];
  If[! (AssociationQ[res] && KeyExistsQ[res, "Results"]),
    If[!firedPreflight, Print["  (3) EvaluateTropicalMCLifted FAILED: ", res]; allPass = False],
    (* res ok: now the accurate value check *)
    If[! cubaHere,
      Print["  (3) value check SKIPPED (CUBA absent; VEGAS needed for the \
heavy-tail oscillatory integrand). Phase codegen still loaded."],
      resV = Quiet @ EvaluateTropicalMCLifted[spec, {{}},
        "LiftRules"->lr, "ComplexExponentMode"->"SplitRealImag",
        "Integrator"->"VEGAS", "NSamples"->10^7, "VegasEpsRel"->1.*^-9,
        "RunChecks"->False, "Verbose"->False, "WorkingDirectory"->$cxWorkDir];
      If[AssociationQ[resV] && KeyExistsQ[resV, "Results"],
        vRe = resV["Results"][[1]]["Re"]; vIm = resV["Results"][[1]]["Im"];
        relErrV = Abs[(vRe + I vIm) - ref]/Abs[ref];
        Print["  (3) VEGAS = ", vRe, " + ", vIm, " I"];
        Print["      ref   = ", N[Re[ref]], " + ", N[Im[ref]], " I   relErr = ", N[relErrV]];
        If[NumericQ[relErrV] && relErrV < 5*^-2,
          Print["  (3) SplitRealImag value: PASS (relErr < 5e-2 vs reference; \
phase incl. monomial factor)"],
          Print["  (3) SplitRealImag value: FAIL (relErr=", N[relErrV], ")"]; allPass = False],
        Print["  (3) VEGAS run FAILED: ", resV]; allPass = False]
    ]
  ];

  Print[]; Print[If[allPass, "24F PASS", "24F FAIL"]];
  allPass
];


(* ============================================================================
   Test 24G: BUG 1 (TMCv2_BUG_LOG.md) — SplitRealImag phase on a COMPLEX base.

   24F uses a real-positive polynomial, so arg P = 0 on the domain and the
   oscillatory-phase magnitude factor exp(-Im(B) arg P) is identically 1 — the
   bug is invisible there.  CALC2's atoms (u x + v y) have COMPLEX coefficients,
   so arg P != 0 and that factor is real and order-one.  The pre-fix codegen
   emitted the phase as exp(i Im(B) log|P|) using the MODULUS log|P|, which keeps
   only the pure phase exp(i Im(B) log|P|) and silently DROPS exp(-Im(B) arg P).
   The fix uses the COMPLEX log std::log(P): log P = log|P| + i arg P, so
   exp(i Im(B) log P) = exp(-Im(B) arg P) * exp(i Im(B) log|P|) restores it.

   Integrand:  P = (1+I) + 1e-4 x1^2 + x2^2 + x1 x2^2 ,  B = -(2+I).
   The O(1) complex CONSTANT term gives arg P ~ pi/4 in the dominant region near
   the origin (Im P == 1, Re P >= 1, so |P| >= sqrt 2 — the integrand is bounded,
   not heavy-tailed, so MC and VEGAS both converge cleanly).  The small REAL
   1e-4 x1^2 coefficient triggers a clean real lift (z0 = 0.01, residual 1); the
   complex constant rides along as an ordinary complex coefficient.  The dropped
   magnitude factor is ~40% of the integral, so the pre-fix value is ~40% off —
   decisively outside the gates below.

   PASS gates:
     (1) generated C++ phase uses the COMPLEX log (cx(0.0,imB[..]) * std::log(P..))
         and contains NO modulus form std::log(std::abs(P..))  [the exact fix;
         deterministic, CUBA-independent].
     (2) plain-MC value matches NIntegrate[P^B] to relErr < 1.2e-1  [behavioral,
         CUBA-independent; pre-fix ~0.41].
     (3) VEGAS value matches NIntegrate[P^B] to relErr < 5e-2       [tight
         behavioral confirmation; needs CUBA, skipped if absent].
   ============================================================================ *)

RunTestCxG[] := Module[
  {cc, poly, vars, B, spec, lr, ref, resMC, mcRe, mcIm, relErrMC,
   src, hasFixed, hasBuggy, cubaHere, resV, vRe, vIm, relErrV, allPass},

  Print["--- Test 24G: BUG 1 — complex-COEFFICIENT base, SplitRealImag phase ---"];
  allPass = True;

  cc   = 1 + I;
  vars = {x[1], x[2]};
  poly = cc + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  B    = -(2 + I);
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{B}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;
  lr   = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

  (* true reference: the actual complex integrand P^B (principal branch).
     Re P >= 1 and Im P == 1 on the domain, so P stays in the right half-plane —
     no branch-cut crossing; std::log (C++) and Log (WL) use the same branch. *)
  ref = Quiet @ NIntegrate[(cc + 10^-4 t1^2 + t2^2 + t1 t2^2)^B,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion->40, PrecisionGoal->8, WorkingPrecision->30];
  Print["  reference NIntegrate[P^B] = ", N[Re[ref]], " + ", N[Im[ref]], " I"];

  (* Generate + run plain MC (needs only g++).  This also writes the C++ we
     inspect in gate (1).  Pre-fix this yields ~ -2.86 - 8.49 I (~41% off);
     post-fix ~ -4.43 - 14.42 I. *)
  resMC = EvaluateTropicalMCLifted[spec, {{}}, "LiftRules"->lr,
    "ComplexExponentMode"->"SplitRealImag", "NSamples"->10^6,
    "RunChecks"->False, "Verbose"->False, "WorkingDirectory"->$cxWorkDir];
  If[!(AssociationQ[resMC] && KeyExistsQ[resMC, "Results"]),
    Print["  EvaluateTropicalMCLifted (plain MC) FAILED (is g++ installed?): ", resMC];
    Print[]; Print["24G FAIL"]; Return[False]
  ];

  (* (1) codegen structure — the exact line BUG 1 fixed.  Two genuinely
     direction-discriminating markers (verified against the pre-fix build):
       fixed:  "cx(0.0, imB[0]) * "  -- the cx is CLOSED (i folded) before the
               multiply, so the phase exponent is cx(0,imB)*(.. + std::log(P..)).
       buggy:  "std::log(std::abs(P" -- the modulus log (pure phase) the fix
               removed; pre-fix the term was cx(0.0, imB[0] * (.. log|P|)).
     (The looser prefix "cx(0.0, imB[" appears in BOTH, since pre-fix wrapped the
     whole real phaseSum as std::exp(cx(0.0, imB[0]*...)) — do not key on it.) *)
  src = Quiet @ Import[FileNameJoin[{$cxWorkDir, "tropical_mc_generated.cpp"}], "Text"];
  hasFixed = StringQ[src] && StringContainsQ[src, "cx(0.0, imB[0]) *"];
  hasBuggy = StringQ[src] && StringContainsQ[src, "std::log(std::abs(P"];
  If[hasFixed && !hasBuggy,
    Print["  (1) codegen: PASS (phase = cx(0.0,imB[..]) * (.. + std::log(P..)); complex log, no log|P|)"],
    Print["  (1) codegen: FAIL (complex-log form present=", hasFixed,
          ", modulus std::log(std::abs(P..)) present=", hasBuggy, ")"];
    allPass = False
  ];

  (* (2) plain-MC value (CUBA-independent behavioral check) *)
  mcRe = resMC["Results"][[1]]["Re"]; mcIm = resMC["Results"][[1]]["Im"];
  relErrMC = Abs[(mcRe + I mcIm) - ref]/Abs[ref];
  Print["  (2) MC 1e6 = ", mcRe, " + ", mcIm, " I   relErr = ", N[relErrMC]];
  If[NumericQ[relErrMC] && relErrMC < 1.2*^-1,
    Print["  (2) MC value: PASS (relErr < 0.12 vs complex P^B; pre-fix ~0.41)"],
    Print["  (2) MC value: FAIL (relErr=", N[relErrMC], ")"]; allPass = False
  ];

  (* (3) VEGAS value (tight; needs CUBA) *)
  cubaHere = AnyTrue[{"/opt/homebrew", "/usr/local", "/usr"},
    FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
    (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
     FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
     FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &];
  If[! cubaHere,
    Print["  (3) VEGAS value SKIPPED (CUBA absent; complex-log codegen already \
verified structurally in gate (1))."],
    resV = Quiet @ EvaluateTropicalMCLifted[spec, {{}}, "LiftRules"->lr,
      "ComplexExponentMode"->"SplitRealImag", "Integrator"->"VEGAS",
      "NSamples"->10^7, "VegasEpsRel"->1.*^-9,
      "RunChecks"->False, "Verbose"->False, "WorkingDirectory"->$cxWorkDir];
    If[AssociationQ[resV] && KeyExistsQ[resV, "Results"],
      vRe = resV["Results"][[1]]["Re"]; vIm = resV["Results"][[1]]["Im"];
      relErrV = Abs[(vRe + I vIm) - ref]/Abs[ref];
      Print["  (3) VEGAS 1e7 = ", vRe, " + ", vIm, " I   relErr = ", N[relErrV]];
      If[NumericQ[relErrV] && relErrV < 5*^-2,
        Print["  (3) VEGAS value: PASS (relErr < 5e-2 vs complex P^B; magnitude factor restored)"],
        Print["  (3) VEGAS value: FAIL (relErr=", N[relErrV], ")"]; allPass = False],
      Print["  (3) VEGAS run FAILED: ", resV]; allPass = False]
  ];

  Print[]; Print[If[allPass, "24G PASS", "24G FAIL"]];
  allPass
];


(* ============================================================================
   Test 24H: BUG 2 (TMCv2_BUG_LOG.md) — auxiliary variable on RAW-SYMBOL bases.

   LiftCoefficients built the aux variable as Head[vars[[1]]][n+1].  For INDEXED
   variables x[1],x[2] the head is x, so x[n+1] is fine.  For RAW symbols (e.g.
   ba1, bb1 — as the CALC2 bubble spec used) Head[ba1] is the literal `Symbol`,
   so it built Symbol[4] — an invalid, NON-atomic expression that trips
   `Symbol::string: String expected at position 1 in Symbol[4].` on every lift.
   The fix mints a guaranteed-fresh inert Unique[] symbol for raw-symbol bases
   (indexed bases keep x[n+1]).

   PASS gates:
     (1) LiftCoefficients does NOT emit Symbol::string.
     (2) AuxVariable is a genuine fresh ATOMIC Symbol (AtomQ True; Symbol[4] is
         non-atomic, AtomQ False), distinct from the base variables.
     (3) the lift identity (lifted /. aux -> z0 == original) still holds.
   ============================================================================ *)

RunTestCxH[] := Module[
  {vars, poly, spec, lr, firedSymStr, lc, ld, av, ls, idOK, allPass},

  Print["--- Test 24H: BUG 2 — raw-symbol base vars, no Symbol[k] aux ---"];
  allPass = True;

  ClearAll[ba1, bb1];
  vars = {ba1, bb1};                                   (* RAW symbols, not x[i] *)
  poly = 1 + 10^-4 ba1^2 + bb1^2 + ba1 bb1^2;
  spec = <|"Polynomials"->{poly}, "MonomialExponents"->{0,0},
           "PolynomialExponents"->{-2}, "Variables"->vars,
           "KinematicSymbols"->{}, "RegulatorSymbol"->None|>;
  lr   = {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>};

  (* (1) LiftCoefficients must NOT emit Symbol::string (the Symbol[4] symptom).
     Check returns the failexpr (True) iff that message fires; otherwise lc is
     the Association. *)
  firedSymStr = False;
  lc = Check[LiftCoefficients[spec, lr], firedSymStr = True, Symbol::string];
  If[!firedSymStr,
    Print["  (1) no Symbol::string warning: PASS"],
    Print["  (1) Symbol::string fired (built Symbol[4]): FAIL"]; allPass = False
  ];
  If[!AssociationQ[lc],
    Print["  LiftCoefficients did not return an Association: FAIL"];
    Print[]; Print["24H FAIL"]; Return[False]
  ];

  (* (2) aux variable is a fresh atomic Symbol, distinct from base vars *)
  ld = lc["LiftData"]; av = ld["AuxVariable"];
  If[AtomQ[av] && Head[av] === Symbol && !MemberQ[vars, av],
    Print["  (2) aux var is a fresh atomic Symbol (", av, "): PASS"],
    Print["  (2) aux var malformed: FAIL -> ", av,
          " (AtomQ=", AtomQ[av], ", Head=", Head[av], ")"]; allPass = False
  ];

  (* (3) lift identity still holds with the fresh aux symbol *)
  ls = lc["LiftedSpec"];
  idOK = Simplify[(ls["Polynomials"][[1]] /. av -> ld["z0"]) - poly] === 0;
  If[idOK,
    Print["  (3) identity (lifted /. aux->z0 == original): PASS"],
    Print["  (3) identity: FAIL"]; allPass = False
  ];

  Print[]; Print[If[allPass, "24H PASS", "24H FAIL"]];
  allPass
];


(* ============================================================================
   Umbrella
   ============================================================================ *)

RunTestCx[] := Module[{pa, pb, pc, pd, pe, pf, pg, ph, all},
  Print[""];
  Print["================================================================"];
  Print["  Test 24: lifting with COMPLEX small-magnitude coefficients"];
  Print["================================================================"];
  Print[""];

  pa = RunTestCxA[]; Print[""];
  pb = RunTestCxB[]; Print[""];
  pc = RunTestCxC[]; Print[""];
  pd = RunTestCxD[]; Print[""];
  pe = RunTestCxE[]; Print[""];
  pf = RunTestCxF[]; Print[""];
  pg = RunTestCxG[]; Print[""];
  ph = RunTestCxH[]; Print[""];

  all = pa && pb && pc && pd && pe && pf && pg && ph;

  Print["================================================================"];
  Print["  Test 24 Summary"];
  Print["  24A (single complex coeff, full C++):      ", If[pa, "PASS", "FAIL"]];
  Print["  24B (very small, automatic k*=3):          ", If[pb, "PASS", "FAIL"]];
  Print["  24C (multi-rule complex residuals):        ", If[pc, "PASS", "FAIL"]];
  Print["  24D (complex exponent -> liftcomplex):     ", If[pd, "PASS", "FAIL"]];
  Print["  24E (Fix 1: real float, no false positive):", If[pe, "PASS", "FAIL"]];
  Print["  24F (Fix 2+3: SplitRealImag mode, C++ MC): ", If[pf, "PASS", "FAIL"]];
  Print["  24G (BUG 1: complex base, log P vs log|P|):", If[pg, "PASS", "FAIL"]];
  Print["  24H (BUG 2: raw-symbol aux, no Symbol[k]):  ", If[ph, "PASS", "FAIL"]];
  Print["  Overall Test 24: ", If[all, "PASS", "FAIL"]];
  Print["================================================================"];
  Print[""];
  all
];


(* --- Execute --- *)
If[RunTestCx[], Print["ALL COMPLEX-LIFT TESTS PASSED"],
   Print["SOME COMPLEX-LIFT TESTS FAILED"]; Exit[1]];
