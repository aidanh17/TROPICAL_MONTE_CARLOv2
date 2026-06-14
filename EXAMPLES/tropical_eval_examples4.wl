(* ============================================================================
   tropical_eval_examples4.wl

   Examples 24-27: auxiliary-variable LIFTING with COMPLEX, small-magnitude
   polynomial coefficients.

   Examples 2/15-19 made the integrand complex through the EXPONENTS A_i, B_j.
   Example 20 lifted a small but REAL coefficient.  This file combines the two
   regimes the package had not previously exercised together: a polynomial
   coefficient that is itself COMPLEX and of very small magnitude, fed through
   the lifting module.

   Why this works.  Lifting replaces an extreme coefficient C of monomial
   C x^alpha by z^k with anchor z0 = |C|^{1/k} (real, by construction).  When C
   is complex the MAGNITUDE goes into the real anchor z0, and the PHASE is
   carried exactly by the residual c = C/z0^k -- an O(1) complex number that
   sits on the lifted monomial.  Everything that must stay real to define the
   sector geometry (effective exponents, the convergence gate, the delta-root
   domain indicator) depends only on the integer Newton-polytope data and the
   real exponents A, B, so it is untouched; only the polynomial values become
   complex, and the C++ pipeline already carries those as std::complex<double>.
   The round-trip identity (lifted poly at z->z0 equals the original) therefore
   holds exactly, and the lifted decomposition reproduces the complex integral.

   Contrast with complex exponents.  Lifting still requires REAL polynomial
   exponents B_j: a complex B makes the effective exponent atilde complex, which
   admits no real-valued domain indicator, so ProcessSectorLifted fires
   TropicalEval::liftcomplex (shown at the end of Example 27).  Complex
   COEFFICIENTS are fine; complex EXPONENTS are not.

   All integrals here are convergent and cross-checked against direct
   NIntegrate.  Example 24 additionally runs the full C++ Monte Carlo pipeline
   (requires g++); 25-27 use the quadrature validator.

   Load the package first:
     SetDirectory[FileNameJoin[{NotebookDirectory[], ".."}]];
     Get["tropical_eval.wl"];
   Or run as a script:
     wolframscript -file EXAMPLES/tropical_eval_examples4.wl

   See EXAMPLES/test_complex_lifted.wl for the same integrands as automated
   PASS/FAIL assertions.
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Print["tropical_eval.wl loaded successfully"];
Print[];


(* ============================================================================
   Example 24: a single small COMPLEX coefficient, end to end

   Integral[0,Inf] dx1 dx2 / (1 + (1+I) 10^-4 x1^2 + x2^2 + x1 x2^2)^2

   Same Newton polytope as the real Case A / Example 20, but the x1^2
   coefficient is complex, |C| = Sqrt[2] 10^-4 ~ 1.41e-4.  The exact integral
   is a complex number near 11.63 - 2.58 I.
   ============================================================================ *)

Print["=== Example 24: lifting a single small complex coefficient ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + (1+I)10^-4 x1^2 + x2^2 + x1 x2^2)^2"];
Print[];

Module[{c, poly, vars, spec, extreme, liftRules, lc, liftedSpec, liftData,
        lverts, liftedFan, vl, ref, resPlain, resLifted, workDir},

  c    = (1 + I) 10^-4;
  poly = 1 + c x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};
  workDir = FileNameJoin[{Directory[], "INTERFILES"}];

  spec = <|
    "Polynomials"         -> {poly},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  (* --- Step 1: detect the extreme coefficient (magnitude via |C|) --- *)
  extreme = DetectExtremeCoefficients[spec, 1000];
  Print["DetectExtremeCoefficients: ", extreme];
  Print["  -> magnitude is |C| (the complex modulus); kStar auto-suggests k = ",
        First[#["SuggestedK"] & /@ extreme]];
  Print[];

  (* --- Step 2: lift with k = 2.  Magnitude -> anchor, phase -> residual --- *)
  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};
  lc = LiftCoefficients[spec, liftRules];
  liftedSpec = lc["LiftedSpec"];
  liftData   = lc["LiftData"];
  Print["Anchor z0   = ", liftData["z0"], "  (= |C|^{1/2}, real) ~ ",
        N[liftData["z0"]]];
  Print["Residual c  = ", liftData["Residuals"][[1]], " = (1+I)/Sqrt[2]"];
  Print["  |residual| = ", N[Abs[liftData["Residuals"][[1]]]],
        "  (unit modulus: the phase, now an O(1) coefficient)"];
  Print["Lifted poly = ", liftedSpec["Polynomials"][[1]]];
  Print["Identity (lifted poly at z -> z0 equals original): ",
        Simplify[(liftedSpec["Polynomials"][[1]] /.
                  liftData["AuxVariable"] -> liftData["z0"]) - poly] === 0];
  Print[];

  (* --- Step 3: lifted fan (automatic; polytope is full-dimensional) --- *)
  lverts = PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                            liftedSpec["Variables"]];
  liftedFan = ComputeDecomposition[lverts, "ShowProgress" -> False];
  Print["Lifted fan: ", Length[liftedFan[[1]]], " rays, ",
        Length[liftedFan[[2]]], " sectors (3 variables)"];
  Print[];

  (* --- Steps 4-5: exactness of the lifted decomposition (Re AND Im) --- *)
  Print["ValidateLiftedDecomposition (delta resolved per sector):"];
  vl = Quiet @ ValidateLiftedDecomposition[spec, liftedSpec, liftedFan,
                                           liftData, {}, 4];
  Print["  Direct NIntegrate: ", vl["DirectResult"]];
  Print["  Lifted sector sum: ", vl["SectorSum"]];
  Print["  Relative error:    ", vl["RelativeError"], "  <- exact to ~1e-5"];
  Print["  Dropped (EmptyDomain) sectors: ",
        Length[Lookup[vl, "DroppedSectors", {}]]];
  Print[];

  (* --- Step 6: plain vs lifted C++ Monte Carlo, same sample count --- *)
  ref = Quiet @ NIntegrate[1/(1 + c t1^2 + t2^2 + t1 t2^2)^2,
    {t1, 0, Infinity}, {t2, 0, Infinity}, MaxRecursion -> 25, PrecisionGoal -> 8];
  Print["NIntegrate reference = ", ref];
  Print[];

  Print["Plain pipeline (no lifting), 2*10^5 samples:"];
  Module[{verts, fanData},
    verts   = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1), vars];
    fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
    resPlain = EvaluateTropicalMC[spec, fanData, {{}},
      "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False,
      "WorkingDirectory" -> workDir];
  ];
  If[AssociationQ[resPlain],
    Print["  plain MC  = ", resPlain["Results"][[1]]["Re"], " + ",
          resPlain["Results"][[1]]["Im"], " I"];
    Print["              +/- (", resPlain["Results"][[1]]["ReErr"], ", ",
          resPlain["Results"][[1]]["ImErr"], ")"];
    Print["    deviation from ref: (", resPlain["Results"][[1]]["Re"] - Re[ref],
          ", ", resPlain["Results"][[1]]["Im"] - Im[ref], ")  <- heavy tail"];,
    Print["  plain pipeline failed (is g++ installed?)"]
  ];
  Print[];

  Print["Lifted pipeline (k=2), same 2*10^5 samples:"];
  resLifted = EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> liftRules,
    "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> workDir];
  If[AssociationQ[resLifted] && KeyExistsQ[resLifted, "Results"],
    Print["  lifted MC = ", resLifted["Results"][[1]]["Re"], " + ",
          resLifted["Results"][[1]]["Im"], " I"];
    Print["              +/- (", resLifted["Results"][[1]]["ReErr"], ", ",
          resLifted["Results"][[1]]["ImErr"], ")"];
    Print["    deviation from ref: (", resLifted["Results"][[1]]["Re"] - Re[ref],
          ", ", resLifted["Results"][[1]]["Im"] - Im[ref], ")"];
    Print["  -> both Re and Im consistent with the reference; the magnitude"];
    Print["     scale sits in the deterministic prefactor, not the random part."];,
    Print["  lifted pipeline failed (is g++ installed?)"]
  ];
];
Print[];


(* ============================================================================
   Example 25: a VERY small complex coefficient, automatic anchor

   c = (1+I) 10^-7  =>  |C| ~ 1.41e-7.  The kStar rule
   k* = Ceiling[|log10|C|| / 3] now returns 3 (not 2), so "LiftRules" ->
   Automatic lifts at k = 3, z0 = |C|^{1/3} ~ 5.2e-3.  The decomposition stays
   exact.  This is the "magnitudes are very small" regime where the anchor
   integer grows with |log|C||.
   ============================================================================ *)

Print["=== Example 25: very small complex coefficient, automatic k* ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + (1+I)10^-7 x1^2 + x2^2 + x1 x2^2)^2"];
Print[];

Module[{c, poly, vars, spec, det, k, liftRules, lc, liftedSpec, liftData,
        lverts, liftedFan, vl},

  c    = (1 + I) 10^-7;
  poly = 1 + c x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};
  spec = <|
    "Polynomials"         -> {poly},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  det = DetectExtremeCoefficients[spec, 1000];
  k   = det[[1]]["SuggestedK"];
  Print["DetectExtremeCoefficients: |C| = ", det[[1]]["Magnitude"],
        ", SuggestedK (kStar) = ", k, "  (grew to 3 as |C| shrank)"];

  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> k|>};
  lc = LiftCoefficients[spec, liftRules];
  liftedSpec = lc["LiftedSpec"]; liftData = lc["LiftData"];
  Print["Anchor z0 = ", N[liftData["z0"]], " = |C|^{1/3};  residual = ",
        N[liftData["Residuals"][[1]]], " (still unit modulus)"];

  lverts = PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                            liftedSpec["Variables"]];
  liftedFan = ComputeDecomposition[lverts, "ShowProgress" -> False];
  Print["Lifted fan: ", Length[liftedFan[[1]]], " rays, ",
        Length[liftedFan[[2]]], " sectors"];

  vl = Quiet @ ValidateLiftedDecomposition[spec, liftedSpec, liftedFan,
                                           liftData, {}, 4];
  Print["ValidateLiftedDecomposition:"];
  Print["  Direct NIntegrate: ", vl["DirectResult"]];
  Print["  Lifted sector sum: ", vl["SectorSum"]];
  Print["  Relative error:    ", vl["RelativeError"], "  <- exact"];
];
Print[];


(* ============================================================================
   Example 26: several small complex coefficients (multi-rule lift)

   P = 1 + (1+I)10^-4 x1^2 + x2^2 + (1-2I)10^-4 x1 x2^2,  B = -2.

   Both extreme monomials are lifted against ONE shared auxiliary variable.
   The anchor z0 is fixed from the primary (most extreme) coefficient; each
   residual c_i = C_i/z0^{k_i} keeps its own O(1) phase, so the round-trip
   stays exact for several lifted monomials at once.
   ============================================================================ *)

Print["=== Example 26: several small complex coefficients (multi-rule) ==="];
Print["P = 1 + (1+I)10^-4 x1^2 + x2^2 + (1-2I)10^-4 x1 x2^2,  B = -2"];
Print[];

Module[{c1, c2, poly, vars, spec, liftRules, lc, liftedSpec, liftData,
        lverts, liftedFan, vl},

  c1 = (1 + I) 10^-4;  c2 = (1 - 2 I) 10^-4;
  poly = 1 + c1 x[1]^2 + x[2]^2 + c2 x[1] x[2]^2;
  vars = {x[1], x[2]};
  spec = <|
    "Polynomials"         -> {poly},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>,
               <|"PolyIndex" -> 1, "ExponentVector" -> {1, 2}, "k" -> 2|>};
  lc = LiftCoefficients[spec, liftRules];
  liftedSpec = lc["LiftedSpec"]; liftData = lc["LiftData"];

  Print["Shared anchor z0 = ", N[liftData["z0"]]];
  Print["Residuals (one per lifted monomial) = ", N[liftData["Residuals"]]];
  Print["  magnitudes = ", N[Abs /@ liftData["Residuals"]], " (both O(1))"];
  Print["Lifted poly = ", liftedSpec["Polynomials"][[1]]];
  Print["Identity (lifted poly at z -> z0 equals original): ",
        Simplify[(liftedSpec["Polynomials"][[1]] /.
                  liftData["AuxVariable"] -> liftData["z0"]) - poly] === 0];
  Print[];

  lverts = PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                            liftedSpec["Variables"]];
  liftedFan = ComputeDecomposition[lverts, "ShowProgress" -> False];
  Print["Lifted fan: ", Length[liftedFan[[1]]], " rays, ",
        Length[liftedFan[[2]]], " sectors"];

  vl = Quiet @ ValidateLiftedDecomposition[spec, liftedSpec, liftedFan,
                                           liftData, {}, 4];
  Print["ValidateLiftedDecomposition:"];
  Print["  Direct NIntegrate: ", vl["DirectResult"]];
  Print["  Lifted sector sum: ", vl["SectorSum"]];
  Print["  Relative error:    ", vl["RelativeError"], "  <- exact"];
];
Print[];


(* ============================================================================
   Example 27: a purely imaginary small coefficient, and the boundary case

   Part A: c = I 10^-4 (pure phase).  |C| = 10^-4, k* = 2, z0 = 10^-2,
           residual = I.  Lifting and the decomposition are exact.

   Part B: the supported envelope.  Complex coefficients are fine; complex
           polynomial EXPONENTS are not.  With B = -(2 + I/2) the effective
           exponents become complex and ProcessSectorLifted fires
           TropicalEval::liftcomplex (the message below is expected).
   ============================================================================ *)

Print["=== Example 27: purely imaginary coefficient + the complex-exponent boundary ==="];
Print[];

Print["Part A: c = I 10^-4 in  1/(1 + I 10^-4 x1^2 + x2^2 + x1 x2^2)^2"];
Module[{c, poly, vars, spec, liftRules, lc, liftedSpec, liftData,
        lverts, liftedFan, vl},
  c    = I 10^-4;
  poly = 1 + c x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
           "PolynomialExponents" -> {-2}, "Variables" -> vars,
           "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};
  lc = LiftCoefficients[spec, liftRules];
  liftedSpec = lc["LiftedSpec"]; liftData = lc["LiftData"];
  Print["  z0 = ", N[liftData["z0"]], ", residual = ", liftData["Residuals"][[1]],
        " (= I; |residual| = 1)"];
  lverts = PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                            liftedSpec["Variables"]];
  liftedFan = ComputeDecomposition[lverts, "ShowProgress" -> False];
  vl = Quiet @ ValidateLiftedDecomposition[spec, liftedSpec, liftedFan,
                                           liftData, {}, 4];
  Print["  Direct NIntegrate: ", vl["DirectResult"]];
  Print["  Lifted sector sum: ", vl["SectorSum"]];
  Print["  Relative error:    ", vl["RelativeError"], "  <- exact"];
];
Print[];

Print["Part B: complex EXPONENT B = -(2+I/2) is NOT liftable (expected message)"];
Module[{poly, vars, spec, liftRules, lc, liftedSpec, liftData,
        lverts, liftedFan, sd, fired},
  poly = 1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
           "PolynomialExponents" -> {-(2 + I/2)}, "Variables" -> vars,
           "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};
  lc = LiftCoefficients[spec, liftRules];
  liftedSpec = lc["LiftedSpec"]; liftData = lc["LiftData"];
  lverts = PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                            liftedSpec["Variables"]];
  liftedFan = ComputeDecomposition[lverts, "ShowProgress" -> False];
  fired = False;
  Check[
    Do[sd = ProcessSectorLifted[liftedSpec, liftedFan[[1]], liftedFan[[2, s]],
                                s, liftData];
       If[sd === $Failed, fired = True], {s, Length[liftedFan[[2]]]}],
    fired = True, TropicalEval::liftcomplex];
  Print["  liftcomplex fired (complex exponents unsupported, as documented): ",
        fired];
];
Print[];

Print["=== Examples 24-27 complete ==="];
