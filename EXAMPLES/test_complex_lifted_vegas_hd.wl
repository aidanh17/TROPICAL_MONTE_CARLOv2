(* ============================================================================
   test_complex_lifted_vegas_hd.wl — Test 26: an 8-DIMENSIONAL generalized Euler
   integral with a TINY COMPLEX polynomial coefficient, evaluated with
   CUBA-Vegas, BOTH WITHOUT and WITH auxiliary-variable lifting.

   This is the high-dimensional, complex-coefficient, VEGAS stress test the
   existing suite lacked (Tests 24/25 are 2-4 variables; test_vegas.wl is n<=4).

   ----------------------------------------------------------------------------
   Integrand (n = 8):

     I = Int_{[0,inf)^8} ( 1 + c*x1^2 + x1 + x2^2 + x3^2 + ... + x8^2 )^(-6) dx
         c = (1 + I) * 1e-6           (tiny: |c| = sqrt(2)*1e-6 ~ 1.41e-6)

   Two structural choices make this a clean, decisive test:

   (a) x1 appears in TWO monomials (x1 and x1^2).  Lifting the c*x1^2 coefficient
       (z^2 * x1^2, k=2) therefore keeps the (n+1)=9-dimensional lifted Newton
       polytope FULL-dimensional (a single extra monomial that already shares its
       support would leave it degenerate).  So BOTH the 8D unlifted and the 9D
       lifted fans exist and are simplicial.

   (b) The integrand factorizes (P is a sum over disjoint variable blocks), so an
       EXACT high-precision reference is available from a single 1D NIntegrate:
         I = (sqrt(pi)/2)^7 * Gamma(5/2)/Gamma(6) * Int_0^inf (1+t+c t^2)^(-5/2) dt
       (integrate the seven x_i^2, i>=2, analytically; the x1 piece by 1D NI).
       The reference METHOD is itself validated against full NIntegrate at n=2,3.

   The tiny coefficient makes the UNLIFTED integrand heavy-tailed (the tropically
   dominant x1^2 vertex carries a ~1e-6 coefficient), which is exactly the
   variance pathology lifting cures -- and where adaptive VEGAS earns its keep.

   PASS gates:
     (0) reference method matches full NIntegrate at n=2,3 (< 1e-6).
     (1) DetectExtremeCoefficients flags c with SuggestedK = 2; |residual| = 1.
     (2) UNLIFTED VEGAS reproduces the reference (relErr < 1e-2).
     (3) LIFTED   VEGAS reproduces the reference (relErr < 5e-3) -- complex value.
     (4) lifted == unlifted to combined accuracy (exactness of lifting).
     (5) [reported] lifting tightens VEGAS here (lifted relErr < unlifted); the
         exact ordering depends on the maxeval/NStart split, so it is informational.

   Notes / dependencies (see the two HIGH-DIMENSION fixes this test exercises):
     * The fans below were computed once via the package's scale-robust path and
       are embedded for speed and determinism.  `FanData -> Automatic` ALSO works
       now: the automatic path retries on a scaled polytope (computeFanScaled),
       fixing the n>=4 failure where a raw ComputeDecomposition leaks $Failed.
     * VEGAS sizing is left Automatic: in dimension 8 the package now scales the
       per-iteration grid up (NStart ~ 8e4) instead of the low-dim default 1000,
       which otherwise makes VEGAS converge to a confidently wrong value.

   Requires: tropical_eval.wl; g++; CUBA (Integrator -> "VEGAS").
   Run:  cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file EXAMPLES/test_complex_lifted_vegas_hd.wl
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
Print["tropical_eval.wl loaded"];
$wd = FileNameJoin[{Directory[], "INTERFILES"}];

(* --- CUBA presence (Integrator -> "VEGAS" needs cuba.h + libcuba) --- *)
cubaAvailable = AnyTrue[{"/opt/homebrew", "/usr/local", "/usr"},
  FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
  (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &];
If[! cubaAvailable,
  Print["test_complex_lifted_vegas_hd.wl requires CUBA for Integrator -> \"VEGAS\"; \
install via `brew install cuba`. Exiting (skip)."];
  Exit[0]];

(* ============================================================================
   Embedded fans (computed once via the scale-robust normal-fan path; the normal
   fan is scale-invariant).  Unlifted: 8D, 9 simplicial cones.  Lifted: 9D, 10
   simplicial cones (full-dimensional => non-degenerate lift).
   ============================================================================ *)
fanU = {{{1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12}, {0, 0, 0, 0, 0, 0, 0, -1},
   {-1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, -1, 0}, {0, -1, 0, 0, 0, 0, 0, 0},
   {0, 0, -1, 0, 0, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0, 0, 0},
   {0, 0, 0, 0, 0, -1, 0, 0}},
  {{2, 3, 4, 5, 6, 7, 8, 9}, {1, 2, 4, 5, 6, 7, 8, 9}, {1, 2, 3, 4, 6, 7, 8, 9},
   {1, 2, 3, 4, 5, 7, 8, 9}, {1, 2, 3, 4, 5, 6, 8, 9}, {1, 2, 3, 4, 5, 6, 7, 9},
   {1, 2, 3, 4, 5, 6, 7, 8}, {1, 2, 3, 5, 6, 7, 8, 9}, {1, 3, 4, 5, 6, 7, 8, 9}}};

fanL = {{{1/6, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, -1/12}, {0, 0, 0, 0, 0, 0, 0, -1, 0},
   {0, 0, 0, 0, 0, 0, 0, 0, -1}, {-1, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, -1, 0, 0},
   {0, -1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 0, 0},
   {0, 0, 0, 0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -1, 0, 0, 0}},
  {{2, 3, 4, 5, 6, 7, 8, 9, 10}, {1, 2, 3, 5, 6, 7, 8, 9, 10}, {1, 2, 3, 4, 5, 7, 8, 9, 10},
   {1, 2, 3, 4, 5, 6, 8, 9, 10}, {1, 2, 3, 4, 5, 6, 7, 9, 10}, {1, 2, 3, 4, 5, 6, 7, 8, 10},
   {1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 2, 3, 4, 6, 7, 8, 9, 10}, {1, 3, 4, 5, 6, 7, 8, 9, 10},
   {1, 2, 4, 5, 6, 7, 8, 9, 10}}};

RunTest26[] := Module[
  {c, vars, P, spec, s, ref, refRe, refIm, lr, det,
   relMethod, fanOK, residual,
   rUV, rLV, rUM, rLM,
   reU, reL, relU, relL, allPass},

  Print["\n================================================================"];
  Print["  Test 26: 8D integral, tiny complex coeff, VEGAS, +/- lifting"];
  Print["================================================================\n"];
  allPass = True;

  c    = (1 + I) 10^-6;
  vars = Table[x[i], {i, 8}];
  P    = 1 + c x[1]^2 + x[1] + Sum[x[i]^2, {i, 2, 8}];
  s    = 6;
  spec = <|"Polynomials" -> {P}, "MonomialExponents" -> ConstantArray[0, 8],
           "PolynomialExponents" -> {-s}, "Variables" -> vars,
           "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  lr   = {<|"PolyIndex" -> 1, "ExponentVector" -> Prepend[ConstantArray[0, 7], 2],
            "k" -> 2|>};

  (* ---- reference: factorize 7 quadratics analytically, x1 by 1D NIntegrate ---- *)
  ref = (Sqrt[Pi]/2)^7 Gamma[s - 7/2]/Gamma[s] *
    NIntegrate[(1 + t + c t^2)^(7/2 - s), {t, 0, Infinity},
      MaxRecursion -> 60, PrecisionGoal -> 14, WorkingPrecision -> 40];
  refRe = Re[ref]; refIm = Im[ref];
  Print["reference I = ", ScientificForm[N[ref, 12]]];

  (* ---- (0) validate the reference METHOD vs full NIntegrate at n=2,3 ---- *)
  relMethod = Max @ Table[
    Module[{full, meth, tv = Table[Unique[], {nn}]},
      full = NIntegrate[
        (1 + tv[[1]] + c tv[[1]]^2 + Sum[tv[[i]]^2, {i, 2, nn}])^(-s),
        Evaluate[Sequence @@ ({#, 0, Infinity} & /@ tv)],
        MaxRecursion -> 30, PrecisionGoal -> 8];
      meth = (Sqrt[Pi]/2)^(nn - 1) Gamma[s - (nn - 1)/2]/Gamma[s] *
        NIntegrate[(1 + t + c t^2)^((nn - 1)/2 - s), {t, 0, Infinity},
          MaxRecursion -> 60, PrecisionGoal -> 12, WorkingPrecision -> 30];
      Abs[full - meth]/Abs[full]],
    {nn, {2, 3}}];
  If[relMethod < 1*^-6,
    Print["  (0) reference method validated vs full NIntegrate (n=2,3): PASS (max reldiff ",
          ScientificForm[N[relMethod], 3], ")"],
    Print["  (0) reference method: FAIL (reldiff ", relMethod, ")"]; allPass = False];

  (* ---- (1) detection + unit residual ---- *)
  det = DetectExtremeCoefficients[spec, 1000];
  residual = LiftCoefficients[spec, lr]["LiftData"]["Residuals"][[1]];
  If[Length[det] == 1 && det[[1]]["SuggestedK"] == 2 &&
     Abs[N[Abs[residual]] - 1] < 1*^-9,
    Print["  (1) detect+lift: PASS (k*=2, |C|=", ScientificForm[det[[1]]["Magnitude"], 4],
          ", |residual|=1)"],
    Print["  (1) detect+lift: FAIL -> ", det]; allPass = False];

  (* sanity on the embedded fans *)
  fanOK = (Length[fanU[[1]]] == 9 && AllTrue[fanU[[2]], Length[#] == 8 &] &&
           Length[fanL[[1]]] == 10 && AllTrue[fanL[[2]], Length[#] == 9 &]);
  Print["  fans: unlifted ", Length[fanU[[2]]], " cones (8D), lifted ",
        Length[fanL[[2]]], " cones (9D) -> ", If[fanOK, "OK", "WRONG SHAPE"]];
  If[! fanOK, allPass = False];

  rel[r_] := Abs[(r["Re"] + I r["Im"]) - ref]/Abs[ref];

  (* ---- (2) UNLIFTED VEGAS (continuous sectors; VEGAS's best case) ---- *)
  Print["\n  Running UNLIFTED VEGAS (8D, 9 sectors, maxeval 5e6) ..."];
  rUV = EvaluateTropicalMC[spec, fanU, {{}}, "Integrator" -> "VEGAS",
    "NSamples" -> 5000000, "VegasEpsRel" -> 1.*^-9,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> $wd];
  If[AssociationQ[rUV],
    reU = rUV["Results"][[1]]["Re"]; relU = N[rel[rUV["Results"][[1]]]];
    Print["  (2) unlifted VEGAS: Re=", ScientificForm[reU, 10], " Im=",
          ScientificForm[rUV["Results"][[1]]["Im"], 5], "  relErr=", ScientificForm[relU, 4]];
    If[relU < 1*^-2, Print["  (2) unlifted vs reference: PASS (relErr < 1e-2)"],
      Print["  (2) unlifted vs reference: FAIL"]; allPass = False],
    Print["  (2) unlifted VEGAS FAILED: ", rUV]; allPass = False; relU = Infinity];

  (* ---- (3) LIFTED VEGAS (the headline: complex tiny coeff + lifting + VEGAS) ---- *)
  Print["\n  Running LIFTED VEGAS (9D fan -> 10 sectors, maxeval 5e6) ..."];
  rLV = EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> lr, "FanData" -> fanL,
    "Integrator" -> "VEGAS", "NSamples" -> 5000000, "VegasEpsRel" -> 1.*^-9,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> $wd];
  If[AssociationQ[rLV],
    reL = rLV["Results"][[1]]["Re"]; relL = N[rel[rLV["Results"][[1]]]];
    Print["  (3) lifted VEGAS:   Re=", ScientificForm[reL, 10], " Im=",
          ScientificForm[rLV["Results"][[1]]["Im"], 5], "  relErr=", ScientificForm[relL, 4]];
    Print["      reference:      Re=", ScientificForm[refRe, 10], " Im=", ScientificForm[refIm, 5]];
    If[relL < 5*^-3, Print["  (3) lifted vs reference: PASS (relErr < 5e-3, complex value)"],
      Print["  (3) lifted vs reference: FAIL"]; allPass = False],
    Print["  (3) lifted VEGAS FAILED: ", rLV]; allPass = False; relL = Infinity];

  (* ---- (4) lifted == unlifted (exactness of the lifting reformulation) ---- *)
  If[NumericQ[reU] && NumericQ[reL],
    If[Abs[reL - reU] < 1*^-2 Abs[ref],
      Print["\n  (4) lifted == unlifted: PASS (|dRe|=", ScientificForm[N[Abs[reL - reU]], 3],
            " < 1% of |I|)"],
      Print["\n  (4) lifted == unlifted: FAIL (|dRe|=", ScientificForm[N[Abs[reL - reU]], 3], ")"];
      allPass = False]];

  (* ---- (5) lifting tightens VEGAS on this heavy-tailed integrand ---- *)
  If[NumericQ[relU] && NumericQ[relL],
    Print["  (5) accuracy: unlifted relErr=", ScientificForm[relU, 3],
          ", lifted relErr=", ScientificForm[relL, 3],
          " -> lifting ", If[relL <= relU, "TIGHTER (PASS)", "not tighter (INFO)"]];
    If[relL > relU, Print["      (note: both already meet their gates; ordering is informational)"]]];

  (* ---- contrast: plain uniform MC (heavy-tailed here; error bars optimistic) ---- *)
  Print["\n  [contrast] plain uniform MC (2e6 samples) -- expected noisier/biased:"];
  rUM = EvaluateTropicalMC[spec, fanU, {{}}, "Integrator" -> "MC", "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> $wd];
  rLM = EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> lr, "FanData" -> fanL,
    "Integrator" -> "MC", "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> $wd];
  If[AssociationQ[rUM], Print["      unlifted MC: relErr=", ScientificForm[N[rel[rUM["Results"][[1]]]], 3]]];
  If[AssociationQ[rLM], Print["      lifted   MC: relErr=", ScientificForm[N[rel[rLM["Results"][[1]]]], 3]]];

  Print["\n================================================================"];
  Print["  Test 26 overall: ", If[allPass, "PASS", "FAIL"]];
  Print["================================================================\n"];
  allPass
];

If[RunTest26[], Print["TEST 26 PASSED"], Print["TEST 26 FAILED"]; Exit[1]];
