(* ============================================================================
   tropical_eval_examples2.wl

   Examples 15-20: complex exponents and small/extreme coefficients.

   Focus of this file:
     - polynomials raised to COMPLEX exponents (Examples 15, 16, 17, 19)
     - polynomials with SMALL coefficients   (Examples 18, 19, 20)
   All integrals here are convergent (TROPICAL_MONTE_CARLOv2 supports
   convergent integrals only) and every example is cross-checked against
   direct NIntegrate.

   Requires: tropical_fan.wl, tropical_eval.wl, Polymake.
   Examples 17 and 20 additionally require g++ (they run the full C++
   Monte Carlo pipeline).

   Load the package first:
     SetDirectory[FileNameJoin[{NotebookDirectory[], ".."}]];
     Get["tropical_eval.wl"];

   Or run as a script:
     wolframscript -file EXAMPLES/tropical_eval_examples2.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["tropical_eval.wl loaded successfully"];
Print[];


(* ============================================================================
   Example 15: Two polynomials, BOTH raised to complex exponents

   Integral[0,Inf] dx1 dx2 (1+x1+x2)^{-(3/2 + I/2)} * (1+x1*x2)^{-(1 - I/3)}

   Example 2 showed a single polynomial with a complex exponent; here a
   PRODUCT of two polynomials carries two different complex exponents.
   Convergence is governed by the real parts alone: the fan is computed
   from the product P1*P2 (coefficients and exponents are irrelevant for
   the fan), and every sector must end up with Re(a_eff) > 0.

   The imaginary parts enter the effective exponents, the flattening, and
   the prefactor, making every sector integrand complex-valued.
   ============================================================================ *)

Print["=== Example 15: Two polynomials with complex exponents ==="];
Print["Integral[0,Inf] dx1 dx2 (1+x1+x2)^{-(3/2+I/2)} (1+x1*x2)^{-(1-I/3)}"];
Print[];

Module[{p1, p2, vars, spec, fanPoly, verts, fanData, vr},

  p1   = 1 + x[1] + x[2];
  p2   = 1 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {p1, p2},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-(3/2 + I/2), -(1 - I/3)},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* Fan from the product of the polynomials (first powers suffice) *)
  fanPoly = p1 * p2;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[1]]], " rays, ",
        Length[fanData[[2]]], " sectors"];

  (* All sectors: complex effective exponents, all with Re > 0 *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", N[sd["NewExponents"]],
              ", div = ", sd["IsDivergent"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  (* Spot-check the flattened integrand magnitude in sector 1: even with
     complex exponents it should be O(1) on [0,1]^2 *)
  Print[];
  Print["Flattening magnitude check, sector 1 (should be O(1)):"];
  Module[{sd, check},
    sd    = ProcessSector[spec, fanData[[1]], fanData[[2, 1]], 1];
    check = CheckFlatteningMagnitude[sd, 100];
    Print["  Min |f| = ", check["Min"], "  Mean = ", check["Mean"],
          "  Max = ", check["Max"]];
  ];

  Print[];
  Print["Validation (complex result):"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 16: Complex monomial exponents AND complex polynomial exponent

   Integral[0,Inf] dx1 dx2 x1^{1/2+I/3} x2^{-1/4-I/5}
                            / (1 + x1^2 + x2^2 + x1*x2)^{2+I}

   The monomial prefactor exponents A_i are complex too.  Only their real
   parts matter for convergence:
     - at x_i -> 0:   need Re(A_i) > -1          (here 1/2 and -1/4)
     - at infinity:   (Re(A)+1) must lie in the interior of
                      Re(-B) * Newt(P)            (polytope condition)
   The imaginary parts produce oscillatory factors x^{I*b} = e^{I b log x}
   that the log-exp evaluation handles without branch-cut issues.
   ============================================================================ *)

Print["=== Example 16: Complex monomial + polynomial exponents ==="];
Print["Integral[0,Inf] dx1 dx2 x1^{1/2+I/3} x2^{-1/4-I/5} / (1+x1^2+x2^2+x1*x2)^{2+I}"];
Print[];

Module[{poly, vars, spec, verts, fanData, vr},

  poly = 1 + x[1]^2 + x[2]^2 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {1/2 + I/3, -1/4 - I/5},
    "PolynomialExponents" -> {-(2 + I)},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[1]]], " rays, ",
        Length[fanData[[2]]], " sectors"];

  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", N[sd["NewExponents"]],
              ", prefactor = ", N[sd["Prefactor"]]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation (complex result):"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 17: Full C++ Monte Carlo pipeline with a complex exponent
               (EvaluateTropicalMC driver + kinematic scan)

   Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^{2+I/2}

   Examples 6-7 ran the C++ pipeline step by step for REAL exponents.
   This example uses the top-level driver EvaluateTropicalMC end-to-end
   with a COMPLEX polynomial exponent: one compilation, five values of
   lam, complex MC estimates (Re, Im) with independent error bars,
   cross-checked against NIntegrate at every point.

   Requires g++.  Generated files go to INTERFILES/.
   ============================================================================ *)

Print["=== Example 17: C++ MC with complex exponent (driver) ==="];
Print["Integral[0,Inf] dx1 dx2 / (1+lam*x1^2+x2^2+x1*x2^2)^{2+I/2}"];
Print["5 values of lam in [0.5, 8]"];
Print[];

Module[{lam, poly, vars, A, spec, verts, fanData,
        lamValues, kinPoints, res},

  lam  = Symbol["lam"];
  A    = 2 + I/2;
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  (* Fan from unit-coefficient proxy (fan is independent of lam) *)
  verts   = PolytopeVertices[(poly /. lam -> 1)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  lamValues = {0.5, 1.0, 2.0, 4.0, 8.0};
  kinPoints = List /@ lamValues;

  res = EvaluateTropicalMC[spec, fanData, kinPoints,
    "NSamples"  -> 200000,
    "RunChecks" -> False,
    "Verbose"   -> False,
    "WorkingDirectory" -> FileNameJoin[{Directory[], "INTERFILES"}]];

  If[AssociationQ[res],
    Print[];
    Print["Results (MC vs NIntegrate, complex values):"];
    Do[
      Module[{r, mc, mcErr, ni, dev},
        r     = res["Results"][[i]];
        mc    = r["Re"] + I r["Im"];
        mcErr = Sqrt[r["ReErr"]^2 + r["ImErr"]^2];
        ni = Quiet@NIntegrate[
          (1 + lamValues[[i]] t1^2 + t2^2 + t1 t2^2)^(-2 - I/2),
          {t1, 0, Infinity}, {t2, 0, Infinity},
          MaxRecursion -> 20, PrecisionGoal -> 5];
        dev = Abs[mc - ni] / Abs[ni];
        Print["  lam = ", lamValues[[i]]];
        Print["    MC         = ", mc, " +/- ", mcErr];
        Print["    NIntegrate = ", ni];
        Print["    |rel dev|  = ", dev];
      ],
      {i, Length[lamValues]}
    ];,
    Print["EvaluateTropicalMC failed (is g++ installed?)"];
  ];
];
Print[];


(* ============================================================================
   Example 18: Small polynomial coefficients — exactness and diagnostics

   Integral[0,Inf] dx1 dx2 / (1 + 10^-4 x1^2 + x2^2 + 10^-3 x1*x2^2)^2

   The tropical fan is computed from the SUPPORT of the polynomial only,
   so coefficients never affect the decomposition: the tropical factoring
   identity P = y^d * Q is exact algebra for ANY coefficients.  What small
   coefficients do change:
     - In sectors whose tropically dominant monomial carries the small
       coefficient, the cleared polynomial Q has a SMALL constant term
       (here 10^-4), so Q^{-2} reaches ~10^8 near the y -> 0 corner.
     - The sector decomposition stays exact, but uniform MC sampling on
       that sector acquires a large per-sample variance (heavy tail).

   This example verifies exactness via NIntegrate and uses
   DetectExtremeCoefficients / CheckFlatteningMagnitude to LOCATE the
   high-variance sectors.  Example 20 shows the mitigation (lifting).
   ============================================================================ *)

Print["=== Example 18: Small coefficients (10^-4, 10^-3) ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + 10^-4 x1^2 + x2^2 + 10^-3 x1*x2^2)^2"];
Print[];

Module[{poly, vars, spec, extreme, verts, fanData, vr},

  poly = 1 + 10^-4 x[1]^2 + x[2]^2 + 10^-3 x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* Identify the extreme coefficients (|C| < 1/100 or > 100) *)
  extreme = DetectExtremeCoefficients[spec, 100];
  Print["DetectExtremeCoefficients (threshold 100):"];
  Do[
    Print["  poly ", e["PolyIndex"], ", monomial x^", e["ExponentVector"],
          ": coeff = ", e["Coefficient"], " (|C| = ", e["Magnitude"], ")"],
    {e, extreme}
  ];
  Print[];

  (* Fan from the unit-coefficient proxy: same support, same fan *)
  verts   = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  Print["Fan: ", Length[fanData[[2]]], " sectors (computed from unit-",
        "coefficient proxy; fan is coefficient-independent)"];
  Print[];

  (* Per-sector magnitude spot check: the spread of Max |f| across sectors
     flags where uniform MC will have high variance *)
  Print["Per-sector flattened-integrand magnitudes (500 random points):"];
  Do[
    Module[{sd, check},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        check = CheckFlatteningMagnitude[sd, 500];
        Print["  Sector ", s, ": Min = ", check["Min"],
              "  Mean = ", check["Mean"], "  Max = ", check["Max"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];
  Print[];

  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 4];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 19: Small coefficients AND a complex exponent

   Integral[0,Inf] dx1 dx2 / (1 + 10^-3 x1^2 + x2^2 + 10^-2 x1*x2^2)^{2+I/2}

   Combines both features of this file: coefficients spanning three orders
   of magnitude and a complex polynomial exponent.  The complex power of a
   positive polynomial is taken as P^{2+I/2} = exp((2+I/2) log P), so the
   integrand is well-defined; convergence is again controlled by Re only.

   Practical note: as coefficients shrink further (here at ~10^-4 and
   below with complex exponents), the sharply peaked oscillatory sector
   integrands start to defeat ADAPTIVE QUADRATURE on both sides of this
   comparison, well before the tropical decomposition itself loses
   anything.  In that regime the C++ Monte Carlo pipeline plus the
   lifting of Example 20 (real exponents) are the robust tools.
   ============================================================================ *)

Print["=== Example 19: Small coefficients + complex exponent ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + 10^-3 x1^2 + x2^2 + 10^-2 x1*x2^2)^{2+I/2}"];
Print[];

Module[{poly, vars, spec, verts, fanData, vr},

  poly = 1 + 10^-3 x[1]^2 + x[2]^2 + 10^-2 x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-(2 + I/2)},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", N[sd["NewExponents"]],
              ", div = ", sd["IsDivergent"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation (complex result):"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 4];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 20: Lifting a small coefficient (variance reduction)

   Integral[0,Inf] dx1 dx2 / (1 + 10^-4 x1^2 + x2^2 + x1*x2^2)^2

   Example 18 showed that small coefficients leave the decomposition exact
   but create heavy-tailed sector integrands.  This example demonstrates
   the auxiliary-variable LIFTING workflow that mitigates it:

     1. DetectExtremeCoefficients finds the 10^-4 coefficient and, via the
        kStar anchor rule (k* = Ceiling[|log|C||/log tau]; see the manual
        section "Choosing the Lift Anchor z0"), auto-suggests k = 2 -- exactly
        the value lifted by hand below.  Step 7 then runs the fully automatic
        path "LiftRules" -> Automatic, which uses that suggestion.
     2. LiftCoefficients replaces 10^-4 x1^2 -> z^2 x1^2 with the anchor
        z0 = (10^-4)^{1/2} = 1/100, exactly:
            I = Int dz delta(z - z0) Int dx [lifted integrand]
     3. The lifted polynomial 1 + x2^2 + x1 x2^2 + x1^2 z^2 has a
        full-dimensional Newton polytope in 3 variables, so the lifted
        fan is computed automatically.
     4. ProcessSectorLifted resolves the delta constraint analytically in
        each 3D sector (pivot substitution), returning 2D sectors with a
        domain constraint.
     5. ValidateLiftedDecomposition confirms the lifted sector sum equals
        the original integral.
     6. EvaluateTropicalMCLifted runs the C++ MC on the lifted sectors;
        compare against the PLAIN pipeline at the same sample count.

   At the time of writing, with 2*10^5 samples (fixed seeds):
        NIntegrate reference:  13.0901
        plain  MC:   7.83 +/- 1.56   <- heavy tail badly undersampled;
                                        the sample error bar UNDERSTATES
                                        the true error
        lifted MC:  12.79 +/- 0.34   <- consistent with the reference,
                                        ~4.6x smaller error bar

   The plain-MC failure mode is exactly why lifting exists: the sector
   whose cleared polynomial has constant term 10^-4 contributes spikes of
   order (10^-4)^{-2} = 10^8 on a tiny region that uniform sampling almost
   never hits, so the estimate AND its error bar come out too small.

   Pushing the coefficient to 10^-6 makes even the lifted integrand
   heavy-tailed at moderate sample counts (the lifted decomposition still
   validates to <1% by quadrature, but finite-sample MC needs >> 10^6
   samples) — increase NSamples, or lift with a larger k, in that regime.

   NOTE: lifting requires REAL polynomial exponents (a complex effective
   exponent admits no real-valued domain indicator; ProcessSectorLifted
   fires TropicalEval::liftcomplex).  Requires g++.
   ============================================================================ *)

Print["=== Example 20: Lifting a small coefficient ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + 10^-4 x1^2 + x2^2 + x1*x2^2)^2"];
Print[];

Module[{poly, vars, spec, extreme, liftRules, lcRes, liftedSpec, liftData,
        verts, fanData, lverts, liftedFan, vl, ref, resPlain, resLifted,
        resAuto, workDir},

  poly = 1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  workDir = FileNameJoin[{Directory[], "INTERFILES"}];

  (* --- Step 1: detect the extreme coefficient (anchor chosen automatically) --- *)
  extreme = DetectExtremeCoefficients[spec, 1000];
  Print["DetectExtremeCoefficients (threshold 1000): ", extreme];
  Print["  kStar anchor rule auto-suggests k = ",
        First[#["SuggestedK"] & /@ extreme], " (z0 = 1/100); ",
        "the explicit rule below uses the same value."];
  Print[];

  (* --- Step 2: lift it with k = 2  ->  z0 = 1/100 (matches SuggestedK) --- *)
  liftRules = {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>};
  lcRes = LiftCoefficients[spec, liftRules];
  If[!AssociationQ[lcRes],
    Print["LiftCoefficients failed: ", lcRes],

    liftedSpec = lcRes["LiftedSpec"];
    liftData   = lcRes["LiftData"];
    Print["Anchor z0 = ", liftData["z0"],
          "   auxiliary variable: ", liftData["AuxVariable"]];
    Print["Lifted polynomial: ", liftedSpec["Polynomials"][[1]]];
    Print["Identity check (lifted poly at z -> z0 equals original): ",
      Expand[(liftedSpec["Polynomials"][[1]] /.
              liftData["AuxVariable"] -> liftData["z0"]) - poly] === 0];
    Print[];

    (* --- Step 3: lifted fan (automatic; polytope is full-dimensional) --- *)
    lverts = PolytopeVertices[
      (Times @@ liftedSpec["Polynomials"])^(-1), liftedSpec["Variables"]];
    liftedFan = ComputeDecomposition[lverts, "ShowProgress" -> False];
    Print["Lifted fan: ", Length[liftedFan[[1]]], " rays, ",
          Length[liftedFan[[2]]], " sectors (3 variables)"];
    Print[];

    (* --- Steps 4-5: exactness of the lifted decomposition --- *)
    Print["ValidateLiftedDecomposition (delta resolved per sector):"];
    vl = Quiet@ValidateLiftedDecomposition[
      spec, liftedSpec, liftedFan, liftData, {}, 4];
    If[AssociationQ[vl],
      Print["  Direct NIntegrate: ", vl["DirectResult"]];
      Print["  Lifted sector sum: ", vl["SectorSum"]];
      Print["  Relative error:    ", vl["RelativeError"]];
      Print["  Dropped (EmptyDomain) sectors: ",
            Length[Lookup[vl, "DroppedSectors", {}]]];
    ];
    Print[];

    (* --- Step 6: MC comparison, plain vs lifted, same sample count --- *)
    ref = Quiet@NIntegrate[
      1/(1 + 10^-4 t1^2 + t2^2 + t1 t2^2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 25, PrecisionGoal -> 8];
    Print["NIntegrate reference = ", ref];
    Print[];

    Print["Plain pipeline (no lifting), 2*10^5 samples:"];
    verts   = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1), vars];
    fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
    resPlain = EvaluateTropicalMC[spec, fanData, {{}},
      "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False,
      "WorkingDirectory" -> workDir];
    If[AssociationQ[resPlain],
      Print["  plain MC  = ", resPlain["Results"][[1]]["Re"],
            " +/- ", resPlain["Results"][[1]]["ReErr"],
            "   (deviation from ref: ",
            resPlain["Results"][[1]]["Re"] - ref, ")"];
      Print["  -> heavy-tailed: the deviation exceeds the sample error bar"];,
      Print["  plain pipeline failed (is g++ installed?)"];
    ];
    Print[];

    Print["Lifted pipeline (k=2), same 2*10^5 samples:"];
    resLifted = EvaluateTropicalMCLifted[spec, {{}},
      "LiftRules" -> liftRules,
      "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False,
      "WorkingDirectory" -> workDir];
    If[AssociationQ[resLifted],
      Print["  lifted MC = ", resLifted["Results"][[1]]["Re"],
            " +/- ", resLifted["Results"][[1]]["ReErr"],
            "   (deviation from ref: ",
            resLifted["Results"][[1]]["Re"] - ref, ")"];,
      Print["  lifted pipeline failed (is g++ installed?)"];
    ];
    Print[];

    (* --- Step 7: the anchor is now chosen automatically --- *)
    Print["Step 7: automatic anchor selection (the kStar rule)."];
    Print["  How DetectExtremeCoefficients sets SuggestedK for this 10^-4 coeff:"];
    Print["    AnchorRule -> \"kStar\" (default) : k = ",
          First[#["SuggestedK"] & /@ DetectExtremeCoefficients[spec, 1000]]];
    Print["    AnchorRule -> \"Unit\" (legacy)   : k = ",
          First[#["SuggestedK"] & /@
            DetectExtremeCoefficients[spec, 1000, "AnchorRule" -> "Unit"]]];
    Print["    BandEdgeGuard -> True (opt-in)   : k = ",
          First[#["SuggestedK"] & /@
            DetectExtremeCoefficients[spec, 1000, "BandEdgeGuard" -> True]]];
    Print[];
    Print["  Fully automatic lift (LiftRules -> Automatic; no explicit k):"];
    resAuto = EvaluateTropicalMCLifted[spec, {{}},
      "LiftRules" -> Automatic,
      "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False,
      "WorkingDirectory" -> workDir];
    If[AssociationQ[resAuto] && KeyExistsQ[resAuto, "Results"],
      Print["    automatic MC = ", resAuto["Results"][[1]]["Re"],
            " +/- ", resAuto["Results"][[1]]["ReErr"],
            "   (deviation from ref: ",
            resAuto["Results"][[1]]["Re"] - ref, ")"];
      Print["    -> same k=2 anchor as the hand-picked lift above."];,
      Print["    automatic lift failed (is g++ installed?)"];
    ];
  ];
];
Print[];


Print["=== Examples 15-20 complete ==="];
