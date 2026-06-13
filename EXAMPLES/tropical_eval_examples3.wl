(* ============================================================================
   tropical_eval_examples3.wl

   Examples 21-23: cross-checks via INDEPENDENT methods.

   Examples 1-20 validate the tropical pipeline against Mathematica's
   NIntegrate.  This file cross-checks against methods that share no code,
   no algorithm, and no CAS with either the tropical pipeline or NIntegrate:

     - an EXACT analytic closed form (Beta/Gamma functions, Example 21);
     - the CUBA library (https://feynarts.de/cuba/), integrating the
       ORIGINAL integrand directly on [0,inf)^n via the compactification
       x_i = t_i/(1-t_i)  (no tropical decomposition, no flattening):
         * Cuhre  -- deterministic globally adaptive cubature
                     (reference-grade accuracy in low dimensions)
         * Vegas  -- adaptive importance-sampling Monte Carlo
                     (independent MC algorithm and RNG)

   The integrands continue the themes of Examples 15-20: polynomials raised
   to complex exponents and/or carrying small coefficients.  All integrals
   are convergent.

   Requires: tropical_fan.wl, tropical_eval.wl, Polymake, g++.
   Optional:  CUBA  -- install with `brew install cuba` (macOS/Homebrew) or
   build from https://feynarts.de/cuba/.  The examples probe /opt/homebrew,
   /usr/local, /usr for cuba.h + libcuba and skip the CUBA parts (with a
   message) if it is not installed.

   Run as a script:
     wolframscript -file EXAMPLES/tropical_eval_examples3.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["tropical_eval.wl loaded successfully"];
Print[];

(* CUBA cross-check helpers (generateCubaSource, runCubaCheck):
   shared with test_cuba.wl — see EXAMPLES/cuba_common.wl for the
   generator documentation. *)
Get[FileNameJoin[{Directory[], "EXAMPLES", "cuba_common.wl"}]];

interfilesDir = FileNameJoin[{Directory[], "INTERFILES"}];


(* ============================================================================
   Example 21: Exact analytic cross-check (Beta/Gamma), complex exponents

   Integral[0,Inf] dx1 dx2 x1^{a1-1} x2^{a2-1} (1 + x1 + x2)^{-b}
       = Gamma[a1] Gamma[a2] Gamma[b-a1-a2] / Gamma[b]

   with  a1 = 5/4 + I/3,  a2 = 3/2 - I/5,  b = 5 + I/2.

   This generalized Beta integral has an EXACT closed form valid for
   complex parameters (convergence: Re(a1), Re(a2) > 0 and
   Re(b - a1 - a2) > 0; here 1.25, 1.5 and 2.25).  It is the strongest
   possible reference: every numerical method can be measured against it.

   Four independent computations of the same number:
     1. exact Gamma-function formula        (analytic)
     2. direct NIntegrate + sector sum      (Mathematica quadrature)
     3. tropical Monte Carlo                (this package's C++ pipeline)
     4. CUBA Cuhre on the direct integrand  (independent C++ cubature)

   Measured at the time of writing: Cuhre agrees with the exact value to
   ~3e-8 relative; the tropical MC at 5*10^5 samples lands well within
   its ~1e-4 error bars.
   ============================================================================ *)

Print["=== Example 21: Exact Beta/Gamma cross-check, complex exponents ==="];
Print["Integral[0,Inf] dx1 dx2 x1^{a1-1} x2^{a2-1} (1+x1+x2)^{-b}"];
Print["a1 = 5/4 + I/3,  a2 = 3/2 - I/5,  b = 5 + I/2"];
Print[];

Module[{a1, a2, b, exact, poly, vars, spec, verts, fanData, vr, mcRes, cuba,
        report},

  a1 = 5/4 + I/3;  a2 = 3/2 - I/5;  b = 5 + I/2;
  exact = N[Gamma[a1] Gamma[a2] Gamma[b - a1 - a2] / Gamma[b]];

  poly = 1 + x[1] + x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {a1 - 1, a2 - 1},   (* x^{a-1} *)
    "PolynomialExponents" -> {-b},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  Print["Exact value: ", exact];
  Print[];

  verts   = PolytopeVertices[poly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  Print["Fan: ", Length[fanData[[1]]], " rays, ",
        Length[fanData[[2]]], " sectors (all convergent, complex effA)"];

  report[label_, val_, err_: None] := Print["  ",
    StringPadRight[label, 22], val,
    If[err =!= None, "  +/- " <> ToString[err, InputForm], ""],
    "   |rel dev| = ", Abs[(val - exact)/exact]];

  (* 2. Mathematica quadrature: direct + tropical sector sum *)
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    report["NIntegrate direct:", vr["DirectResult"]];
    report["Sector sum (NInt):", vr["SectorSum"]];
  ];

  (* 3. Tropical Monte Carlo (this package's C++ pipeline) *)
  mcRes = EvaluateTropicalMC[spec, fanData, {{}},
    "NSamples" -> 500000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> interfilesDir];
  If[AssociationQ[mcRes],
    Module[{r = mcRes["Results"][[1]]},
      report["Tropical MC:", r["Re"] + I r["Im"],
             {r["ReErr"], r["ImErr"]}];
    ],
    Print["  Tropical MC failed (is g++ installed?)"];
  ];

  (* 4. CUBA Cuhre on the direct integrand *)
  cuba = runCubaCheck[spec, "ex21", "RunVegas" -> False];
  If[AssociationQ[cuba] && KeyExistsQ[cuba, "CUHRE"],
    report["CUBA Cuhre (direct):", cuba["CUHRE"]["Value"],
           cuba["CUHRE"]["Error"]];
  ];
];
Print[];


(* ============================================================================
   Example 22: CUBA cross-check of the small-coefficient integral

   Integral[0,Inf] dx1 dx2 / (1 + 10^-4 x1^2 + x2^2 + x1*x2^2)^2

   The integrand of Example 20, where plain uniform tropical MC at 2*10^5
   samples was shown to underestimate badly (7.83 +/- 1.56 vs ~13.09) due
   to a heavy-tailed sector.  Here the same number is computed by three
   methods that have nothing in common algorithmically:

     - CUBA Cuhre   (deterministic adaptive cubature, direct integrand):
                    reference-grade, agrees with NIntegrate to ~1e-8
     - CUBA Vegas   (adaptive importance-sampling MC, direct integrand):
                    an INDEPENDENT MC that adapts its sampling density,
                    so it copes with the small coefficient where plain
                    uniform sampling fails
     - lifted tropical MC (k=2 auxiliary-variable lifting, Example 20)

   Note on the lifted MC error bars: the lifted integrand still carries
   some residual tail for this integrand (the sample error decreases
   slower than 1/sqrt(N) between 2*10^5 and 10^6 samples), so expect
   agreement at the 1-2 sigma level rather than perfectly calibrated
   errors.  Cuhre is the referee.
   ============================================================================ *)

Print["=== Example 22: CUBA cross-check, small coefficient (10^-4) ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + 10^-4 x1^2 + x2^2 + x1*x2^2)^2"];
Print[];

Module[{poly, vars, spec, cuba, resLifted, cuhreVal},

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

  (* CUBA: Vegas + Cuhre on the direct integrand *)
  cuba = runCubaCheck[spec, "ex22"];
  cuhreVal = None;
  If[AssociationQ[cuba],
    If[KeyExistsQ[cuba, "CUHRE"],
      cuhreVal = Re[cuba["CUHRE"]["Value"]];
      Print["  CUBA Cuhre (direct):  ", cuhreVal,
            " +/- ", Re[cuba["CUHRE"]["Error"]],
            "   (", cuba["CUHRE"]["NEval"], " evals, deterministic)"];
    ];
    If[KeyExistsQ[cuba, "VEGAS"],
      Print["  CUBA Vegas (direct):  ", Re[cuba["VEGAS"]["Value"]],
            " +/- ", Re[cuba["VEGAS"]["Error"]],
            "   (", cuba["VEGAS"]["NEval"], " evals, importance sampling)"];
    ];
  ];

  (* Lifted tropical MC (cf. Example 20) *)
  resLifted = EvaluateTropicalMCLifted[spec, {{}},
    "LiftRules" -> {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
    "NSamples" -> 200000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> interfilesDir];
  If[AssociationQ[resLifted],
    Module[{r = resLifted["Results"][[1]]},
      Print["  Lifted tropical MC:   ", r["Re"], " +/- ", r["ReErr"],
            "   (2*10^5 samples)"];
      If[NumericQ[cuhreVal],
        Print["    deviation from Cuhre = ", r["Re"] - cuhreVal,
              "  (", Abs[r["Re"] - cuhreVal]/r["ReErr"], " sigma)"]];
    ],
    Print["  Lifted tropical MC failed (is g++ installed?)"];
  ];

  Print[];
  Print["  (Plain uniform tropical MC at this sample count underestimates"];
  Print["   badly -- see Example 20.  Vegas succeeds on the same direct"];
  Print["   integrand because it adapts its sampling density.)"];
];
Print[];


(* ============================================================================
   Example 23: CUBA cross-check with complex exponent + small coefficients

   Integral[0,Inf] dx1 dx2 / (1 + 10^-3 x1^2 + x2^2 + 10^-2 x1*x2^2)^{2+I/2}

   The integrand of Example 19 (small coefficients AND a complex
   polynomial exponent).  Lifting is not available here (complex
   exponents have no real-valued domain indicator), so the tropical MC
   runs unlifted -- its error bars are honest but sizable in the
   heavy-tail sector.  CUBA Cuhre integrates the direct integrand with
   two components (Re, Im) and provides the independent reference,
   agreeing with NIntegrate (Example 19) to ~1e-7.
   ============================================================================ *)

Print["=== Example 23: CUBA cross-check, complex exponent + small coeffs ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + 10^-3 x1^2 + x2^2 + 10^-2 x1*x2^2)^{2+I/2}"];
Print[];

Module[{poly, vars, spec, verts, fanData, cuba, cuhreVal, mcRes},

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

  (* CUBA Cuhre, complex result *)
  cuba = runCubaCheck[spec, "ex23", "RunVegas" -> False];
  cuhreVal = None;
  If[AssociationQ[cuba] && KeyExistsQ[cuba, "CUHRE"],
    cuhreVal = cuba["CUHRE"]["Value"];
    Print["  CUBA Cuhre (direct):  ", cuhreVal];
    Print["                 +/-    ", ToString[cuba["CUHRE"]["Error"], InputForm],
          "   (", cuba["CUHRE"]["NEval"], " evals)"];
  ];

  (* Tropical MC, unlifted (complex exponent forbids lifting) *)
  verts   = PolytopeVertices[(1 + x[1]^2 + x[2]^2 + x[1] x[2]^2)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  mcRes = EvaluateTropicalMC[spec, fanData, {{}},
    "NSamples" -> 2000000, "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> interfilesDir];
  If[AssociationQ[mcRes],
    Module[{r = mcRes["Results"][[1]], mc, sigma},
      mc = r["Re"] + I r["Im"];
      sigma = Sqrt[r["ReErr"]^2 + r["ImErr"]^2];
      Print["  Tropical MC:          ", mc];
      Print["                 +/-    (", r["ReErr"], ", ", r["ImErr"],
            ")   (2*10^6 samples)"];
      If[cuhreVal =!= None,
        Print["    |deviation| from Cuhre = ", Abs[mc - cuhreVal],
              "  (", Abs[mc - cuhreVal]/sigma, " sigma-combined)"]];
    ],
    Print["  Tropical MC failed (is g++ installed?)"];
  ];
];
Print[];


Print["=== Examples 21-23 complete ==="];
