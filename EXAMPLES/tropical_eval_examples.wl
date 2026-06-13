(* ============================================================================
   tropical_eval_examples.wl

   Examples demonstrating the tropical_eval evaluation pipeline.
   Requires: tropical_fan.wl, tropical_eval.wl, Polymake, g++ with OpenMP.

   Load the package first:
     SetDirectory[FileNameJoin[{NotebookDirectory[], ".."}]];
     Get["tropical_eval.wl"];

   Or run as a script:
     wolframscript -file EXAMPLES/tropical_eval_examples.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["tropical_eval.wl loaded successfully"];
Print[];


(* ============================================================================
   Example 1: ProcessSector — basic 2D convergent integral

   Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3

   This is the simplest case: one polynomial with positive coefficients,
   real negative exponent, no kinematic parameters, no epsilon regulator.
   ============================================================================ *)

Print["=== Example 1: Basic 2D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3"];
Print[];

Module[{poly, vars, spec, verts, fanData, dualVerts, simplices},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};

  (* IntegrandSpec: the standard input format for all pipeline functions *)
  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},          (* no x_i^{A_i} prefactor *)
    "PolynomialExponents" -> {-3},            (* P^{-3} *)
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},               (* no parameters *)
    "RegulatorSymbol"    -> None              (* no epsilon *)
  |>;

  (* Step 1: Compute the tropical fan *)
  verts    = PolytopeVertices[poly^(-3), vars];
  fanData  = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dualVerts, simplices} = fanData;

  Print["Fan: ", Length[dualVerts], " rays, ",
        Length[simplices], " sectors"];

  (* Step 2: Process each sector *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, dualVerts, simplices[[s]], s,
                         "Verbose" -> True];
      If[AssociationQ[sd],
        Print["  Effective exponents: ", sd["NewExponents"]];
        Print["  Prefactor: ", sd["Prefactor"]];
        Print["  # monomials in Q: ",
              Length /@ sd["FlattenedPolys"]];
        Print["  Dominant monomial (should have exps=0): ",
              SelectFirst[sd["FlattenedPolys"][[1]],
                AllTrue[#[[2]], TrueQ[# == 0] &] &]];
        Print[];
      ];
    ],
    {s, Length[simplices]}
  ];

  (* Step 3: Validate the decomposition against NIntegrate *)
  Print["Validating against NIntegrate..."];
  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 4];
    Print["  Direct NIntegrate:  ", vr["DirectResult"]];
    Print["  Sector sum:         ", vr["SectorSum"]];
    Print["  Relative error:     ", vr["RelativeError"]];
    Print["  Per-sector results: ", vr["SectorResults"]];
  ];
];
Print[];


(* ============================================================================
   Example 2: Complex exponents and the flattening check

   Integral[0,Inf] dx1 dx2 / (1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2)^{2+i}

   Demonstrates complex polynomial exponents. The polynomial always
   evaluates to a positive real (all coefficients positive), so
   P^{2+i} = exp((2+i) * log(P)) is well-defined.
   ============================================================================ *)

Print["=== Example 2: Complex exponents ==="];
Print["Integral[0,Inf] dx1 dx2 / (1+2x1^2+x2^2+x1*x2^2+3x1^2*x2)^{2+I}"];
Print[];

Module[{poly, vars, A, spec, verts, fanData, sectorData},

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};
  A    = 2 + I;

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  (* Process sector 1 and spot-check the integrand magnitude *)
  sectorData = ProcessSector[spec, fanData[[1]], fanData[[2, 1]], 1];

  Print["Sector 1 details:"];
  Print["  Effective exponents (complex): ", sectorData["NewExponents"]];
  Print["  Prefactor (complex): ", sectorData["Prefactor"]];

  Print[];
  Print["Flattening magnitude check (should be O(1)):"];
  Module[{check},
    check = CheckFlatteningMagnitude[sectorData, 10];
    Print["  Mean |integrand|: ", check["Mean"]];
    Print["  Min:  ", check["Min"]];
    Print["  Max:  ", check["Max"]];
  ];

  Print[];
  Print["Full validation:"];
  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 3: Kinematic-dependent integral

   Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2

   The polynomial coefficient 'lam' is a kinematic parameter that varies.
   This demonstrates how to set up kinematic symbols and evaluate at
   multiple parameter values.
   ============================================================================ *)

Print["=== Example 3: Kinematic-dependent integral ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2"];
Print[];

Module[{lam, poly, vars, spec, verts, fanData, lamValues},

  lam  = Symbol["lam"];
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  (* Compute fan at lam=1 (topology doesn't change with lam) *)
  verts   = PolytopeVertices[(poly /. lam -> 1)^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  (* Process one sector with symbolic lam to see the structure *)
  Module[{sd},
    sd = ProcessSector[spec, fanData[[1]], fanData[[2, 1]], 1];
    Print["Sector 1 (symbolic):"];
    Print["  Effective exponents: ", sd["NewExponents"]];
    Print["  Prefactor: ", sd["Prefactor"]];
    Print["  First monomial coeff: ", sd["FlattenedPolys"][[1, 1, 1]]];
    Print[];
  ];

  (* Evaluate at several values of lam *)
  lamValues = {0.5, 1.0, 2.0, 5.0};

  Print["Validation at multiple kinematic points:"];
  Do[
    Module[{kinRules, vr},
      kinRules = {lam -> lamVal};
      vr = Quiet@ValidateDecomposition[spec, fanData, kinRules, 3];
      If[AssociationQ[vr],
        Print["  lam = ", lamVal,
              ":  direct = ", vr["DirectResult"],
              "  sector sum = ", vr["SectorSum"],
              "  rel err = ", vr["RelativeError"]];
      ];
    ],
    {lamVal, lamValues}
  ];
];
Print[];


(* ============================================================================
   Example 4: Inspecting the tropical factoring

   This example shows what happens inside ProcessSector step by step:
   the monomial substitution, the tropical factoring (extracting the
   dominant monomial), and the resulting cleared polynomials.

   Polynomial: P = 1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2
   Rays: rho_1 = (0,1), rho_2 = (1,1)  (sector 5 of the fan)
   ============================================================================ *)

Print["=== Example 4: Tropical factoring step by step ==="];
Print[];

Module[{poly, vars, spec, dualVerts, simplex, sd},

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};

  dualVerts = {{0, 1}, {1, 1}};
  simplex   = {1, 2};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-Symbol["A"]},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  sd = ProcessSector[spec, dualVerts, simplex, 1, "Verbose" -> True];

  Print[];
  Print["Step-by-step breakdown:"];

  Print["  Rays: rho1 = ", dualVerts[[1]], ", rho2 = ", dualVerts[[2]]];
  Print["  M = -Transpose[rays] = ", sd["RayMatrix"]];
  Print["  det(M) = ", sd["DetM"]];
  Print[];

  Print["  Raw a (from monomial exponents only):"];
  Print["    a = (monoExps+1) . M = ", sd["RawExponents"]];
  Print[];

  Print["  Transformed polynomial monomials (in y-space):"];
  Do[
    Print["    ", mono[[1]], " * y^", mono[[2]]],
    {mono, sd["TransformedPolys"][[1]]}
  ];
  Print[];

  Print["  Min exponents (dominant monomial): ", sd["MinExponents"]];
  Print[];

  Print["  Cleared polynomial Q (non-negative exponents):"];
  Do[
    Print["    ", mono[[1]], " * y^", mono[[2]]],
    {mono, sd["ClearedPolys"][[1]]}
  ];
  Print[];

  Print["  Effective exponents:"];
  Print["    a_eff = rawA + B*minExp = ", sd["NewExponents"]];
  Print["    (For A=2: a_eff = ", sd["NewExponents"] /. Symbol["A"] -> 2, ")"];
  Print["    (For A=3: a_eff = ", sd["NewExponents"] /. Symbol["A"] -> 3, ")"];
  Print[];

  Print["  Flattened polynomial (exponents / a_eff):"];
  Do[
    Print["    ", mono[[1]], " * y'^", mono[[2]]],
    {mono, sd["FlattenedPolys"][[1]]}
  ];
];
Print[];


(* ============================================================================
   Example 5: MmaToC — Mathematica to C++ conversion

   Demonstrates the expression converter used by the C++ code generator.
   ============================================================================ *)

Print["=== Example 5: MmaToC expression converter ==="];
Print[];

Module[{lam, paramMap},

  lam = Symbol["lam"];
  paramMap = <|lam -> "params[0]"|>;

  Print["Integer:      ", MmaToC[42]];
  Print["Rational:     ", MmaToC[3/7]];
  Print["Real:         ", MmaToC[3.14159]];
  Print["Complex:      ", MmaToC[2 + 3 I]];
  Print["Power:        ", MmaToC[Power[x, 3]]];
  Print["Sqrt:         ", MmaToC[Sqrt[x]]];
  Print["Log:          ", MmaToC[Log[x]]];
  Print["Exp:          ", MmaToC[Exp[x]]];
  Print["Sum:          ", MmaToC[a + b + c]];
  Print["Product:      ", MmaToC[a * b * c]];
  Print["With params:  ", MmaToC[1 + 2 lam, paramMap]];
  Print["Complex expr: ", MmaToC[Exp[(2 + I) * Log[1 + lam]], paramMap]];
];
Print[];


(* ============================================================================
   Example 6: C++ code generation (convergent, no kinematics)

   Generates a C++ Monte Carlo source file for a simple convergent integral,
   compiles it, and runs it.

   Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3  =  Pi/8
   ============================================================================ *)

Print["=== Example 6: C++ code generation and execution ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3 = Pi/8"];
Print[];

Module[{poly, vars, spec, verts, fanData, allSectors,
        convergent, cppResult, cppFile, exact},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  exact = Pi/8;

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  (* Process all sectors *)
  allSectors = Table[
    ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
    {s, Length[fanData[[2]]]}
  ];

  convergent = Select[allSectors,
    AssociationQ[#] && !#["IsDivergent"] &];

  Print["Processed ", Length[convergent], " convergent sectors"];

  (* Generate C++ *)
  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile = FileNameJoin[{Directory[], "INTERFILES", "example6_mc.cpp"}];
  cppResult = GenerateCppMonteCarlo[
    convergent, {},   (* no divergent sectors *)
    spec, cppFile,
    "NSamples" -> 500000
  ];

  If[AssociationQ[cppResult],
    Print[];

    (* Compile *)
    Module[{binary, kinFile, resultFile, runResult},
      binary     = FileNameJoin[{Directory[], "INTERFILES", "example6_mc"}];
      kinFile    = FileNameJoin[{Directory[], "INTERFILES", "example6_kin.txt"}];
      resultFile = FileNameJoin[{Directory[], "INTERFILES", "example6_results.txt"}];

      (* No kinematics: write count = 1 *)
      Export[kinFile, "1\n", "Text"];

      If[CompileCpp[cppFile, binary, False] =!= $Failed,
        Print[];
        runResult = RunProcess[{binary, kinFile, resultFile,
                                "500000", "2"}];
        If[runResult["ExitCode"] == 0,
          Module[{output, vals, mcVal, mcErr, relErr},
            output = Import[resultFile, "Text"];
            vals   = ToExpression /@ StringSplit[StringTrim[output]];
            mcVal  = vals[[1]] + I vals[[2]];
            mcErr  = vals[[3]] + I vals[[4]];

            Print["Results:"];
            Print["  Exact   = ", N[exact]];
            Print["  MC      = ", mcVal, " +/- ", Abs[mcErr]];
            Print["  Rel err = ",
                  Abs[(Re[mcVal] - N[exact]) / N[exact]]];
          ];,
          Print["Run failed: ", runResult["StandardError"]];
        ];

        (* Clean up *)
        Quiet[DeleteFile /@ {cppFile, binary, kinFile, resultFile}];
      ];
    ];
  ];
];
Print[];


(* ============================================================================
   Example 7: C++ code generation with kinematic parameters

   Generates C++ code for a kinematic-dependent integral and evaluates
   at 10 values of lambda.

   Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2
   ============================================================================ *)

Print["=== Example 7: C++ with kinematic scan ==="];
Print["Integral[0,Inf] dx1 dx2 / (1+lam*x1^2+x2^2+x1*x2^2)^2"];
Print["10 values of lam in [0.5, 5]"];
Print[];

Module[{lam, poly, vars, spec, verts, fanData, allSectors,
        convergent, cppResult, cppFile,
        lamValues, kinPoints},

  lam  = Symbol["lam"];
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  lamValues = Table[0.5 + 4.5 (i - 1)/9, {i, 10}];
  kinPoints = List /@ lamValues;

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[(poly /. lam -> 1)^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  allSectors = Table[
    ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
    {s, Length[fanData[[2]]]}
  ];

  convergent = Select[allSectors,
    AssociationQ[#] && !#["IsDivergent"] &];

  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile = FileNameJoin[{Directory[], "INTERFILES", "example7_mc.cpp"}];
  cppResult = GenerateCppMonteCarlo[
    convergent, {}, spec, cppFile,
    "NSamples" -> 200000
  ];

  If[AssociationQ[cppResult],
    Module[{binary, kinFile, resultFile, runResult},
      binary     = FileNameJoin[{Directory[], "INTERFILES", "example7_mc"}];
      kinFile    = FileNameJoin[{Directory[], "INTERFILES", "example7_kin.txt"}];
      resultFile = FileNameJoin[{Directory[], "INTERFILES", "example7_results.txt"}];

      (* Write kinematic data: one lam value per line *)
      Export[kinFile,
        StringRiffle[ToString[CForm[#]] & /@ lamValues, "\n"] <> "\n",
        "Text"
      ];

      If[CompileCpp[cppFile, binary, False] =!= $Failed,
        Print[];
        runResult = RunProcess[{binary, kinFile, resultFile,
                                "200000", "2"}];
        If[runResult["ExitCode"] == 0,
          Module[{lines, parsed},
            lines  = Select[
              StringSplit[Import[resultFile, "Text"], "\n"],
              StringLength[StringTrim[#]] > 0 &
            ];
            parsed = (ToExpression /@ StringSplit[#]) & /@ lines;

            Print["Results (MC vs NIntegrate):"];
            Print[StringPadRight["  lam", 10],
                  StringPadRight["MC Re", 16],
                  StringPadRight["MC err", 14],
                  "NIntegrate"];
            Do[
              Module[{mcRe, mcErr, niResult},
                mcRe  = parsed[[i, 1]];
                mcErr = parsed[[i, 3]];
                niResult = Quiet@NIntegrate[
                  1 / (1 + lamValues[[i]] t1^2 + t2^2 + t1 t2^2)^2,
                  {t1, 0, Infinity}, {t2, 0, Infinity},
                  PrecisionGoal -> 5
                ];
                Print["  ",
                  StringPadRight[ToString@NumberForm[lamValues[[i]], {4,2}], 8],
                  StringPadRight[ToString@NumberForm[mcRe, {8,6}], 16],
                  StringPadRight[ToString[mcErr], 14],
                  ToString@NumberForm[niResult, {8,6}]];
              ],
              {i, Length[lamValues]}
            ];
          ];,
          Print["Run failed: ", runResult["StandardError"]];
        ];

        Quiet[DeleteFile /@ {cppFile, binary, kinFile, resultFile}];
      ];
    ];
  ];
];
Print[];


(* ============================================================================
   Example 8: ParsePolynomial utility

   Shows how ParsePolynomial extracts {coefficient, exponentVector} pairs
   from a symbolic polynomial.
   ============================================================================ *)

Print["=== Example 8: ParsePolynomial utility ==="];
Print[];

Module[{lam, poly, vars, parsed},

  lam  = Symbol["lam"];
  poly = 1 + 3 lam x[1]^2 + x[2]^3 + 2 x[1] x[2];
  vars = {x[1], x[2]};

  parsed = ParsePolynomial[poly, vars];

  Print["Polynomial: ", poly];
  Print["Variables:  ", vars];
  Print[];
  Print["Parsed monomials {coefficient, exponents}:"];
  Do[
    Print["  ", mono[[1]], " * x^", mono[[2]]],
    {mono, parsed}
  ];
];
Print[];


(* ============================================================================
   Example 9: Multiple polynomials (product of two factors)

   Integral[0,Inf] dx1 dx2 x1^{1/2} / [(1+x1+x2)^2 * (1+x1*x2)]

   Shows how to handle a product of polynomials with separate exponents,
   plus a non-trivial monomial prefactor x1^{1/2}.
   The two factors have different Newton polytopes, ensuring a
   full-dimensional combined polytope.
   ============================================================================ *)

Print["=== Example 9: Multiple polynomials ==="];
Print["Integral[0,Inf] dx1 dx2 x1^{1/2} / [(1+x1+x2)^2 * (1+x1*x2)]"];
Print[];

Module[{p1, p2, vars, spec, integrand, verts, fanData, vr},

  p1   = 1 + x[1] + x[2];
  p2   = 1 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {p1, p2},
    "MonomialExponents"  -> {1/2, 0},        (* x1^{1/2} prefactor *)
    "PolynomialExponents" -> {-2, -1},        (* P1^{-2} * P2^{-1} *)
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* Combined integrand for fan computation *)
  integrand = p1^(-2) * p2^(-1);
  verts     = PolytopeVertices[integrand, vars];
  fanData   = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  (* Process sectors *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 10: Non-trivial numerator — single numerator polynomial

   Integral[0,Inf] dx1 dx2  (1 + x1)^{3/2} / (1 + x1^2 + x2^2 + x1*x2)^3

   The integrand has a POSITIVE polynomial exponent (numerator) and a
   NEGATIVE one (denominator).  The IntegrandSpec lists both polynomials
   with their signed exponents.

   For the Newton polytope / fan computation, we use the product of all
   polynomials raised to |B_j| (absolute exponents), since the tropical
   fan only depends on which monomials dominate, not on the sign of the
   power.
   ============================================================================ *)

Print["=== Example 10: Numerator polynomial (1+x1)^{3/2} ==="];
Print["Integral[0,Inf] dx1 dx2 (1+x1)^{3/2} / (1+x1^2+x2^2+x1*x2)^3"];
Print[];

Module[{pNum, pDen, vars, spec, fanPoly, verts, fanData, vr},

  pNum = 1 + x[1];
  pDen = 1 + x[1]^2 + x[2]^2 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {pNum, pDen},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {3/2, -3},        (* numerator^{3/2} * denom^{-3} *)
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* For the fan: the normal fan of a Minkowski sum a*Newt(P1) + b*Newt(P2)
     is the same for any positive a, b.  So just multiply the polynomials
     once -- all that matters is that every monomial appears. *)
  fanPoly = pNum * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  (* Process and validate *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"],
              ", prefactor = ", sd["Prefactor"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 11: Non-trivial numerator — two numerator factors

   Integral[0,Inf] dx1 dx2  (1 + x1 + x2)^{1/2} * (x1 + x2)^{1/3}
                             / (1 + 2*x1^2 + x2^2 + x1*x2^2)^2

   Three polynomials: two in the numerator (positive exponents 1/2, 1/3)
   and one in the denominator (negative exponent -2).
   ============================================================================ *)

Print["=== Example 11: Two numerator factors ==="];
Print["Integral[0,Inf] dx1 dx2 (1+x1+x2)^{1/2}*(x1+x2)^{1/3} / (1+2x1^2+x2^2+x1*x2^2)^2"];
Print[];

Module[{pN1, pN2, pDen, vars, spec, fanPoly, verts, fanData, vr},

  pN1  = 1 + x[1] + x[2];
  pN2  = x[1] + x[2];
  pDen = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {pN1, pN2, pDen},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {1/2, 1/3, -2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* Fan from product of all polynomials *)
  fanPoly = pN1 * pN2 * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 12: Numerator with kinematic-dependent coefficients

   Integral[0,Inf] dx1 dx2  (1 + lam*x1)^{3/2} / (1 + x1^2 + x2^2)^3

   The numerator coefficient 'lam' varies.  One compile, multiple evaluations.
   ============================================================================ *)

Print["=== Example 12: Numerator with kinematic coefficients ==="];
Print["Integral[0,Inf] dx1 dx2 (1+lam*x1)^{3/2} / (1+x1^2+x2^2)^3"];
Print[];

Module[{lam, pNum, pDen, vars, spec, fanPoly, verts, fanData},

  lam  = Symbol["lam"];
  pNum = 1 + lam x[1];
  pDen = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {pNum, pDen},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {3/2, -3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  (* Fan: use unit-coefficient polynomial since fan is independent of lam *)
  fanPoly = (1 + x[1]) * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];
  Print[];

  (* Validate at several lambda values *)
  Print["Validation at multiple kinematic points:"];
  Do[
    Module[{kinRules, vr},
      kinRules = {lam -> lamVal};
      vr = Quiet@ValidateDecomposition[spec, fanData, kinRules, 3];
      If[AssociationQ[vr],
        Print["  lam = ", lamVal,
              ":  direct = ", vr["DirectResult"],
              "  sector sum = ", vr["SectorSum"],
              "  rel err = ", vr["RelativeError"]];
      ];
    ],
    {lamVal, {0.5, 1., 2., 5.}}
  ];
];
Print[];


(* ============================================================================
   Example 13: 3D convergent integral

   Integral[0,Inf] dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3

   First example in dimension > 2.  Tests that ProcessSector, tropical
   factoring, and ValidateDecomposition all work correctly in 3D.
   ============================================================================ *)

Print["=== Example 13: 3D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3"];
Print[];

Module[{poly, vars, spec, verts, fanData, vr},

  poly = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  vars = {x[1], x[2], x[3]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[1]]], " rays, ", Length[fanData[[2]]], " sectors"];

  (* Show sector details *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"],
              ", prefactor = ", sd["Prefactor"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 14: 4D convergent integral

   Integral[0,Inf] dx1 dx2 dx3 dx4
       / (1 + x1^2 + x2^2 + x3^2 + x4^2 + x1*x2 + x3*x4)^4

   Tests the pipeline in 4 dimensions.  The Newton polytope is richer
   (mix of degree-1 and degree-2 monomials), producing more sectors.
   ============================================================================ *)

Print["=== Example 14: 4D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 dx3 dx4 / (1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4)^4"];
Print[];

Module[{poly, vars, spec, verts, fanData, allSectors, convergent, vr},

  poly = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1] x[2] + x[3] x[4];
  vars = {x[1], x[2], x[3], x[4]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0, 0, 0},
    "PolynomialExponents" -> {-4},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[1]]], " rays, ", Length[fanData[[2]]], " sectors"];

  (* Process all sectors *)
  allSectors = Table[
    ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
    {s, Length[fanData[[2]]]}
  ];
  convergent = Select[allSectors, AssociationQ[#] && !#["IsDivergent"] &];
  Print["Convergent sectors: ", Length[convergent], " / ", Length[fanData[[2]]]];

  (* Show a few sectors *)
  Do[
    Print["Sector ", convergent[[i]]["ConeIndex"],
          ": effA = ", convergent[[i]]["NewExponents"],
          ", prefactor = ", convergent[[i]]["Prefactor"]],
    {i, Min[4, Length[convergent]]}
  ];
  If[Length[convergent] > 4, Print["  ... (", Length[convergent] - 4, " more)"]];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 2];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


Print["=== All examples complete ==="];


(* ============================================================================
   MODULE 5: VALIDATION TESTS
   ============================================================================ *)

RunAllTests[] := Module[
  {results, nPass, nFail},

  results = {};
  nPass = 0;
  nFail = 0;

  Print[""];
  Print["================================================================"];
  Print["  TropicalEval Validation Suite (v2: convergent integrals only)"];
  Print["  Tests: 1, 2, 3v2, 5, 6, 7"];
  Print["================================================================"];
  Print[""];

  Module[{pass},
    pass = RunTest1[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 1", pass}];
  ];

  Module[{pass},
    pass = RunTest2[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 2", pass}];
  ];

  Module[{pass},
    pass = RunTest3v2[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 3v2", pass}];
  ];

  Module[{pass},
    pass = RunTest5[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 5", pass}];
  ];

  Module[{pass},
    pass = RunTest6[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 6", pass}];
  ];

  Module[{pass},
    pass = RunTest7[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 7", pass}];
  ];

  Print[""];
  Print["================================================================"];
  Print["  Results: ", nPass, " PASSED, ", nFail, " FAILED"];
  Print["================================================================"];

  results
];

(* --------------------------------------------------------------------------
   Test 1: Convergent 2D, real exponents
   -------------------------------------------------------------------------- *)

RunTest1[] := Module[
  {poly, vars, integrandSpec, verts, fanData,
   testAValues, allPass},

  Print["--- Test 1: Convergent 2D, real exponents ---"];
  Print["Integral[0,Inf] dx1 dx2 / (1+2x1^2+x2^2+x1*x2^2+3x1^2*x2)^A"];

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};

  testAValues = {2, 3};
  allPass = True;

  Do[
    Module[{A, spec, directResult, mcResult, relErr, pass},
      A = testA;

      spec = <|
        "Polynomials"       -> {poly},
        "MonomialExponents" -> {0, 0},
        "PolynomialExponents" -> {-A},
        "Variables"         -> vars,
        "KinematicSymbols"  -> {},
        "RegulatorSymbol"   -> None
      |>;

      directResult = NIntegrate[
        1 / (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^A,
        {t1, 0, Infinity}, {t2, 0, Infinity},
        MaxRecursion -> 20, PrecisionGoal -> 6
      ];

      verts = PolytopeVertices[poly^(-A), vars];
      fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

      Module[{vr},
        vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
        If[AssociationQ[vr],
          mcResult = vr["SectorSum"];
          relErr   = vr["RelativeError"];,
          mcResult = 0;
          relErr = Infinity;
        ];
      ];

      pass = NumericQ[relErr] && (relErr < 0.02);
      Print["  A = ", A, ":"];
      Print["    NIntegrate = ", directResult];
      Print["    Sector sum = ", mcResult];
      Print["    Rel error  = ", relErr];
      Print["    ", If[pass, "PASS", "FAIL"]];

      If[!pass, allPass = False];
    ],
    {testA, testAValues}
  ];

  allPass
];

(* --------------------------------------------------------------------------
   Test 2: Convergent 2D, complex exponents
   -------------------------------------------------------------------------- *)

RunTest2[] := Module[
  {poly, vars, A, spec, verts, fanData,
   directResult, mcResult, relErr, pass},

  Print["--- Test 2: Convergent 2D, complex exponents ---"];
  Print["Same integral, A = 2 + 0.5I"];

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};
  A    = 2 + 0.5 I;

  spec = <|
    "Polynomials"       -> {poly},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> vars,
    "KinematicSymbols"  -> {},
    "RegulatorSymbol"   -> None
  |>;

  directResult = NIntegrate[
    1 / (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^A,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 5
  ];

  verts   = PolytopeVertices[poly^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
    If[AssociationQ[vr],
      mcResult = vr["SectorSum"];
      relErr   = Abs[(mcResult - directResult) / directResult];,
      mcResult = 0;
      relErr = Infinity;
    ];
  ];

  pass = NumericQ[relErr] && (relErr < 0.02);
  Print["  NIntegrate = ", directResult];
  Print["  Sector sum = ", mcResult];
  Print["  Rel error  = ", relErr];
  Print["  ", If[pass, "PASS", "FAIL"]];

  pass
];

(* --------------------------------------------------------------------------
   Test 3v2 (negative test): regulator guard and divergent-input guard
   -------------------------------------------------------------------------- *)

RunTest3v2[] := Module[
  {allPass, eps},

  Print["--- Test 3v2: Negative tests (regulator and divergent input) ---"];
  allPass = True;

  (* Part A: eps-regulated spec -> must get $Failed with TropicalEval::noregulator *)
  Module[{pass, result},
    eps = Symbol["epsTest3v2"];
    Module[{spec, verts, fanData},
      spec = <|
        "Polynomials"        -> {1 + x[1]^2},
        "MonomialExponents"  -> {2 eps - 1},
        "PolynomialExponents" -> {-1},
        "Variables"          -> {x[1]},
        "KinematicSymbols"   -> {},
        "RegulatorSymbol"    -> eps
      |>;
      verts = PolytopeVertices[(1 + x[1]^2)^(-1), {x[1]}];
      fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
      result = Quiet[
        EvaluateTropicalMC[spec, fanData, {{}},
          "RunChecks" -> False, "Verbose" -> False,
          "NSamples" -> 100],
        TropicalEval::noregulator
      ];
    ];
    pass = (result === $Failed);
    Print["  Part A (eps-regulated spec -> $Failed): ",
          If[pass, "PASS", "FAIL (got " <> ToString[result] <> ")"]];
    If[!pass, allPass = False];
  ];

  (* Part B: regulator-free but genuinely divergent spec -> $Failed with divergentinput *)
  Module[{pass, result, spec, verts, fanData},
    spec = <|
      "Polynomials"        -> {1 + x[1] + x[2]},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-1},
      "Variables"          -> {x[1], x[2]},
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;
    verts = PolytopeVertices[(1 + x[1] + x[2])^(-1), {x[1], x[2]}];
    fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
    result = Quiet[
      EvaluateTropicalMC[spec, fanData, {{}},
        "RunChecks" -> False, "Verbose" -> False,
        "NSamples" -> 100],
      {TropicalEval::divergentinput}
    ];
    pass = (result === $Failed);
    Print["  Part B (genuinely divergent spec -> $Failed): ",
          If[pass, "PASS", "FAIL (got " <> ToString[result] <> ")"]];
    If[!pass, allPass = False];
  ];

  Print["  Test 3v2 ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 5: End-to-end kinematic scan
   -------------------------------------------------------------------------- *)

RunTest5[] := Module[
  {poly, vars, lam, A, spec, verts, fanData,
   lamValues, nPoints,
   allPass, maxRelErr},

  Print["--- Test 5: End-to-end kinematic scan ---"];
  Print["Int dx1 dx2 / (1+lam*x1^2+x2^2+x1*x2^2)^{2+0.5I}"];
  Print["100 values of lam in [0.1, 10]"];

  lam = Symbol["lam"];
  A   = 2 + 0.5 I;
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  nPoints  = 100;
  lamValues = Table[0.1 + (10.0 - 0.1) (i - 1)/(nPoints - 1),
                    {i, nPoints}];

  spec = <|
    "Polynomials"       -> {poly},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> vars,
    "KinematicSymbols"  -> {lam},
    "RegulatorSymbol"   -> None
  |>;

  verts   = PolytopeVertices[(poly /. lam -> 1)^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Computing NIntegrate reference values..."];
  Module[{testIndices, refResults},
    testIndices = {1, 25, 50, 75, 100};
    refResults = Table[
      Module[{lamVal, result},
        lamVal = lamValues[[idx]];
        result = Quiet@NIntegrate[
          1 / (1 + lamVal t1^2 + t2^2 + t1 t2^2)^A,
          {t1, 0, Infinity}, {t2, 0, Infinity},
          MaxRecursion -> 20, PrecisionGoal -> 5
        ];
        {lamVal, result}
      ],
      {idx, testIndices}
    ];

    allPass = True;
    maxRelErr = 0;

    Do[
      Print["  lam = ", refResults[[i, 1]], ": NIntegrate = ",
            refResults[[i, 2]]];,
      {i, Length[refResults]}
    ];

    (* Sector decomposition check *)
    Do[
      Module[{kinRules, vr, relErr},
        kinRules = {lam -> refResults[[i, 1]]};
        vr = Quiet@ValidateDecomposition[spec, fanData, kinRules, 3];
        If[AssociationQ[vr],
          relErr = vr["RelativeError"];
          If[NumericQ[relErr],
            If[relErr > maxRelErr, maxRelErr = relErr];
            If[relErr > 0.05, allPass = False];
            Print["  lam = ", refResults[[i, 1]],
                  ": sector rel err = ", relErr,
                  " ", If[relErr < 0.05, "PASS", "FAIL"]];,
            Print["  lam = ", refResults[[i, 1]],
                  ": non-numeric error, FAIL"];
            allPass = False;
          ];
        ];
      ],
      {i, Length[refResults]}
    ];
  ];

  Print["  Max relative error: ", maxRelErr];
  Print["  ", If[allPass, "PASS", "FAIL"]];

  allPass
];

(* --------------------------------------------------------------------------
   Test 6: Large coefficients
   Verifies that the tropical decomposition and sector integrals remain
   correct when polynomial coefficients span many orders of magnitude.
   The tropically dominant monomial (min exponents) may NOT be the
   numerically largest monomial, but the factoring is an exact algebraic
   identity so the result must still agree with direct NIntegrate.
   -------------------------------------------------------------------------- *)

RunTest6[] := Module[
  {poly, vars, spec, verts, fanData, allPass,
   testCases, t1, t2},

  Print["--- Test 6: Large polynomial coefficients ---"];
  Print["Verifies correctness when coefficients span many orders of magnitude"];

  allPass = True;

  (* Test case A: coefficients O(10^6)
     P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2
     The 10^6 term dominates numerically but is NOT the tropically
     dominant monomial in most sectors. *)
  Module[{polyA, specA, vertsA, fanA, directA, vrA, relErrA},
    Print[];
    Print["  Case A: P = 1 + 10^6 x1^2 + x2^2 + x1*x2^2, exponent -2"];
    polyA = 1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2;
    vars  = {x[1], x[2]};

    specA = <|
      "Polynomials"        -> {polyA},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directA = Quiet@NIntegrate[
      1 / (1 + 10^6 t1^2 + t2^2 + t1 t2^2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    vertsA = PolytopeVertices[polyA^(-2), vars];
    fanA   = ComputeDecomposition[vertsA, "ShowProgress" -> False];

    vrA = Quiet@ValidateDecomposition[specA, fanA, {}, 3];
    If[AssociationQ[vrA],
      relErrA = vrA["RelativeError"];
      Print["    NIntegrate = ", directA];
      Print["    Sector sum = ", vrA["SectorSum"]];
      Print["    Rel error  = ", relErrA];
      If[!NumericQ[relErrA] || relErrA > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  (* Test case B: mixed large and small coefficients
     P = 10^(-4) + 10^4 x1^2 + 10^(-4) x2^2 + 10^4 x1 x2^2 + x1^2 x2
     Coefficients span 8 orders of magnitude. *)
  Module[{polyB, specB, vertsB, fanB, directB, vrB, relErrB},
    Print[];
    Print["  Case B: coefficients from 10^-4 to 10^4, exponent -2"];
    polyB = 10^(-4) + 10^4 x[1]^2 + 10^(-4) x[2]^2 +
            10^4 x[1] x[2]^2 + x[1]^2 x[2];
    vars  = {x[1], x[2]};

    specB = <|
      "Polynomials"        -> {polyB},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directB = Quiet@NIntegrate[
      1 / (10^(-4) + 10^4 t1^2 + 10^(-4) t2^2 +
           10^4 t1 t2^2 + t1^2 t2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    vertsB = PolytopeVertices[polyB^(-2), vars];
    fanB   = ComputeDecomposition[vertsB, "ShowProgress" -> False];

    vrB = Quiet@ValidateDecomposition[specB, fanB, {}, 3];
    If[AssociationQ[vrB],
      relErrB = vrB["RelativeError"];
      Print["    NIntegrate = ", directB];
      Print["    Sector sum = ", vrB["SectorSum"]];
      Print["    Rel error  = ", relErrB];
      If[!NumericQ[relErrB] || relErrB > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  (* Test case C: large coefficient with higher exponent
     P = 1 + 10^8 x1^3 x2 + x2^3, exponent -3
     The 10^8 monomial has degree 4 and large coefficient, ensuring
     it numerically dominates even though the constant term is tropically
     dominant in its cone. *)
  Module[{polyC, specC, vertsC, fanC, directC, vrC, relErrC},
    Print[];
    Print["  Case C: P = 1 + 10^8 x1^3*x2 + x2^3, exponent -3"];
    polyC = 1 + 10^8 x[1]^3 x[2] + x[2]^3;
    vars  = {x[1], x[2]};

    specC = <|
      "Polynomials"        -> {polyC},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-3},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directC = Quiet@NIntegrate[
      1 / (1 + 10^8 t1^3 t2 + t2^3)^3,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 5
    ];

    vertsC = PolytopeVertices[polyC^(-3), vars];
    fanC   = ComputeDecomposition[vertsC, "ShowProgress" -> False];

    vrC = Quiet@ValidateDecomposition[specC, fanC, {}, 2];
    If[AssociationQ[vrC],
      relErrC = vrC["RelativeError"];
      Print["    NIntegrate = ", directC];
      Print["    Sector sum = ", vrC["SectorSum"]];
      Print["    Rel error  = ", relErrC];
      If[!NumericQ[relErrC] || relErrC > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  Print[];
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 7: Higher-dimensional integrals (3D and 4D)
   -------------------------------------------------------------------------- *)

RunTest7[] := Module[
  {allPass = True, vars3, poly3, spec3, verts3, fan3, vr3,
   vars4, poly4, spec4, verts4, fan4, vr4},

  Print["--- Test 7: Higher-dimensional integrals (3D and 4D) ---"];
  Print["Tests that the pipeline works correctly in dimensions > 2"];

  (* 3D: Int dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3 *)
  Print[];
  Print["  Case A (3D): P = 1 + x1^2 + x2^2 + x3^2 + x1*x2*x3, exponent -3"];

  poly3 = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  vars3 = {x[1], x[2], x[3]};

  spec3 = <|
    "Polynomials"        -> {poly3},
    "MonomialExponents"  -> {0, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars3,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts3 = PolytopeVertices[poly3^(-1), vars3];
  fan3   = ComputeDecomposition[verts3, "ShowProgress" -> False];
  Print["    Fan: ", Length[fan3[[1]]], " rays, ", Length[fan3[[2]]], " sectors"];

  vr3 = Quiet@ValidateDecomposition[spec3, fan3, {}, 3];
  If[AssociationQ[vr3],
    Print["    NIntegrate = ", vr3["DirectResult"]];
    Print["    Sector sum = ", vr3["SectorSum"]];
    Print["    Rel error  = ", vr3["RelativeError"]];
    If[!NumericQ[vr3["RelativeError"]] || vr3["RelativeError"] > 0.01,
      Print["    FAIL"];
      allPass = False;,
      Print["    PASS"];
    ];,
    Print["    ValidateDecomposition returned non-association, FAIL"];
    allPass = False;
  ];

  (* 4D: Int dx1 dx2 dx3 dx4 / (1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4)^4 *)
  Print[];
  Print["  Case B (4D): P = 1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4, exponent -4"];

  poly4 = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1] x[2] + x[3] x[4];
  vars4 = {x[1], x[2], x[3], x[4]};

  spec4 = <|
    "Polynomials"        -> {poly4},
    "MonomialExponents"  -> {0, 0, 0, 0},
    "PolynomialExponents" -> {-4},
    "Variables"          -> vars4,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts4 = PolytopeVertices[poly4^(-1), vars4];
  fan4   = ComputeDecomposition[verts4, "ShowProgress" -> False];
  Print["    Fan: ", Length[fan4[[1]]], " rays, ", Length[fan4[[2]]], " sectors"];

  vr4 = Quiet@ValidateDecomposition[spec4, fan4, {}, 2];
  If[AssociationQ[vr4],
    Print["    NIntegrate = ", vr4["DirectResult"]];
    Print["    Sector sum = ", vr4["SectorSum"]];
    Print["    Rel error  = ", vr4["RelativeError"]];
    If[!NumericQ[vr4["RelativeError"]] || vr4["RelativeError"] > 0.01,
      Print["    FAIL"];
      allPass = False;,
      Print["    PASS"];
    ];,
    Print["    ValidateDecomposition returned non-association, FAIL"];
    allPass = False;
  ];

  Print[];
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

