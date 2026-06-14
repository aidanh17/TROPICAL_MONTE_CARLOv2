(* ::Package:: *)

(* ============================================================================
   tropical_eval.wl

   Numerical evaluation of generalized Euler integrals via tropical
   decomposition.  Consumes the fan output of tropical_fan.wl and produces
   C++ Monte-Carlo code.

   Pipeline:
     Module 1 - ProcessSector        (symbolic coordinate transform + flattening)
     Module 2 - [removed in v2: divergence regulation; use TROPICAL_MONTE_CARLO for divergent integrals]
     Module 3 - C++ code generation   (MmaToC + GenerateCppMonteCarlo)
     Module 4 - EvaluateTropicalMC    (driver: fan -> C++ -> results)
     Module 5 - RunAllTests           (validation suite)

   Dependencies: tropical_fan.wl (loaded automatically from same directory)
   ============================================================================ *)

(* Load tropical_fan.wl BEFORE BeginPackage so TropicalFan` is on $ContextPath.
   Set $SkipPolymakeLoad = True before loading this file to skip the polymake
   dependency (e.g. on cluster nodes without polymake installed). *)
If[!TrueQ[$SkipPolymakeLoad],
  Get[FileNameJoin[{DirectoryName[$InputFileName], "tropical_fan.wl"}]];
  BeginPackage["TropicalEval`", {"TropicalFan`"}],
  BeginPackage["TropicalEval`"]
];

(* ---- Public symbols ---- *)

ProcessSector::usage =
  "ProcessSector[integrandSpec, dualVertices, simplex, coneIndex] performs \
monomial change of variables and flattening for one simplicial cone.";

FlattenSector::usage =
  "FlattenSector[clearedPolys, effectiveAVals, prefactorBase] performs the \
divergence check, flattening, and prefactor computation for a set of \
tropically-cleared polynomials.  Returns an Association with keys \
IsDivergent, DivergentVariable, FlattenedPolys (list or None), and \
Prefactor (prefactorBase/(Times@@effectiveAVals) or None).  Works for any \
dimension n = Length[effectiveAVals]; ProcessSector calls it with \
n-dimensional inputs, and the lifted path (ProcessSectorLifted) will call \
it with the n post-delta-resolution effective exponents.";

CheckFlatteningMagnitude::usage =
  "CheckFlatteningMagnitude[sectorData, nSamples] spot-checks the flattened \
integrand magnitude at random points.  When sectorData contains a \
DomainConstraint key (from a lifted sector), feasible points are drawn by \
rejection sampling and the returned association gains a FeasibleFraction key. \
The no-constraint path is byte-identical to the pre-lifting behavior.";

ValidateLiftedDecomposition::usage =
  "ValidateLiftedDecomposition[originalSpec, liftedSpec, liftedFanData, \
liftData, testKinematics, precisionGoal:3] cross-checks the lifted sector sum \
against a direct NIntegrate of originalSpec on [0,Inf)^n.  Sectors via \
ProcessSectorLifted; EmptyDomain sectors contribute 0 and are listed under \
DroppedSectors.  Returns <|DirectResult, SectorSum, RelativeError, \
SectorResults, DroppedSectors|>.";

DetectExtremeCoefficients::usage =
  "DetectExtremeCoefficients[integrandSpec, threshold:1000, opts] scans all \
polynomials in integrandSpec (via ParsePolynomial) and returns a list of \
Associations <|\"PolyIndex\"->j, \"ExponentVector\"->alpha, \
\"Coefficient\"->C, \"Magnitude\"->Abs[C], \"SuggestedK\"->kStar|> for every \
monomial whose NUMERIC coefficient satisfies Abs[C] < 1/threshold or \
Abs[C] > threshold.  SuggestedK automates the anchor z0 = Abs[C]^(1/k): \
choosing the integer k FIXES z0, so picking k IS picking the anchor.  By \
default (\"AnchorRule\"->\"kStar\") it returns the smallest k that pulls z0 \
back inside the non-extreme band [1/threshold, threshold], \
kStar = Max[1, Ceiling[Abs[Log[Abs[C]]]/Log[threshold]]] -- this reproduces \
every hand-picked k in the v2 benchmarks (e.g. k=2 for |C|=10^6, 10^4, \
10^-4).  Options: \"AnchorRule\" -> \"kStar\" (default) | \"Unit\" (legacy \
SuggestedK -> 1) | a function f[mag, threshold]; \"BandEdgeGuard\" -> False \
(default) | True (adds +1 when Abs[C] is within one decade of the band edge, \
where the residual spread z0^(De_p/m_p) is still > 1 at kStar; an opt-in \
accuracy refinement, see Manual sec:anchor).  Symbolic or kinematic-symbol \
coefficients are SKIPPED because their magnitude is unknown at lift time; \
lifting those monomials remains possible by constructing explicit liftRules \
and calling LiftCoefficients directly.";

LiftCoefficients::usage =
  "LiftCoefficients[integrandSpec, liftRules] applies the auxiliary-variable \
lifting of §3.1 to integrandSpec.  Each rule has the form \
<|\"PolyIndex\"->j, \"ExponentVector\"->alpha, \"k\"->k|>.  Returns \
<|\"LiftedSpec\"->..., \"LiftData\"-><|\"z0\" (exact anchor), \
\"AuxVariable\"->x[n+1], \"AuxIndex\"->n+1, \"Rules\"->liftRules, \
\"Residuals\"->{residual coefficients c_i}, \
\"OriginalSpec\"->integrandSpec|>|>.  The primary rule (largest |log|C||) \
sets z0 = Abs[C_primary]^(1/k_primary) exactly; every rule i gets a \
residual c_i = C_i/z0^{k_i} so that the lifted polynomial reduces to the \
original under substitution z -> z0.";

ProcessSectorLifted::usage =
  "ProcessSectorLifted[liftedSpec, dualVertices, simplex, coneIndex, \
liftData, opts] runs the full delta-resolution pipeline (§3.2-§3.5) for one \
sector of a lifted integrand.  Returns a convergent SectorData Association \
augmented with DomainConstraint, LiftData, PivotIndex, ZRow, AugmentedA, \
HasConstantTerm; or <|\"EmptyDomain\"->True, \"ConeIndex\"->...|> when the \
delta root lies outside the unit cube; or \$Failed with a diagnostic message.";


ValidateDecomposition::usage =
  "ValidateDecomposition[integrandSpec, fanData, testKinematics, precisionGoal] \
cross-checks sector sum against direct NIntegrate.";

MmaToC::usage =
  "MmaToC[expr, paramMap] converts a Mathematica expression to a C++ string.";

GenerateCppMonteCarlo::usage =
  "GenerateCppMonteCarlo[convergentSectors, divergentSectors, integrandSpec, \
outputFile] generates self-contained C++ Monte Carlo source.  Option \
\"Integrator\" -> \"MC\" (default) emits the plain Monte-Carlo sampler; \
\"Integrator\" -> \"VEGAS\" emits a CUBA-Vegas sampler (adaptive importance \
sampling; requires the CUBA library at compile time).  The \"VEGAS\" branch \
keeps the per-kinematic-point 4-number output contract (re im reErr imErr, \
summed over sectors) and uses ncomp = 2 (Re/Im) so complex integrands are \
first-class.  VEGAS tuning options (used only when Integrator -> \"VEGAS\"): \
\"VegasEpsRel\" -> Automatic | number, \"VegasSeed\" -> 0 (0 = Sobol \
low-discrepancy), \"VegasNStart\" -> 1000, \"VegasNIncrease\" -> 500, \
\"VegasNBatch\" -> 1000, \"VegasMinEval\" -> 0.  The runtime n_samples argv is \
interpreted as maxeval per sector.  Returns an Association whose keys include \
\"Integrator\" and \"NeedsCuba\".";

EvaluateTropicalMC::usage =
  "EvaluateTropicalMC[integrandSpec, fanData, kinematicPoints, opts] runs \
the full tropical Monte Carlo pipeline.  When the option \"LiftData\" is set \
to a lift-data association (from LiftCoefficients), the first argument is \
treated as the LIFTED spec; the fan must be the (n+1)-dimensional lifted fan, \
and all sector processing is routed through ProcessSectorLifted.  Default \
\"LiftData\" -> None leaves all behavior identical to the pre-lifting code.  \
Option \"Integrator\" -> \"MC\" (default) | \"VEGAS\" selects the numerical \
integrator: \"VEGAS\" routes the flattened sectors through CUBA-Vegas \
(requires CUBA; if not found, EvaluateTropicalMC issues TropicalEval::nocuba \
and returns $Failed -- there is no silent fallback to MC).  VEGAS tuning: \
\"VegasEpsRel\" -> Automatic (= 10^-PrecisionGoal, clamped to [1e-12, 1e-2]) | \
number, \"VegasSeed\", \"VegasNStart\", \"VegasNIncrease\", \"VegasNBatch\" \
(all Automatic by default -- resolved from the max sector dimension d: the \
historical 1000/500/1000 for d<=4, scaled up ~1000*3^(d-4) for d>=5 so that \
high-dimensional VEGAS does not silently converge to a wrong value on a too-coarse \
grid; explicit integers are honored), \"VegasMinEval\".  In high dimension also \
raise \"NSamples\" (Vegas maxeval/sector) accordingly, or TropicalEval::vegasbudget \
warns.  With Integrator -> \"MC\" the generated source, compile flags, and results \
are byte-identical to the pre-VEGAS pipeline.";

EvaluateTropicalMCLifted::usage =
  "EvaluateTropicalMCLifted[integrandSpec, kinematicPoints, opts] is a \
convenience wrapper that automatically detects extreme coefficients, lifts \
the integrand, builds the lifted fan, and calls EvaluateTropicalMC with \
\"LiftData\" set appropriately.  The \"Integrator\" and \"Vegas*\" options are \
passed through to EvaluateTropicalMC (so Integrator -> \"VEGAS\" works on \
lifted sectors too).  Options: \"LiftRules\" -> Automatic | \
explicit list (list of <|\"PolyIndex\"->j, \"ExponentVector\"->alpha, \
\"k\"->k|> associations); \"Threshold\" -> 1000, \"AnchorRule\" -> \"kStar\", \
\"BandEdgeGuard\" -> False (all passed to DetectExtremeCoefficients when \
\"LiftRules\" -> Automatic, automating the anchor z0 = |C|^(1/k) -- see \
DetectExtremeCoefficients); \
\"FanData\" -> Automatic | explicit {dualVertices, simplexList} for the \
LIFTED fan.  All other options are passed through to EvaluateTropicalMC. \
If Automatic detection finds no extreme coefficients, falls back to plain \
EvaluateTropicalMC on the original spec.";

RunAllTests::usage =
  "RunAllTests[] runs the validation suite (6 tests: 1, 2, 3v2, 5, 6, 7) with structured reporting. \
Defined in tropical_eval_examples.wl.";

ParsePolynomial::usage =
  "ParsePolynomial[poly, vars] parses a polynomial into a list of \
{coefficient, exponentVector} pairs.";

CompileCpp::usage =
  "CompileCpp[srcFile, outputBinary, debug] compiles generated C++ \
Monte Carlo source code. debug=True adds -DTROPICAL_MC_DEBUG.";

(* ---- Error messages ---- *)

TropicalEval::degenerate = "Sector `1`: degenerate cone, det(M) = 0.";
TropicalEval::notsimplicial = "Sector `1`: `2` rays for `3` variables (not simplicial).";
TropicalEval::divergent = "Sector `1`: variable y_`2` divergent, a_`2` = `3`.";
TropicalEval::nested = "Sector `1`: multiple divergent variables (`2`). Nested subtraction not implemented.";
TropicalEval::badck = "Sector `1`: c_k = 0 for divergent variable y_`2`. Higher-order pole.";
TropicalEval::validate = "Validation `1`: relative error `2` exceeds tolerance `3`.";
TropicalEval::divergentinput = "Sectors `1` have divergent exponents `2`; TROPICAL_MONTE_CARLOv2 supports convergent integrals only — use the original TROPICAL_MONTE_CARLO package for eps-regulated/divergent integrals.";
TropicalEval::noregulator = "RegulatorSymbol must be None (or absent) in TROPICAL_MONTE_CARLOv2; this version supports convergent integrals only — use the original TROPICAL_MONTE_CARLO package for eps-regulated/divergent integrals.";
TropicalEval::nodivergent = "GenerateCppMonteCarlo: divergent sectors passed; not supported in TROPICAL_MONTE_CARLOv2.";
TropicalEval::liftidentity = "LiftCoefficients: round-trip identity check FAILED for polynomial `1`; lifted poly at z->z0 does not match original.";
TropicalEval::liftnopivot = "ProcessSectorLifted: cone `1` — no admissible pivot found.  z-row m=`2`, per-pivot atilde=`3`.  Try a different k in the lift rules; alternatively, the domain constraint may cut off all divergent regions (log-space remap, future work).";
TropicalEval::liftcomplex = "ProcessSectorLifted: cone `1` — all candidate pivots produce complex atilde; cannot emit a real-valued domain indicator.  Lift with a different k or check that polynomial exponents B are real.";
TropicalEval::liftcomplexexponents =
  "EvaluateTropicalMCLifted: `1` of `2` polynomial exponents (B_k / \[Gamma]_k) are \
complex.  The current lifting algorithm requires real exponents to construct a \
real-valued domain indicator.  Options: (a) use \"ComplexExponentMode\"->\"SplitRealImag\" \
to fold imaginary parts into an oscillatory phase weight; \
(b) condition the exponent representation before lifting.";
TropicalEval::liftdivergent = "ProcessSectorLifted: cone `1` — atilde `2` has a non-positive component after delta resolution; the lifted sector is divergent.  TROPICAL_MONTE_CARLOv2 supports convergent integrals only.";
TropicalEval::liftfandim = "EvaluateTropicalMC with LiftData: the fan dimension is `1` (Length[dualVertices[[1]]]) but n+1 = `2` is required (n = Length[liftData[\"OriginalSpec\"][\"Variables\"]]).  Supply the (n+1)-dimensional lifted fan.";
TropicalEval::liftdegenerate = "EvaluateTropicalMCLifted: the lifted Newton polytope is lower-dimensional; automatic fan construction is not possible — supply an explicit complete simplicial fan via the \"FanData\" option.";
TropicalEval::nocuba = "Integrator -> \"VEGAS\" requires the CUBA library (cuba.h + libcuba). Not found under /opt/homebrew, /usr/local, or /usr. Install via `brew install cuba` or https://feynarts.de/cuba/, or use Integrator -> \"MC\".";
TropicalEval::badintegrator = "Unknown Integrator `1`; expected \"MC\" or \"VEGAS\".";
TropicalEval::vegasbudget = "Integrator -> \"VEGAS\" in dimension `1`: the maxeval budget NSamples = `2` is small relative to the (dimension-aware) per-iteration grid VegasNStart = `3`.  VEGAS may stop before its grid has adapted and return a confidently wrong value with a deceptively tight error bar.  Increase \"NSamples\" to at least ~`4`, or set \"VegasNStart\" explicitly.";

(* ============================================================================
   PRIVATE IMPLEMENTATION
   ============================================================================ *)

Begin["`Private`"]

(* tropical_fan.wl is loaded before BeginPackage above *)

(* --------------------------------------------------------------------------
   Helper: ParsePolynomial
   Converts a polynomial into {coefficient, exponentVector} pairs.
   -------------------------------------------------------------------------- *)

ParsePolynomial[poly_, vars_List] := Module[
  {expanded, terms, result},
  expanded = Expand[poly];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  result = Table[
    Module[{coeff, exps},
      exps = Exponent[term, #] & /@ vars;
      coeff = term / (Times @@ MapThread[Power, {vars, exps}]);
      {Simplify[coeff], exps}
    ],
    {term, terms}
  ];
  result
];

(* --------------------------------------------------------------------------
   Helper: TransformExponents
   Given original exponent vector and matrix M, compute new exponent vector.
   -------------------------------------------------------------------------- *)

TransformExponents[expVec_List, mMatrix_List] := expVec . mMatrix;

(* --------------------------------------------------------------------------
   FlattenSector
   Divergence check + flattening + prefactor computation.
   Extracted from ProcessSector so that the lifted path (ProcessSectorLifted)
   can reuse it with n post-delta-resolution effective exponents.
   -------------------------------------------------------------------------- *)

FlattenSector[clearedPolys_List, effectiveAVals_List, prefactorBase_] :=
Module[
  {n, isDivergent, divVar, a0vals, flattenedPolys, prefactor},

  n = Length[effectiveAVals];
  isDivergent = False;
  divVar = 0;

  a0vals = effectiveAVals;
  Do[
    If[TrueQ[Re[a0vals[[i]]] <= 0] ||
       (NumericQ[a0vals[[i]]] && Re[a0vals[[i]]] <= 0),
      isDivergent = True;
      divVar = i;
    ],
    {i, n}
  ];

  If[isDivergent,
    Return[<|
      "IsDivergent"      -> True,
      "DivergentVariable" -> divVar,
      "FlattenedPolys"   -> None,
      "Prefactor"        -> None
    |>]
  ];

  (* Convergent: flatten y_i -> (y_i')^{1/a_i^eff} *)
  flattenedPolys = Table[
    Table[
      {mono[[1]],
       MapThread[#1/#2 &, {mono[[2]], effectiveAVals}]},
      {mono, clearedPolys[[j]]}
    ],
    {j, Length[clearedPolys]}
  ];

  prefactor = prefactorBase / (Times @@ effectiveAVals);

  <|
    "IsDivergent"       -> False,
    "DivergentVariable" -> 0,
    "FlattenedPolys"    -> flattenedPolys,
    "Prefactor"         -> prefactor
  |>
];

(* --------------------------------------------------------------------------
   MODULE 1: ProcessSector

   Key insight (tropical factoring):
   After the monomial substitution x_i = prod y_j^{M_ij}, each polynomial
   P_k becomes a sum of monomials in y with SIGNED exponents.  The dominant
   monomial (the one the tropical fan says dominates in this cone) has the
   minimum exponents.  We factor it out:
       P_k(y) = prod_j y_j^{d_{k,j}} * Q_k(y)
   where d_{k,j} = min_m (transformed exponent of y_j in monomial m),
   and Q_k has all non-negative y-exponents with a constant term.

   The effective monomial prefactor then becomes:
       a_j^eff = rawA_j + sum_k B_k * d_{k,j}

   For a properly constructed tropical fan and convergent integral,
   a_j^eff > 0 for all j, and we can flatten using these effective exponents.
   -------------------------------------------------------------------------- *)

Options[ProcessSector] = {"Verbose" -> False};

ProcessSector[integrandSpec_Association, dualVertices_List,
              simplex_List, coneIndex_Integer, OptionsPattern[]] :=
Module[
  {polys, monoExps, polyExps, vars, kinSyms,
   selectedRays, mMatrix, detM, n,
   transformedPolys, clearedPolys, minExponents,
   rawAVals, effectiveAVals,
   flattenedPolys, prefactor,
   isDivergent, divVar, verbose, sectorData,
   parsedPolys},

  If[!MatchQ[integrandSpec["RegulatorSymbol"], None | _Missing],
    Message[TropicalEval::noregulator];
    Return[$Failed]
  ];

  verbose = OptionValue["Verbose"];

  (* Extract fields from integrand spec *)
  polys    = integrandSpec["Polynomials"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  vars     = integrandSpec["Variables"];
  kinSyms  = integrandSpec["KinematicSymbols"];
  n        = Length[vars];

  (* --- Step 1: Monomial change of variables --- *)

  (* Extract ray vectors for this simplex (rows of dualVertices) *)
  selectedRays = dualVertices[[#]] & /@ simplex;

  (* Check: simplicial condition *)
  If[Length[selectedRays] != n,
    Message[TropicalEval::notsimplicial, coneIndex,
            Length[selectedRays], n];
    Return[$Failed]
  ];

  (* M_{ij} = -rho_j[i], i.e. M = -Transpose[selectedRays] *)
  mMatrix = -Transpose[selectedRays];

  (* Check: non-degenerate *)
  detM = Det[mMatrix];
  If[detM === 0 || TrueQ[detM == 0],
    Message[TropicalEval::degenerate, coneIndex];
    Return[$Failed]
  ];

  (* Raw exponents from monomial part + Jacobian:
     rawA_i = sum_k (A_k + 1) * M_{ki} = (monoExps + 1) . M *)
  rawAVals = (monoExps + 1) . mMatrix;

  (* --- Transform polynomials --- *)
  parsedPolys = ParsePolynomial[#, vars] & /@ polys;

  transformedPolys = Table[
    Table[
      Module[{coeff, origExp, newExp},
        coeff   = mono[[1]];
        origExp = mono[[2]];
        newExp  = TransformExponents[origExp, mMatrix];
        {coeff, newExp}
      ],
      {mono, parsedPolys[[j]]}
    ],
    {j, Length[polys]}
  ];

  (* --- Step 1b: Tropical factoring --- *)
  (* For each polynomial, find min exponents and factor them out *)

  minExponents = Table[
    Table[
      Min[#[[2, i]] & /@ transformedPolys[[j]]],
      {i, n}
    ],
    {j, Length[polys]}
  ];

  (* Cleared polynomials: shift exponents so minimum is 0 *)
  clearedPolys = Table[
    Table[
      {mono[[1]], mono[[2]] - minExponents[[j]]},
      {mono, transformedPolys[[j]]}
    ],
    {j, Length[polys]}
  ];

  (* Effective exponents: a_i^eff = rawA_i + sum_j B_j * d_{j,i} *)
  effectiveAVals = rawAVals + Total[
    Table[polyExps[[j]] * minExponents[[j]], {j, Length[polys]}]
  ];

  If[verbose,
    Print["Sector ", coneIndex, ": det(M) = ", detM,
          ", rawA = ", rawAVals,
          ", minExp = ", minExponents,
          ", effA = ", effectiveAVals]
  ];

  (* --- Step 2: Flattening using effective exponents --- *)
  (* Delegate divergence check + flattening + prefactor to FlattenSector *)
  Module[{fsResult},
    fsResult = FlattenSector[clearedPolys, effectiveAVals, Abs[detM]];
    isDivergent = fsResult["IsDivergent"];
    divVar      = fsResult["DivergentVariable"];

    If[isDivergent,
      (* Divergent: return pre-flattened data for the lifted path *)
      sectorData = <|
        "ConeIndex"           -> coneIndex,
        "RayMatrix"           -> mMatrix,
        "DetM"                -> detM,
        "SelectedRays"        -> selectedRays,
        "RawExponents"        -> rawAVals,
        "NewExponents"        -> effectiveAVals,
        "MinExponents"        -> minExponents,
        "TransformedPolys"    -> transformedPolys,
        "ClearedPolys"        -> clearedPolys,
        "Prefactor"           -> Abs[detM],
        "IsDivergent"         -> True,
        "DivergentVariable"   -> divVar,
        "Dimension"           -> n,
        "PolynomialExponents" -> polyExps,
        "MonomialExponents"   -> monoExps
      |>;
      If[verbose,
        Print["  -> Divergent in variable y_", divVar]
      ];
      Return[sectorData]
    ];

    flattenedPolys = fsResult["FlattenedPolys"];
    prefactor      = fsResult["Prefactor"];
  ];

  sectorData = <|
    "ConeIndex"            -> coneIndex,
    "RayMatrix"            -> mMatrix,
    "DetM"                 -> detM,
    "SelectedRays"         -> selectedRays,
    "RawExponents"         -> rawAVals,
    "NewExponents"         -> effectiveAVals,
    "MinExponents"         -> minExponents,
    "TransformedPolys"     -> transformedPolys,
    "ClearedPolys"         -> clearedPolys,
    "FlattenedPolys"       -> flattenedPolys,
    "Prefactor"            -> prefactor,
    "IsDivergent"          -> False,
    "DivergentVariable"    -> 0,
    "Dimension"            -> n,
    "PolynomialExponents"  -> polyExps,
    "MonomialExponents"    -> monoExps,
    (* Per-polynomial log of the tropical monomial factor y^{d_k} that was
       cleared out of P_k, expressed in the FLATTENED sample coords y'
       (y_j = y'_j^{1/a_eff,j} => log y_j = log_y'[j]/a_eff,j).  Needed by the
       SplitRealImag oscillatory phase exp(i Im(B_k) log|P_k|), since
       log|P_k| = Const_k + sum_i Coeffs_{k,i} log_y'[i] + log|Q_k|, and the
       cleared-polynomial value the integrand evaluates is only |Q_k|. *)
    "MonoFactorLog"        -> Table[
      <|"Const" -> 0,
        "Coeffs" -> Table[minExponents[[k, i]]/effectiveAVals[[i]], {i, n}]|>,
      {k, Length[minExponents]}]
  |>;

  If[verbose,
    Print["  -> Convergent, prefactor = ", prefactor]
  ];

  sectorData
];

(* --------------------------------------------------------------------------
   CheckFlatteningMagnitude
   Spot-check that the flattened integrand is O(1) at random points.
   -------------------------------------------------------------------------- *)

CheckFlatteningMagnitude[sectorData_Association, nSamples_Integer: 20,
                         testKinematics_List: {}] :=
Module[
  {flatPolys, polyExps, prefactor, dim, mags, y, polyVals, integrandVal,
   kinRules, dc},

  If[sectorData["IsDivergent"],
    Print["CheckFlatteningMagnitude: sector ", sectorData["ConeIndex"],
          " is divergent, skipping."];
    Return[Null]
  ];

  flatPolys = sectorData["FlattenedPolys"];
  polyExps  = sectorData["PolynomialExponents"];
  prefactor = sectorData["Prefactor"];
  dim       = sectorData["Dimension"];
  kinRules  = If[testKinematics === {}, {}, testKinematics];
  dc        = Lookup[sectorData, "DomainConstraint", None];

  (* B6: lifted-sector path — rejection sampling for domain-constrained sectors *)
  If[dc =!= None,
    Module[{logZ0num, mpNum, icNum, feasibleMags, totalDraws, feasibleCount,
            maxDraws, isFeasible, logYpStar, mag},

      logZ0num = N[dc["LogZ0"]];
      mpNum    = N[dc["MP"]];
      icNum    = N[dc["IndicatorCoeffs"]];
      maxDraws = 50 * nSamples;
      feasibleMags  = {};
      totalDraws    = 0;
      feasibleCount = 0;

      While[feasibleCount < nSamples && totalDraws < maxDraws,
        y = RandomReal[{0.01, 0.99}, dim];
        totalDraws++;

        (* Check domain constraint *)
        logYpStar = (logZ0num - Total[icNum * Log[y]]) / mpNum;
        isFeasible = (logYpStar <= 0);

        If[isFeasible,
          polyVals = Table[
            Total[
              Table[
                Module[{coeff, alphas, logY2},
                  coeff  = mono[[1]] /. kinRules;
                  alphas = mono[[2]] /. kinRules;
                  logY2  = Log[y];
                  coeff * Exp[Total[alphas * logY2]]
                ],
                {mono, flatPolys[[j]]}
              ]
            ],
            {j, Length[flatPolys]}
          ];
          integrandVal = (prefactor /. kinRules) *
            Times @@ MapThread[
              Exp[#2 * Log[#1]] &,
              {polyVals, polyExps /. kinRules}
            ];
          mag = Abs[integrandVal];
          AppendTo[feasibleMags, mag];
          feasibleCount++
        ]
      ];

      If[feasibleCount == 0,
        Print["WARNING: CheckFlatteningMagnitude sector ",
              sectorData["ConeIndex"],
              ": ZERO feasible points found in ", totalDraws, " draws."];
        Return[<|"Mean" -> 0, "Max" -> 0, "Min" -> 0,
                 "Samples" -> {},
                 "FeasibleFraction" -> 0|>]
      ];

      If[feasibleCount < nSamples / 10,
        Print["WARNING: CheckFlatteningMagnitude sector ",
              sectorData["ConeIndex"],
              ": only ", feasibleCount, " feasible points in ",
              totalDraws, " draws (fraction ",
              N[feasibleCount / totalDraws], "); continuing with what was found."]
      ];

      Module[{meanMag, maxMag, minMag, ff},
        meanMag = Mean[feasibleMags];
        maxMag  = Max[feasibleMags];
        minMag  = Min[feasibleMags];
        ff      = N[feasibleCount / totalDraws];
        If[maxMag > 10^3 || minMag < 10^(-6),
          Print["WARNING: Sector ", sectorData["ConeIndex"],
                " flattening check: min=", minMag, " max=", maxMag,
                " mean=", meanMag, " feasibleFrac=", ff]
        ];
        <|"Mean" -> meanMag, "Max" -> maxMag, "Min" -> minMag,
          "Samples" -> feasibleMags,
          "FeasibleFraction" -> ff|>
      ]
    ]
    ,
    (* No DomainConstraint: original code path, byte-identical *)
    mags = Table[
      y = RandomReal[{0.01, 0.99}, dim];
      polyVals = Table[
        Total[
          Table[
            Module[{coeff, alphas, logY},
              coeff  = mono[[1]] /. kinRules;
              alphas = mono[[2]] /. kinRules;
              logY   = Log[y];
              coeff * Exp[Total[alphas * logY]]
            ],
            {mono, flatPolys[[j]]}
          ]
        ],
        {j, Length[flatPolys]}
      ];
      integrandVal = (prefactor /. kinRules) *
        Times @@ MapThread[
          Exp[#2 * Log[#1]] &,
          {polyVals, polyExps /. kinRules}
        ];
      Abs[integrandVal],
      {nSamples}
    ];

    Module[{meanMag, maxMag, minMag},
      meanMag = Mean[mags];
      maxMag  = Max[mags];
      minMag  = Min[mags];
      If[maxMag > 10^3 || minMag < 10^(-6),
        Print["WARNING: Sector ", sectorData["ConeIndex"],
              " flattening check: min=", minMag, " max=", maxMag,
              " mean=", meanMag]
      ];
      <|"Mean" -> meanMag, "Max" -> maxMag, "Min" -> minMag,
        "Samples" -> mags|>
    ]
  ]
];

(* ============================================================================
   MODULE 1b: DETECTION AND LIFTING  (B2)
   ============================================================================ *)

(* --------------------------------------------------------------------------
   DetectExtremeCoefficients
   Scan every polynomial in integrandSpec for numeric coefficients outside
   [1/threshold, threshold].  Returns a list of flagged-monomial associations.
   Symbolic/kinematic coefficients are silently skipped.
   -------------------------------------------------------------------------- *)

Options[DetectExtremeCoefficients] = {
  "AnchorRule"    -> "kStar",
  "BandEdgeGuard" -> False
};

(* Threshold defaults to 1000.  Two definitions (rather than threshold_:1000)
   so the optional numeric threshold is never confused with a trailing option
   rule: DetectExtremeCoefficients[spec, "AnchorRule"->...] routes here. *)
DetectExtremeCoefficients[integrandSpec_Association, opts : OptionsPattern[]] :=
  DetectExtremeCoefficients[integrandSpec, 1000, opts];

DetectExtremeCoefficients[integrandSpec_Association, threshold_?NumericQ,
                          opts : OptionsPattern[]] :=
Module[
  {polys, vars, result, parsedPoly, coeff, mag,
   anchorRule, bandEdgeGuard, logTau, suggestK},

  polys  = integrandSpec["Polynomials"];
  vars   = integrandSpec["Variables"];
  result = {};

  anchorRule    = OptionValue["AnchorRule"];
  bandEdgeGuard = TrueQ[OptionValue["BandEdgeGuard"]];
  logTau        = Log[N[threshold]];

  (* SuggestedK automates the anchor z0 = |C|^(1/k).  The anchor is forced to
     |C|^(1/k), so choosing the integer k IS choosing z0.  The "kStar" rule
     picks the smallest k that pulls z0 back inside the detector's own
     non-extreme band [1/threshold, threshold]:
        kStar = max(1, ceil(|log|C|| / log threshold)).
     This minimizes the residual coefficient spread left inside the sampled
     integrand while keeping the lifted geometry well-conditioned, and it
     reproduces every hand-tuned k in the v2 benchmarks (see Manual
     sec:anchor and AUXT/z0_sweep_results.md).  The optional band-edge guard
     adds +1 when |C| sits within one decade of the band edge, where the
     residual spread is still > 1 at kStar (the moderate-coefficient case).
     "Unit" recovers the legacy SuggestedK -> 1; a function is called as
     anchorRule[mag, threshold]. *)
  suggestK[m_] := Switch[anchorRule,
    "Unit",  1,
    "kStar",
      Module[{absLog = Abs[Log[N[m]]], kS},
        kS = Max[1, Ceiling[absLog / logTau]];
        If[bandEdgeGuard && absLog <= logTau + Log[10.] (1 + 1.*^-8),
          kS = kS + 1
        ];
        kS
      ],
    _, anchorRule[m, threshold]
  ];

  Do[
    parsedPoly = ParsePolynomial[polys[[j]], vars];
    Do[
      coeff = mono[[1]];
      (* Skip symbolic / kinematic coefficients *)
      If[NumericQ[N[coeff]],
        mag = Abs[N[coeff]];
        If[mag < 1/threshold || mag > threshold,
          AppendTo[result, <|
            "PolyIndex"     -> j,
            "ExponentVector"-> mono[[2]],
            "Coefficient"   -> coeff,
            "Magnitude"     -> mag,
            "SuggestedK"    -> suggestK[mag]
          |>]
        ]
      ],
      {mono, parsedPoly}
    ],
    {j, Length[polys]}
  ];

  result
];


(* --------------------------------------------------------------------------
   LiftCoefficients
   Apply auxiliary-variable lifting per §3.1.
   liftRules = { <|"PolyIndex"->j, "ExponentVector"->alpha, "k"->k|>, ... }
   Returns <|"LiftedSpec"->..., "LiftData"->...|>.
   -------------------------------------------------------------------------- *)

LiftCoefficients[integrandSpec_Association, liftRules_List] :=
Module[
  {polys, vars, monoExps, polyExps, kinSyms,
   n, auxIdx, auxVar, newVars, newMonoExps,
   parsedPolys,
   (* primary rule selection *)
   ruleCoeffs, ruleLogMags, primaryIdx, primaryRule,
   Cprimary, kprimary, z0,
   (* residuals *)
   residuals,
   (* lifted polynomials *)
   liftedPolys, j, polyParsed,
   liftedSpec, liftData,
   (* identity check *)
   liftedSubbed, original},

  polys    = integrandSpec["Polynomials"];
  vars     = integrandSpec["Variables"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  kinSyms  = integrandSpec["KinematicSymbols"];
  n        = Length[vars];
  auxIdx   = n + 1;
  auxVar   = Head[vars[[1]]][auxIdx];   (* e.g. x[n+1] *)

  (* ---- Find the primary rule (max |Log[|C|]|) ---- *)
  parsedPolys = ParsePolynomial[#, vars] & /@ polys;

  (* Look up the coefficient for each rule *)
  ruleCoeffs = Table[
    Module[{j0 = r["PolyIndex"], alpha = r["ExponentVector"], matchPos, C0},
      matchPos = Position[parsedPolys[[j0, All, 2]], alpha];
      If[matchPos === {},
        Print["LiftCoefficients: monomial ", alpha,
              " not found in polynomial ", j0, "."];
        Return[$Failed]
      ];
      parsedPolys[[j0, matchPos[[1, 1]], 1]]
    ],
    {r, liftRules}
  ];
  If[MemberQ[ruleCoeffs, $Failed], Return[$Failed]];

  ruleLogMags = Abs[Log[Abs[N[#]]]] & /@ ruleCoeffs;
  primaryIdx  = First@Ordering[ruleLogMags, -1];
  primaryRule = liftRules[[primaryIdx]];
  Cprimary    = ruleCoeffs[[primaryIdx]];
  kprimary    = primaryRule["k"];

  (* z0 = Abs[C_primary]^(1/k_primary), exact.  Sign goes into residual. *)
  z0 = If[IntegerQ[Abs[Cprimary]^(1/kprimary)],
    Abs[Cprimary]^(1/kprimary),
    Power[Abs[Cprimary], 1/kprimary]
  ];

  (* ---- Compute residuals: c_i = C_i / z0^{k_i} ---- *)
  residuals = Table[
    Simplify[ruleCoeffs[[i]] / z0^liftRules[[i]]["k"]],
    {i, Length[liftRules]}
  ];

  (* ---- Build lifted polynomials ---- *)
  (* For each polynomial, replace the affected monomials *)
  liftedPolys = Table[
    Module[{poly = polys[[polyJ]], acc = polys[[polyJ]]},
      (* Apply each lift rule that targets this polynomial *)
      Do[
        Module[{r = liftRules[[ri]], alpha, ki, Ci, ci,
                termToReplace, replacement},
          If[r["PolyIndex"] == polyJ,
            alpha = r["ExponentVector"];
            ki    = r["k"];
            ci    = residuals[[ri]];
            (* Find C_i from parsedPolys *)
            Ci    = ruleCoeffs[[ri]];
            (* C_i * prod x[j]^alpha_j *)
            termToReplace = Ci * Times @@ MapThread[
              Function[{v, e}, If[e == 0, 1, Power[v, e]]],
              {vars, alpha}
            ];
            (* c_i * auxVar^ki * prod x[j]^alpha_j *)
            replacement = ci * auxVar^ki * Times @@ MapThread[
              Function[{v, e}, If[e == 0, 1, Power[v, e]]],
              {vars, alpha}
            ];
            acc = Expand[acc - termToReplace + replacement];
          ]
        ],
        {ri, Length[liftRules]}
      ];
      acc
    ],
    {polyJ, Length[polys]}
  ];

  newVars     = Append[vars, auxVar];
  newMonoExps = Append[monoExps, 0];

  liftedSpec = <|
    "Polynomials"         -> liftedPolys,
    "MonomialExponents"   -> newMonoExps,
    "PolynomialExponents" -> polyExps,
    "Variables"           -> newVars,
    "KinematicSymbols"    -> kinSyms,
    "RegulatorSymbol"     -> Lookup[integrandSpec, "RegulatorSymbol", None]
  |>;

  (* ---- Identity check: liftedPolys /. auxVar -> z0 === originalPolys ---- *)
  (* Coefficient-wise relative magnitude check avoids false positives when z0 is
     a machine float (e.g. from a real extreme coefficient): Simplify cannot
     reduce the O(eps_mach) residuals to the symbol 0, but they are harmless. *)
  Do[
    liftedSubbed = Expand[liftedPolys[[j]] /. auxVar -> z0];
    original     = Expand[polys[[j]]];
    If[!TrueQ[liftedSubbed === original],
      Module[{parsedDiff, scaleOrig, maxResid},
        parsedDiff = ParsePolynomial[Expand[liftedSubbed - original], vars];
        scaleOrig  = Max[1., Max[Abs[N[#[[1]]]] & /@ ParsePolynomial[Expand[original], vars]]];
        maxResid   = Max[Abs[N[#[[1]]]] & /@ parsedDiff];
        If[NumericQ[maxResid] && maxResid > 1*^-8 * scaleOrig,
          Message[TropicalEval::liftidentity, j];
          Print["  liftedSubbed = ", liftedSubbed];
          Print["  original     = ", original];
          (* Do not abort — residuals allow exact cancellation; try to continue *)
        ]
      ]
    ],
    {j, Length[polys]}
  ];

  liftData = <|
    "z0"           -> z0,
    "AuxVariable"  -> auxVar,
    "AuxIndex"     -> auxIdx,
    "Rules"        -> liftRules,
    "Residuals"    -> residuals,
    "OriginalSpec" -> integrandSpec,
    "ImagExponents"-> Lookup[integrandSpec, "ImagExponents", None]
  |>;

  <|"LiftedSpec" -> liftedSpec, "LiftData" -> liftData|>
];


(* ============================================================================
   MODULE 1c: ProcessSectorLifted  (B3)
   ============================================================================ *)

(* --------------------------------------------------------------------------
   ProcessSectorLifted
   Delta-resolution pipeline (§3.2–§3.5) for one sector of a lifted integrand.
   -------------------------------------------------------------------------- *)

Options[ProcessSectorLifted] = {"Verbose" -> False};

ProcessSectorLifted[liftedSpec_Association, dualVertices_List,
                   simplex_List, coneIndex_Integer,
                   liftData_Association, OptionsPattern[]] :=
Module[
  {sdAug, z0, auxIdx, a, clearedPolys, detM, mMatrix, polyExps, n1, n,
   mVec, verbose,
   tryPivot, candidates, bestPivot,
   (* pivot result fields *)
   pivotP, mp, ap, remainIdx, mOtherVec, atilde,
   reclearedPolys, atildeVals,
   (* domain constraint *)
   allOtherZero, domainClass, logZ0,
   (* flatten & assemble *)
   fsResult, flattenedPolys, prefactor,
   prefactorBase},

  verbose = OptionValue["Verbose"];

  (* --- Step 1: run standard (n+1)-dimensional ProcessSector --- *)
  sdAug = ProcessSector[liftedSpec, dualVertices, simplex, coneIndex];
  If[sdAug === $Failed, Return[$Failed]];

  (* Extract from sdAug — both IsDivergent branches expose these keys *)
  z0        = liftData["z0"];
  auxIdx    = liftData["AuxIndex"];
  a            = sdAug["NewExponents"];     (* length n+1 *)
  clearedPolys = sdAug["ClearedPolys"];
  detM         = sdAug["DetM"];
  mMatrix      = sdAug["RayMatrix"];
  polyExps     = sdAug["PolynomialExponents"];

  n1 = Length[a];   (* n+1 *)
  n  = n1 - 1;      (* original n *)

  (* z-row: row auxIdx of M *)
  mVec = mMatrix[[auxIdx]];   (* length n+1 *)

  (* ----------------------------------------------------------------
     tryPivot[p]: compute monomial substitution + re-clear + atilde.
     Returns an Association or $Failed.
     ---------------------------------------------------------------- *)
  tryPivot[p_] := Module[
    {mpLocal, apLocal, rIdx, mOther, aOther,
     subPolys, rcMin, newPolyList, atildeRaw, hasConst},

    mpLocal = mVec[[p]];
    If[mpLocal == 0, Return[$Failed]];
    apLocal = a[[p]];
    rIdx    = DeleteCases[Range[n1], p];   (* indices of the n remaining vars *)
    mOther  = mVec[[rIdx]];               (* m_j for j != p *)
    aOther  = a[[rIdx]];                  (* a_j for j != p *)

    (* §3.2 step 4: monomial substitution *)
    subPolys = Table[
      Map[Function[{cmono},
        Module[{ep = cmono[[2, p]]},
          {cmono[[1]] * z0^(ep / mpLocal),
           Table[cmono[[2, rIdx[[jj]]]] - ep * mOther[[jj]] / mpLocal, {jj, n}]}
        ]
      ], clearedPolys[[k]]],
      {k, Length[clearedPolys]}
    ];

    (* §3.3 re-clear: componentwise minima *)
    rcMin = Table[
      Table[Min[#[[2, jj]] & /@ subPolys[[k]]], {jj, n}],
      {k, Length[subPolys]}
    ];

    (* subtract minima *)
    newPolyList = Table[
      Map[Function[{cmono}, {cmono[[1]], cmono[[2]] - rcMin[[k]]}],
          subPolys[[k]]],
      {k, Length[subPolys]}
    ];

    (* §3.2 step 5: atilde_j = a_j - a_p * m_j/m_p + Sum_k B_k * dtilde_{k,j} *)
    atildeRaw = Table[
      aOther[[jj]] - apLocal * mOther[[jj]] / mpLocal,
      {jj, n}
    ] + Total[Table[polyExps[[k]] * rcMin[[k]], {k, Length[polyExps]}]];

    (* HasConstantTerm *)
    hasConst = And @@ Table[
      AnyTrue[newPolyList[[k]], (#[[2]] === ConstantArray[0, n]) &],
      {k, Length[newPolyList]}
    ];

    <|"pivot"      -> p,
      "mp"         -> mpLocal,
      "ap"         -> apLocal,
      "remainIdx"  -> rIdx,
      "mOther"     -> mOther,
      "atilde"     -> atildeRaw,
      "newPolys"   -> newPolyList,
      "hasConst"   -> hasConst,
      "rcMin"      -> rcMin
    |>
  ];

  (* ---- §3.5 pivot search with AMENDED ranking ---- *)
  (* Ranking (owner-approved):
     (1) HasConstantTerm -> True
     (2) |mp| == 1
     (3) tie-break: max of min_j Re[atilde_j]
     Note: this deviates from the plan's written order (plan says |mp|=1 first)
     but is explicitly approved per the B3 specification footnote. *)
  candidates = {};
  Do[
    If[mVec[[p]] != 0,
      Module[{res = tryPivot[p]},
        If[res =!= $Failed &&
           And @@ Table[
             With[{av = N[res["atilde"][[jj]]]},
               NumericQ[av] && Abs[Im[av]] < 10^-12 && Re[av] > 0
             ],
             {jj, n}
           ],
          AppendTo[candidates, res]
        ]
      ]
    ],
    {p, n1}
  ];

  If[candidates === {},
    (* Check if we had ANY complex-atilde candidates *)
    Module[{anyComplex = False},
      Do[
        If[mVec[[p]] != 0,
          Module[{res = tryPivot[p]},
            If[res =!= $Failed &&
               AnyTrue[N[res["atilde"]],
                 (NumericQ[#] && Abs[Im[#]] >= 10^-12) &],
              anyComplex = True
            ]
          ]
        ],
        {p, n1}
      ];
      If[anyComplex,
        Message[TropicalEval::liftcomplex, coneIndex];
        Return[$Failed]
      ]
    ];
    (* Build per-pivot atilde summary for the error message *)
    Module[{pivotSummary},
      pivotSummary = Table[
        If[mVec[[p]] != 0,
          Module[{res = tryPivot[p]},
            If[res =!= $Failed,
              {p, mVec[[p]], N[res["atilde"]]}
            ]
          ],
          Nothing
        ],
        {p, n1}
      ];
      Message[TropicalEval::liftnopivot, coneIndex, mVec, pivotSummary]
    ];
    Return[$Failed]
  ];

  (* AMENDED ranking: (1) HasConstantTerm; (2) |mp|=1; (3) max min Re[atilde] *)
  candidates = SortBy[candidates,
    {-Boole[#["hasConst"]],
     -Boole[Abs[#["mp"]] == 1],
     -Min[Re[N[#["atilde"]]]]} &
  ];
  bestPivot = candidates[[1]];

  pivotP      = bestPivot["pivot"];
  mp          = bestPivot["mp"];
  ap          = bestPivot["ap"];
  mOtherVec   = bestPivot["mOther"];   (* length n *)
  atildeVals  = bestPivot["atilde"];   (* length n *)
  reclearedPolys = bestPivot["newPolys"];

  If[verbose,
    Print["  PSL cone ", coneIndex, ": pivot p=", pivotP,
          " mp=", mp, " atilde=", N[atildeVals],
          " HasConstantTerm=", bestPivot["hasConst"]]
  ];

  (* ---- Step 4: domain constraint classification (§3.4) ---- *)
  logZ0       = Log[z0];
  allOtherZero = And @@ (# == 0 & /@ mOtherVec);

  If[allOtherZero,
    (* Constant-root case: y_p* = z0^(1/mp) *)
    If[N[z0^(1/mp)] > 1,
      Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
    ];
    domainClass = None;
    ,
    If[mp > 0,
      If[And @@ (#>= 0 & /@ mOtherVec) && N[z0] > 1,
        Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
      ];
      If[And @@ (#<= 0 & /@ mOtherVec) && N[z0] <= 1,
        domainClass = None;
        ,
        domainClass = <|
          "LogZ0"           -> logZ0,
          "MP"              -> mp,
          "IndicatorCoeffs" -> Table[mOtherVec[[jj]] / atildeVals[[jj]], {jj, n}]
        |>
      ],
      (* mp < 0 *)
      If[And @@ (#<= 0 & /@ mOtherVec) && N[z0] < 1,
        Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
      ];
      If[And @@ (#>= 0 & /@ mOtherVec) && N[z0] >= 1,
        domainClass = None;
        ,
        domainClass = <|
          "LogZ0"           -> logZ0,
          "MP"              -> mp,
          "IndicatorCoeffs" -> Table[mOtherVec[[jj]] / atildeVals[[jj]], {jj, n}]
        |>
      ]
    ]
  ];

  (* ---- Step 5: flatten via FlattenSector ---- *)
  prefactorBase = (Abs[detM] / Abs[mp]) * z0^(ap / mp - 1);
  fsResult = FlattenSector[reclearedPolys, atildeVals, prefactorBase];

  If[fsResult["IsDivergent"],
    Message[TropicalEval::liftdivergent, coneIndex, atildeVals];
    Return[$Failed]
  ];

  flattenedPolys = fsResult["FlattenedPolys"];
  prefactor      = fsResult["Prefactor"];

  (* ---- Step 6: assemble full SectorData ---- *)
  <|
    "ConeIndex"           -> coneIndex,
    "RayMatrix"           -> mMatrix,
    "DetM"                -> detM,
    "SelectedRays"        -> sdAug["SelectedRays"],
    "RawExponents"        -> sdAug["RawExponents"],
    "NewExponents"        -> atildeVals,
    "MinExponents"        -> bestPivot["rcMin"],
    "TransformedPolys"    -> clearedPolys,   (* pre-substitution cleared polys *)
    "ClearedPolys"        -> reclearedPolys,
    "FlattenedPolys"      -> flattenedPolys,
    "Prefactor"           -> prefactor,
    "IsDivergent"         -> False,
    "DivergentVariable"   -> 0,
    "Dimension"           -> n,
    "PolynomialExponents" -> polyExps,
    "MonomialExponents"   -> liftData["OriginalSpec"]["MonomialExponents"],
    (* Lifted-path extras *)
    "DomainConstraint"    -> domainClass,
    "LiftData"            -> liftData,
    "PivotIndex"          -> pivotP,
    "ZRow"                -> mVec,
    "AugmentedA"          -> a,
    "HasConstantTerm"     -> bestPivot["hasConst"],
    (* Log of the tropical monomial factor cleared out of each P_k, in the
       flattened sample coords, for the SplitRealImag oscillatory phase (see
       ProcessSector).  In the lifted sector the factor comes from BOTH the
       augmented clearing y^{d^aug_k} and the pivot substitution
       y_p = z0^{1/m_p} prod_{j!=p} y_j^{-m_j/m_p}:
         log|P_k| = (d^aug_{k,p}/m_p) log z0
                  + sum_j [ (d^aug_{k,j} - d^aug_{k,p} m_j/m_p) + rcMin_{k,j} ] / atilde_j * log_y'[j]
                  + log|Q_k| .                                                    *)
    "MonoFactorLog"       -> With[
      {dAug = sdAug["MinExponents"], rIdx = bestPivot["remainIdx"],
       rcM = bestPivot["rcMin"], logZ0v = Log[z0]},
      Table[
        <|"Const" -> (dAug[[k, pivotP]]/mp) logZ0v,
          "Coeffs" -> Table[
            ((dAug[[k, rIdx[[jj]]]] - dAug[[k, pivotP]] mVec[[rIdx[[jj]]]]/mp)
              + rcM[[k, jj]]) / atildeVals[[jj]],
            {jj, n}]|>,
        {k, Length[dAug]}]
    ]
  |>
];


(* --------------------------------------------------------------------------
   ValidateDecomposition
   Cross-check sector sum against direct NIntegrate.
   -------------------------------------------------------------------------- *)

ValidateDecomposition[integrandSpec_Association, fanData_List,
                      testKinematics_List, precisionGoal_Integer: 3] :=
Module[
  {polys, monoExps, polyExps, vars, n, dualVertices, simplexList,
   directIntegrand, directResult, sectorResults, sectorSum,
   relError, kinRules, allSectorData},

  polys    = integrandSpec["Polynomials"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  vars     = integrandSpec["Variables"];
  n        = Length[vars];
  kinRules = testKinematics;

  {dualVertices, simplexList} = fanData;

  (* Direct NIntegrate of the original integrand *)
  directIntegrand = (Times @@ MapThread[Power, {vars, monoExps}]) *
    (Times @@ MapThread[Power, {polys, polyExps}]) /. kinRules;

  directResult = NIntegrate[
    directIntegrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 20,
    PrecisionGoal -> precisionGoal + 1,
    Method -> "GlobalAdaptive"
  ];

  (* Process all sectors *)
  allSectorData = Table[
    ProcessSector[integrandSpec /. kinRules, dualVertices,
                  simplexList[[s]], s],
    {s, Length[simplexList]}
  ];

  (* For each sector, evaluate via NIntegrate on [0,1]^n *)
  sectorResults = Table[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, integrand,
            polyVals, result},
      sd = allSectorData[[s]];
      If[sd === $Failed, 0,
        If[sd["IsDivergent"],
          Message[TropicalEval::divergentinput, {s}, {sd["NewExponents"]}];
          Return[$Failed],
          flatPolys = sd["FlattenedPolys"];
          pExps     = sd["PolynomialExponents"] /. kinRules;
          pf        = sd["Prefactor"] /. kinRules;
          dim       = sd["Dimension"];
          yVars     = Table[Unique["yv"], {dim}];

          (* Evaluate using log-exp form to avoid numerical issues *)
          polyVals = Table[
            Total[
              Table[
                Module[{coeff, alphas, logY},
                  coeff  = mono[[1]] /. kinRules;
                  alphas = mono[[2]] /. kinRules;
                  logY   = Log /@ yVars;
                  coeff * Exp[Total[alphas * logY]]
                ],
                {mono, flatPolys[[j]]}
              ]
            ],
            {j, Length[flatPolys]}
          ];

          integrand = pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals, pExps}
            ];

          result = Quiet@NIntegrate[
            integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 15,
            PrecisionGoal -> precisionGoal,
            Method -> "GlobalAdaptive"
          ];
          result
        ]
      ]
    ],
    {s, Length[simplexList]}
  ];

  sectorSum = Total[sectorResults];
  relError  = Abs[(sectorSum - directResult) / directResult];

  If[relError > 10^(-precisionGoal + 1),
    Message[TropicalEval::validate, "Decomposition", relError,
            10^(-precisionGoal + 1)]
  ];

  <|"DirectResult" -> directResult, "SectorSum" -> sectorSum,
    "RelativeError" -> relError, "SectorResults" -> sectorResults|>
];

(* --------------------------------------------------------------------------
   ValidateLiftedDecomposition  (B4)
   Sibling of ValidateDecomposition for lifted integrands.
   Direct integral from originalSpec; sectors via ProcessSectorLifted.
   EmptyDomain sectors contribute 0 and are listed in DroppedSectors.
   -------------------------------------------------------------------------- *)

ValidateLiftedDecomposition[originalSpec_Association,
                            liftedSpec_Association,
                            liftedFanData_List,
                            liftData_Association,
                            testKinematics_List,
                            precisionGoal_Integer: 3] :=
Module[
  {polys, monoExps, polyExps, vars, n,
   dualVertices, simplexList, kinRules, pg,
   directIntegrand, directResult,
   sectorResults, droppedSectors, sectorSum, relError},

  polys    = originalSpec["Polynomials"];
  monoExps = originalSpec["MonomialExponents"];
  polyExps = originalSpec["PolynomialExponents"];
  vars     = originalSpec["Variables"];
  n        = Length[vars];
  kinRules = testKinematics;
  pg       = precisionGoal;

  {dualVertices, simplexList} = liftedFanData;

  (* Direct NIntegrate of the ORIGINAL integrand on [0,Inf)^n *)
  directIntegrand = (Times @@ MapThread[Power, {vars, monoExps}]) *
    (Times @@ MapThread[Power, {polys, polyExps}]) /. kinRules;

  directResult = NIntegrate[
    directIntegrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 20,
    PrecisionGoal -> pg + 1,
    Method -> "GlobalAdaptive"
  ];

  droppedSectors = {};

  (* Per-sector NIntegrate via ProcessSectorLifted *)
  sectorResults = Table[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, integrand,
            polyVals, dc, result, logYpStar, logZ0num, mpNum, icNum},

      sd = ProcessSectorLifted[liftedSpec, dualVertices,
                               simplexList[[s]], s, liftData];

      Which[
        sd === $Failed,
          Return[$Failed],

        AssociationQ[sd] && KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"],
          AppendTo[droppedSectors, sd["ConeIndex"]];
          0,

        True,
          flatPolys = sd["FlattenedPolys"];
          pExps     = sd["PolynomialExponents"] /. kinRules;
          pf        = sd["Prefactor"] /. kinRules;
          dim       = sd["Dimension"];
          dc        = sd["DomainConstraint"];
          yVars     = Table[Unique["yv"], {dim}];

          (* log-exp evaluation copied from ValidateDecomposition *)
          polyVals = Table[
            Total[
              Table[
                Module[{coeff, alphas, logY},
                  coeff  = mono[[1]] /. kinRules;
                  alphas = mono[[2]] /. kinRules;
                  logY   = Log /@ yVars;
                  coeff * Exp[Total[alphas * logY]]
                ],
                {mono, flatPolys[[j]]}
              ]
            ],
            {j, Length[flatPolys]}
          ];

          integrand = pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals, pExps}
            ];

          (* Multiply by domain indicator when DomainConstraint is present *)
          If[dc =!= None,
            logZ0num = N[dc["LogZ0"]];
            mpNum    = N[dc["MP"]];
            icNum    = N[dc["IndicatorCoeffs"]];
            logYpStar = (logZ0num - Total[icNum * (Log /@ yVars)]) / mpNum;
            integrand = integrand * Boole[logYpStar <= 0]
          ];

          result = Quiet@NIntegrate[
            integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 15,
            PrecisionGoal -> pg,
            Method -> "GlobalAdaptive"
          ];
          result
      ]
    ],
    {s, Length[simplexList]}
  ];

  (* Propagate $Failed from any sector *)
  If[MemberQ[sectorResults, $Failed], Return[$Failed]];

  sectorSum = Total[sectorResults];
  relError  = Abs[(sectorSum - directResult) / directResult];

  If[relError > 10^(-pg + 1),
    Message[TropicalEval::validate, "LiftedDecomposition", relError,
            10^(-pg + 1)]
  ];

  <|"DirectResult"  -> directResult,
    "SectorSum"     -> sectorSum,
    "RelativeError" -> relError,
    "SectorResults" -> sectorResults,
    "DroppedSectors"-> droppedSectors|>
];

(* --------------------------------------------------------------------------
   2D benchmark unit test (hard-coded)
   P = 1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2
   rho_1 = (0,1), rho_2 = (1,1), A_1 = A_2 = 0, polynomial exponent -A
   Expected effective: a_1 = 2A - 1, a_2 = 3A - 2
   -------------------------------------------------------------------------- *)

RunBenchmark2D[] := Module[
  {dualVerts, simplex, integrandSpec, sd, aExpected, aGot, A, pass},

  dualVerts = {{0, 1}, {1, 1}};
  simplex   = {1, 2};

  A = Symbol["Abench"];

  integrandSpec = <|
    "Polynomials"       -> {1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> {x[1], x[2]},
    "KinematicSymbols"  -> {},
    "RegulatorSymbol"   -> None
  |>;

  sd = ProcessSector[integrandSpec, dualVerts, simplex, 1];

  (* Expected effective exponents:
     Raw a = (monoExps+1).M = {1,1}.{{0,-1},{-1,-1}} = {-1,-2}
     Transformed polynomial monomials in y have exponents:
       1         -> {0,0}.M = {0,0}
       2 x1^2    -> {2,0}.M = {0,-2}
       x2^2      -> {0,2}.M = {-2,-2}
       x1 x2^2   -> {1,2}.M = {-2,-3}     <- dominant (minimum)
       3 x1^2 x2 -> {2,1}.M = {-1,-3}
     minExp = {-2, -3}
     effA = rawA + B*minExp = {-1,-2} + (-A)*{-2,-3} = {2A-1, 3A-2} *)

  aExpected = {2 A - 1, 3 A - 2};
  aGot      = sd["NewExponents"];

  pass = TrueQ[Simplify[aGot - aExpected] === {0, 0}] ||
         (aGot === aExpected);
  If[!pass,
    Print["BENCHMARK 2D FAILED: expected a_eff = ", aExpected,
          " got ", aGot];
    ,
    Print["Benchmark 2D: PASSED (a_eff = ", aGot, ")"];
  ];
  pass
];



(* ============================================================================
   MODULE 3: C++ CODE GENERATION
   ============================================================================ *)

(* --------------------------------------------------------------------------
   MmaToC: Convert Mathematica expression to C++ string
   Recursive pattern-matching converter (not CForm-based).
   -------------------------------------------------------------------------- *)

MmaToC[expr_, paramMap_Association: <||>] := mmaToCInternal[expr, paramMap];

(* Integers -> double literals *)
mmaToCInternal[n_Integer, _] := ToString[n] <> ".0";

(* Rationals -> (a.0/b.0) *)
mmaToCInternal[r_Rational, _] :=
  "(" <> ToString[Numerator[r]] <> ".0/" <>
  ToString[Denominator[r]] <> ".0)";

(* Reals -> string *)
mmaToCInternal[r_Real, _] := ToString[CForm[r]];

(* Complex numbers -> cx(re, im) *)
mmaToCInternal[Complex[re_, im_], pm_] :=
  "cx(" <> mmaToCInternal[re, pm] <> ", " <>
  mmaToCInternal[im, pm] <> ")";

(* Symbols -> parameter lookup or literal *)
mmaToCInternal[s_Symbol, pm_] :=
  If[KeyExistsQ[pm, s],
    pm[s],
    ToString[s]
  ];

(* Subscripted variables like x[i] *)
mmaToCInternal[s_Symbol[i_Integer], pm_] :=
  If[KeyExistsQ[pm, s[i]],
    pm[s[i]],
    ToString[s] <> "[" <> ToString[i] <> "]"
  ];

(* Power *)
mmaToCInternal[Power[base_, -1], pm_] :=
  "(1.0/" <> mmaToCInternal[base, pm] <> ")";

mmaToCInternal[Power[base_, 1/2], pm_] :=
  "std::sqrt(" <> mmaToCInternal[base, pm] <> ")";

mmaToCInternal[Power[base_, -1/2], pm_] :=
  "(1.0/std::sqrt(" <> mmaToCInternal[base, pm] <> "))";

mmaToCInternal[Power[E, exp_], pm_] :=
  "std::exp(" <> mmaToCInternal[exp, pm] <> ")";

mmaToCInternal[Power[base_, n_Integer], pm_] :=
  "std::pow(" <> mmaToCInternal[base, pm] <> ", " <>
  ToString[n] <> ".0)";

mmaToCInternal[Power[base_, exp_], pm_] :=
  "std::pow(" <> mmaToCInternal[base, pm] <> ", " <>
  mmaToCInternal[exp, pm] <> ")";

(* Log *)
mmaToCInternal[Log[arg_], pm_] :=
  "std::log(" <> mmaToCInternal[arg, pm] <> ")";

(* Exp *)
mmaToCInternal[Exp[arg_], pm_] :=
  "std::exp(" <> mmaToCInternal[arg, pm] <> ")";

(* Abs *)
mmaToCInternal[Abs[arg_], pm_] :=
  "std::abs(" <> mmaToCInternal[arg, pm] <> ")";

(* Re, Im *)
mmaToCInternal[Re[arg_], pm_] :=
  "(" <> mmaToCInternal[arg, pm] <> ").real()";

mmaToCInternal[Im[arg_], pm_] :=
  "(" <> mmaToCInternal[arg, pm] <> ").imag()";

(* Plus -> infix with parens.
   HoldPattern is needed because Mathematica evaluates Plus[args__] and
   Times[args__] before storing the DownValue inside BeginPackage. *)
mmaToCInternal[HoldPattern[Plus[args__]], pm_] :=
  "(" <> StringRiffle[mmaToCInternal[#, pm] & /@ {args}, " + "] <> ")";

(* Times: handle negation and general products *)
mmaToCInternal[HoldPattern[Times[-1, rest__]], pm_] :=
  "(-" <> mmaToCInternal[Times[rest], pm] <> ")";

mmaToCInternal[HoldPattern[Times[args__]], pm_] :=
  "(" <> StringRiffle[mmaToCInternal[#, pm] & /@ {args}, " * "] <> ")";

(* Fallback: use CForm and warn *)
mmaToCInternal[expr_, pm_] := Module[{str},
  str = ToString[CForm[expr]];
  str = StringReplace[str, {
    "Power(" ~~ a__ ~~ "," ~~ b__ ~~ ")" :>
      "std::pow(" <> a <> ", " <> b <> ")",
    "Sqrt(" ~~ a__ ~~ ")" :> "std::sqrt(" <> a <> ")"
  }];
  str
];

(* --------------------------------------------------------------------------
   Helper: EmitCoeff
   Convert a (possibly complex) coefficient to a C++ literal string.

   A genuinely complex NUMERIC coefficient (Im != 0) is emitted directly as
   cx(re, im) from its numeric value.  This is needed because lifting a
   complex coefficient C carries its phase into a residual c = C/z0^k (e.g.
   (1+I)/Sqrt[2]); Simplify may leave such residuals in a symbolic form whose
   head the recursive MmaToC converter does not handle, which would otherwise
   reach the CForm fallback and print an invalid C++ token (Complex(...), Pi,
   E, ...).  Emitting the numeric value as cx(re, im) is robust regardless of
   the symbolic form.

   Real numeric coefficients and symbolic (kinematic) coefficients are routed
   through mmaToCInternal unchanged, so the real-coefficient code path is
   byte-for-byte identical to before.
   -------------------------------------------------------------------------- *)

EmitCoeff[coeff_, paramMap_Association] :=
  If[NumericQ[coeff] && N[Im[coeff]] != 0,
    Module[{cn = N[coeff]},
      "cx(" <> ToString[CForm[Re[cn]]] <> ", " <>
              ToString[CForm[Im[cn]]] <> ")"
    ],
    mmaToCInternal[coeff, paramMap]
  ];

(* --------------------------------------------------------------------------
   Helper: Generate C++ for one monomial sum (polynomial evaluation)
   -------------------------------------------------------------------------- *)

GenerateMonomialSumCpp[flatPolys_List, polyIndex_Integer,
                       paramMap_Association, dim_Integer,
                       varPrefix_String: "log_y"] :=
Module[{lines, polyVar},
  polyVar = "P" <> ToString[polyIndex];
  lines = {"    cx " <> polyVar <> "(0.0, 0.0);"};

  Do[
    Module[{coeff, alphas, coeffStr, expTerms, expStr},
      coeff    = mono[[1]];
      alphas   = mono[[2]];
      coeffStr = EmitCoeff[coeff, paramMap];

      expTerms = Table[
        If[TrueQ[alphas[[i]] == 0],
          Nothing,
          mmaToCInternal[alphas[[i]], paramMap] <> " * " <>
            varPrefix <> "[" <> ToString[i - 1] <> "]"
        ],
        {i, dim}
      ];

      expStr = If[Length[expTerms] == 0,
        "0.0",
        StringRiffle[expTerms, " + "]
      ];

      AppendTo[lines,
        "    " <> polyVar <> " += " <> coeffStr <>
        " * std::exp(" <> expStr <> ");"
      ];
    ],
    {mono, flatPolys}
  ];

  StringRiffle[lines, "\n"]
];

(* --------------------------------------------------------------------------
   GenerateCppMonteCarlo
   Main code generation function.
   -------------------------------------------------------------------------- *)

Options[GenerateCppMonteCarlo] = {
  "NSamples"       -> 1000000,
  "MaxDim"         -> 20,
  "SeedBase"       -> 42,
  "Integrator"     -> "MC",        (* "MC" | "VEGAS" *)
  (* VEGAS tuning -- used only when Integrator == "VEGAS" *)
  "VegasEpsRel"    -> Automatic,    (* Automatic => 1e-4 here; the driver
                                       resolves it from PrecisionGoal first *)
  "VegasSeed"      -> 0,            (* 0 => Sobol low-discrepancy (REPORT.md) *)
  "VegasNStart"    -> 1000,
  "VegasNIncrease" -> 500,
  "VegasNBatch"    -> 1000,
  "VegasMinEval"   -> 0,
  "ImagExponents"  -> None          (* list of Im[B_k] for oscillatory phase, or None *)
};

GenerateCppMonteCarlo[convergentSectors_List, divergentSectors_List,
                      integrandSpec_Association, outputFile_String,
                      OptionsPattern[]] :=
Module[
  {kinSyms, paramMap, nParams,
   code, integrandFuncs, integrandDims, nIntegrands,
   nConvergent,
   maxDim, seedBase, nSamples,
   integrator, isVegas, vegasEpsRel, vegasEpsRelStr,
   vegasSeed, vegasNStart, vegasNIncrease, vegasNBatch, vegasMinEval,
   imagExpsOpt, hasImagPhase},

  If[divergentSectors =!= {},
    Message[TropicalEval::nodivergent];
    Return[$Failed]
  ];

  (* --- Integrator selection (default "MC" leaves all output byte-identical) --- *)
  integrator = OptionValue["Integrator"];
  If[!MatchQ[integrator, "MC" | "VEGAS"],
    Message[TropicalEval::badintegrator, integrator];
    Return[$Failed]
  ];
  isVegas = (integrator === "VEGAS");

  (* Resolve VEGAS tuning to C++ literals.  VegasEpsRel == Automatic here
     (e.g. a direct call, not via the driver) falls back to 1e-4. *)
  vegasEpsRel    = OptionValue["VegasEpsRel"];
  If[vegasEpsRel === Automatic, vegasEpsRel = 1.*^-4];
  vegasEpsRelStr = ToString[CForm[N[vegasEpsRel, 17]]];
  vegasSeed      = Round[OptionValue["VegasSeed"]];
  vegasNStart    = Round[OptionValue["VegasNStart"]];
  vegasNIncrease = Round[OptionValue["VegasNIncrease"]];
  vegasNBatch    = Round[OptionValue["VegasNBatch"]];
  vegasMinEval   = Round[OptionValue["VegasMinEval"]];

  kinSyms  = integrandSpec["KinematicSymbols"];
  nParams  = Length[kinSyms];
  maxDim   = OptionValue["MaxDim"];
  seedBase = OptionValue["SeedBase"];
  nSamples = OptionValue["NSamples"];

  imagExpsOpt  = OptionValue["ImagExponents"];
  hasImagPhase = imagExpsOpt =!= None &&
    AnyTrue[imagExpsOpt, (Abs[N[#]] > 1*^-15 &)];

  (* Build parameter map: kinematic symbol -> params[i] *)
  paramMap = Association @@ Table[
    kinSyms[[i]] -> ("params[" <> ToString[i - 1] <> "]"),
    {i, nParams}
  ];

  integrandFuncs = {};
  integrandDims  = {};
  nConvergent = 0;

  (* --- Generate convergent sector integrands --- *)
  Do[
    Module[{sd, flatPolys, polyExps, prefactor, dim, funcName,
            funcCode, polyCode, prodCode},
      sd        = convergentSectors[[s]];
      flatPolys = sd["FlattenedPolys"];
      polyExps  = sd["PolynomialExponents"];
      prefactor = sd["Prefactor"];
      dim       = sd["Dimension"];
      funcName  = "integrand_conv_" <> ToString[s - 1];

      funcCode = "inline cx " <> funcName <>
        "(const double* y, const double* params) {\n";
      funcCode = funcCode <>
        "    // Convergent sector " <> ToString[sd["ConeIndex"]] <> "\n";
      funcCode = funcCode <>
        "    double log_y[" <> ToString[dim] <> "];\n";
      funcCode = funcCode <>
        "    for (int i = 0; i < " <> ToString[dim] <>
        "; i++)\n";
      funcCode = funcCode <>
        "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

      (* B5: domain indicator for lifted sectors *)
      Module[{dc},
        dc = Lookup[sd, "DomainConstraint", None];
        If[dc =!= None,
          Module[{logZ0str, mpStr, icList, icTerms, sumStr},
            logZ0str = mmaToCInternal[N[dc["LogZ0"]], paramMap];
            mpStr    = mmaToCInternal[N[dc["MP"]], paramMap];
            icList   = N[dc["IndicatorCoeffs"]];
            icTerms  = Table[
              mmaToCInternal[icList[[i]], paramMap] <>
              " * log_y[" <> ToString[i - 1] <> "]",
              {i, Length[icList]}
            ];
            sumStr = If[Length[icTerms] == 0,
              "0.0",
              StringRiffle[icTerms, " + "]
            ];
            funcCode = funcCode <>
              "    // lifted-sector domain indicator\n";
            funcCode = funcCode <>
              "    double log_ypstar = (" <> logZ0str <>
              " - (" <> sumStr <> ")) * (1.0/" <> mpStr <> ");\n";
            funcCode = funcCode <>
              "    if (log_ypstar > 0.0) return cx(0.0, 0.0);\n\n";
          ]
        ]
      ];

      Do[
        funcCode = funcCode <>
          GenerateMonomialSumCpp[flatPolys[[j]], j - 1, paramMap, dim] <>
          "\n\n";,
        {j, Length[flatPolys]}
      ];

      funcCode = funcCode <> "    cx result = " <>
        EmitCoeff[prefactor, paramMap] <> ";\n";

      Do[
        funcCode = funcCode <>
          "    result *= std::exp(" <>
          mmaToCInternal[polyExps[[j]], paramMap] <>
          " * std::log(P" <> ToString[j - 1] <> "));\n";,
        {j, Length[polyExps]}
      ];

      (* Oscillatory phase from complex polynomial exponents (SplitRealImag mode):
         exp(i * sum_k imB[k] * log|P_k|).  Only emitted when any imB is non-zero.
         log|P_k| is NOT just log|Q_k| (the cleared-polynomial value evaluated
         above): the tropical factoring removed a monomial factor y^{d_k} from
         P_k, and for a lifted sector the pivot substitution added another.  That
         factor's log is the per-sector "MonoFactorLog" linear form (Const_k +
         sum_i Coeffs_{k,i} log_y[i]); omitting it (the pre-fix behavior) gives a
         wrong oscillatory phase whenever any d_k != 0. *)
      If[hasImagPhase,
        Module[{imBStr, mfl, phaseTerms, phaseSum},
          imBStr = StringRiffle[
            ToString[CForm[N[#]]] & /@ imagExpsOpt,
            ", "
          ];
          funcCode = funcCode <>
            "    const double imB[] = {" <> imBStr <> "};\n";
          mfl = Lookup[sd, "MonoFactorLog", None];
          phaseTerms = Table[
            Module[{constK, coeffsK, monoTerms, logQ},
              logQ = "std::log(std::abs(P" <> ToString[j - 1] <> "))";
              If[mfl === None,
                (* fallback: cleared-polynomial value only (legacy behavior) *)
                "imB[" <> ToString[j - 1] <> "] * (" <> logQ <> ")",
                constK  = N[mfl[[j]]["Const"]];
                coeffsK = N[mfl[[j]]["Coeffs"]];
                monoTerms = Table[
                  If[TrueQ[coeffsK[[i]] == 0.], Nothing,
                    "(" <> ToString[CForm[coeffsK[[i]]]] <> ") * log_y[" <>
                    ToString[i - 1] <> "]"],
                  {i, Length[coeffsK]}];
                "imB[" <> ToString[j - 1] <> "] * (" <>
                  StringRiffle[
                    Join[
                      If[TrueQ[constK == 0.], {},
                        {"(" <> ToString[CForm[constK]] <> ")"}],
                      monoTerms, {logQ}],
                    " + "] <> ")"
              ]
            ],
            {j, Length[imagExpsOpt]}
          ];
          phaseSum = StringRiffle[phaseTerms, " + "];
          funcCode = funcCode <>
            "    result *= std::exp(cx(0.0, " <> phaseSum <> "));\n"
        ]
      ];

      funcCode = funcCode <> "    return result;\n}\n";

      AppendTo[integrandFuncs, funcCode];
      AppendTo[integrandDims, dim];
      nConvergent++;
    ],
    {s, Length[convergentSectors]}
  ];

  nIntegrands = Length[integrandFuncs];

  (* --- Assemble the full C++ file --- *)
  code = "// Auto-generated by TropicalEval`GenerateCppMonteCarlo\n";
  code = code <> "// " <> ToString[nConvergent] <> " convergent integrands\n\n";

  code = code <> "#include <complex>\n";
  code = code <> "#include <cmath>\n";
  code = code <> "#include <random>\n";
  code = code <> "#include <fstream>\n";
  code = code <> "#include <vector>\n";
  code = code <> "#include <iostream>\n";
  code = code <> "#include <string>\n";
  code = code <> "#include <cassert>\n";
  code = code <> "#include <cstdlib>\n";
  code = code <> "#include <array>\n";
  code = code <> "#ifdef _OPENMP\n";
  code = code <> "#include <omp.h>\n";
  code = code <> "#endif\n\n";

  code = code <> "using cx = std::complex<double>;\n\n";

  Do[
    code = code <> integrandFuncs[[i]] <> "\n";,
    {i, nIntegrands}
  ];

  code = code <>
    "// Function pointer type\n" <>
    "using IntegrandFunc = cx(*)(const double*, const double*);\n\n";

  code = code <> "IntegrandFunc integrand_table[] = {\n";
  Module[{allNames},
    allNames = {};
    Do[AppendTo[allNames, "integrand_conv_" <> ToString[i - 1]],
       {i, nConvergent}];
    code = code <> "    " <>
      StringRiffle[allNames, ",\n    "] <> "\n";
  ];
  code = code <> "};\n\n";

  code = code <> "int integrand_dim[] = {" <>
    StringRiffle[ToString /@ integrandDims, ", "] <> "};\n";
  code = code <> "const int N_INTEGRANDS = " <>
    ToString[nIntegrands] <> ";\n";
  code = code <> "const int N_PARAMS = " <>
    ToString[nParams] <> ";\n";
  code = code <> "const int MAX_DIM = " <>
    ToString[maxDim] <> ";\n\n";

  (* --- VEGAS support block (emitted only for Integrator -> "VEGAS").
     Placed here so IntegrandFunc, integrand_table/_dim, N_*, and MAX_DIM are
     all in scope for vegas_wrap.  Nothing is emitted in the MC case, so the
     MC source stays byte-identical. --- *)
  If[isVegas,
    code = code <> "// ---- CUBA Vegas support (Integrator -> \"VEGAS\") ----\n";
    code = code <> "extern \"C\" {\n#include <cuba.h>\n}\n";
    code = code <> "static const double VEGAS_EPSREL    = " <> vegasEpsRelStr <> ";\n";
    code = code <> "static const double VEGAS_EPSABS    = 1e-300;\n";
    code = code <> "static const int    VEGAS_SEED      = " <> ToString[vegasSeed] <> ";\n";
    code = code <> "static const int    VEGAS_NSTART    = " <> ToString[vegasNStart] <> ";\n";
    code = code <> "static const int    VEGAS_NINCREASE = " <> ToString[vegasNIncrease] <> ";\n";
    code = code <> "static const int    VEGAS_NBATCH    = " <> ToString[vegasNBatch] <> ";\n";
    code = code <> "static const int    VEGAS_MINEVAL   = " <> ToString[vegasMinEval] <> ";\n\n";
    code = code <> "struct SectorCtx { IntegrandFunc fn; const double* params; };\n";
    code = code <> "static int vegas_wrap(const int* ndim, const cubareal x[], const int* /*ncomp*/,\n";
    code = code <> "                      cubareal ff[], void* ud) {\n";
    code = code <> "    const SectorCtx* ctx = static_cast<const SectorCtx*>(ud);\n";
    code = code <> "    double y[MAX_DIM];\n";
    code = code <> "    for (int i = 0; i < *ndim; ++i) y[i] = (double)x[i];\n";
    code = code <> "    cx v = ctx->fn(y, ctx->params);\n";
    code = code <> "    ff[0] = v.real();\n";
    code = code <> "    ff[1] = v.imag();\n";
    code = code <> "    return 0;\n";
    code = code <> "}\n\n";
  ];

  (* Main function *)
  code = code <> "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [n_samples] [n_threads] [seed_base]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    int n_samples = (argc > 3) ? std::atoi(argv[3]) : " <>
    ToString[nSamples] <> ";\n";
  code = code <> "    int n_threads = (argc > 4) ? std::atoi(argv[4]) : 1;\n";
  code = code <> "    // Optional runtime seed override (argv[5]); falls back to compile-time SeedBase.\n";
  code = code <> "    uint64_t seed_base = (argc > 5) ? std::strtoull(argv[5], nullptr, 10) : " <>
    ToString[seedBase] <> "ULL;\n";
  code = code <> "#ifdef _OPENMP\n";
  code = code <> "    if (n_threads == 1) n_threads = omp_get_max_threads();\n";
  code = code <> "    omp_set_num_threads(n_threads);\n";
  code = code <> "#endif\n\n";

  code = code <> "    // Read kinematic data\n";
  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) {\n";
  code = code <> "        std::cerr << \"Cannot open \" << input_file << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        // No kinematic parameters: read count from file, default 1\n";
  code = code <> "        int count = 1;\n";
  code = code <> "        fin >> count;\n";
  code = code <> "        if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++)\n";
  code = code <> "            kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val;\n";
  code = code <> "        std::vector<double> row;\n";
  code = code <> "        while (fin >> val) {\n";
  code = code <> "            row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) {\n";
  code = code <> "                kinematic_data.push_back(row);\n";
  code = code <> "                row.clear();\n";
  code = code <> "            }\n";
  code = code <> "        }\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points\" << std::endl;\n\n";

  code = code <> "    std::vector<std::array<double, 4>> results(n_kp);\n\n";

  If[isVegas,
    (* ====================================================================
       VEGAS compute region.  One Vegas call per (kp, sector); ncomp = 2
       (component 0 = Re, 1 = Im).  Errors summed in quadrature over sectors.
       Produces exactly the same results[kp] 4-tuple as MC, so the write loop
       below is shared.  seed = 0 (Sobol) => each kp is deterministic and
       independent of thread count (see VEGAS_PLAN.md §2.3-2.4).
       ==================================================================== *)
    code = code <> "    (void)seed_base;  // VEGAS uses VEGAS_SEED (Sobol) below, not the MC seed\n";
    code = code <> "    const int zero = 0;\n";
    code = code <> "    cubacores(&zero, &zero);   // single process: deterministic, macOS-safe\n\n";

    code = code <> "    #pragma omp parallel for schedule(dynamic)\n";
    code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
    code = code <> "        static double dummy_params[1] = {0.0};\n";
    code = code <> "        const double* params = (N_PARAMS > 0) ? kinematic_data[kp].data() : dummy_params;\n\n";

    code = code <> "        double total_re = 0.0, total_im = 0.0;\n";
    code = code <> "        double total_var_re = 0.0, total_var_im = 0.0;\n\n";

    code = code <> "        for (int s = 0; s < N_INTEGRANDS; s++) {\n";
    code = code <> "            int dim = integrand_dim[s];\n";
    code = code <> "            SectorCtx ctx{ integrand_table[s], params };\n";
    code = code <> "            int neval = 0, fail = 0;\n";
    code = code <> "            cubareal integ[2], err[2], prob[2];\n";
    code = code <> "            Vegas(dim, 2, vegas_wrap, &ctx, 1,\n";
    code = code <> "                  VEGAS_EPSREL, VEGAS_EPSABS, 0 /*flags*/, VEGAS_SEED,\n";
    code = code <> "                  VEGAS_MINEVAL, n_samples /*maxeval per sector*/,\n";
    code = code <> "                  VEGAS_NSTART, VEGAS_NINCREASE, VEGAS_NBATCH,\n";
    code = code <> "                  0 /*gridno*/, nullptr, nullptr,\n";
    code = code <> "                  &neval, &fail, integ, err, prob);\n";
    code = code <> "            total_re     += integ[0]; total_im     += integ[1];\n";
    code = code <> "            total_var_re += err[0]*err[0]; total_var_im += err[1]*err[1];\n";
    code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
    code = code <> "            if (kp == 0)\n";
    code = code <> "                std::cerr << \"Sector \" << s << \" (Vegas): est=(\" << integ[0] << \",\"\n";
    code = code <> "                          << integ[1] << \") err=(\" << err[0] << \",\" << err[1]\n";
    code = code <> "                          << \") fail=\" << fail << \" prob=\" << prob[0]\n";
    code = code <> "                          << \" neval=\" << neval << std::endl;\n";
    code = code <> "#endif\n";
    code = code <> "        }\n\n";

    code = code <> "        results[kp] = {total_re, total_im,\n";
    code = code <> "                       std::sqrt(total_var_re), std::sqrt(total_var_im)};\n";
    code = code <> "    }\n\n";
    ,
    (* ====================================================================
       MC compute region -- byte-identical to the pre-VEGAS pipeline.
       ==================================================================== *)
    code = code <> "    #pragma omp parallel for schedule(dynamic)\n";
    code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
    code = code <> "        const double* params = kinematic_data[kp].data();\n";
    code = code <> "        uint64_t seed = seed_base + (uint64_t)kp;\n";
    code = code <> "        std::mt19937_64 rng(seed);\n";
    code = code <> "        std::uniform_real_distribution<double> dist(0.0, 1.0);\n\n";

    code = code <> "        double total_re = 0.0, total_im = 0.0;\n";
    code = code <> "        double total_var_re = 0.0, total_var_im = 0.0;\n\n";

    code = code <> "        for (int s = 0; s < N_INTEGRANDS; s++) {\n";
    code = code <> "            int dim = integrand_dim[s];\n";
    code = code <> "            double mean_re = 0.0, mean_im = 0.0;\n";
    code = code <> "            double M2_re = 0.0, M2_im = 0.0;\n";

    code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
    code = code <> "            int nan_count = 0;\n";
    code = code <> "            double max_mag = 0.0;\n";
    code = code <> "#endif\n\n";

    code = code <> "            for (int k = 0; k < n_samples; k++) {\n";
    code = code <> "                double y[MAX_DIM];\n";
    code = code <> "                for (int i = 0; i < dim; i++) y[i] = dist(rng);\n\n";

    code = code <> "                cx val = integrand_table[s](y, params);\n\n";

    code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
    code = code <> "                if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {\n";
    code = code <> "                    nan_count++;\n";
    code = code <> "                    if (nan_count <= 5) {\n";
    code = code <> "                        std::cerr << \"NaN/Inf in sector \" << s << \" kp=\" << kp << \" y=[\";";
    code = code <> "\n                        for (int i = 0; i < dim; i++) std::cerr << y[i] << \" \";\n";
    code = code <> "                        std::cerr << \"]\" << std::endl;\n";
    code = code <> "                    }\n";
    code = code <> "                    continue;\n";
    code = code <> "                }\n";
    code = code <> "                double mag = std::abs(val);\n";
    code = code <> "                if (mag > max_mag) max_mag = mag;\n";
    code = code <> "#endif\n\n";

    code = code <> "                double d_re = val.real() - mean_re;\n";
    code = code <> "                mean_re += d_re / (k + 1);\n";
    code = code <> "                M2_re += d_re * (val.real() - mean_re);\n";
    code = code <> "                double d_im = val.imag() - mean_im;\n";
    code = code <> "                mean_im += d_im / (k + 1);\n";
    code = code <> "                M2_im += d_im * (val.imag() - mean_im);\n";
    code = code <> "            }\n\n";

    code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
    code = code <> "            if (kp == 0) {\n";
    code = code <> "                std::cerr << \"Sector \" << s << \": mean=(\" << mean_re << \",\" << mean_im\n";
    code = code <> "                          << \") max_mag=\" << max_mag;\n";
    code = code <> "                if (nan_count > 0)\n";
    code = code <> "                    std::cerr << \" NaN_count=\" << nan_count;\n";
    code = code <> "                double mean_mag = std::sqrt(mean_re*mean_re + mean_im*mean_im);\n";
    code = code <> "                if (mean_mag > 0 && max_mag / mean_mag > 1000)\n";
    code = code <> "                    std::cerr << \" WARNING: large fluctuations\";\n";
    code = code <> "                std::cerr << std::endl;\n";
    code = code <> "            }\n";
    code = code <> "            if ((double)nan_count / n_samples > 0.001)\n";
    code = code <> "                std::cerr << \"WARNING: >0.1%% NaN in sector \" << s << \" kp=\" << kp << std::endl;\n";
    code = code <> "#endif\n\n";

    code = code <> "            total_re += mean_re;\n";
    code = code <> "            total_im += mean_im;\n";
    code = code <> "            total_var_re += M2_re / ((double)n_samples * (n_samples - 1));\n";
    code = code <> "            total_var_im += M2_im / ((double)n_samples * (n_samples - 1));\n";
    code = code <> "        }\n\n";

    code = code <> "        results[kp] = {total_re, total_im,\n";
    code = code <> "                       std::sqrt(total_var_re), std::sqrt(total_var_im)};\n";

    code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
    code = code <> "        if (kp == 0) {\n";
    code = code <> "            std::cerr << \"KP 0 total: (\" << total_re << \", \" << total_im\n";
    code = code <> "                      << \") +/- (\" << std::sqrt(total_var_re) << \", \"\n";
    code = code <> "                      << std::sqrt(total_var_im) << \")\" << std::endl;\n";
    code = code <> "        }\n";
    code = code <> "#endif\n";

    code = code <> "    }\n\n";
  ];

  code = code <> "    // Write results\n";
  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) {\n";
  code = code <> "        std::cerr << \"Cannot open \" << output_file << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        fout << results[kp][0] << \" \" << results[kp][1] << \" \"\n";
  code = code <> "             << results[kp][2] << \" \" << results[kp][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n\n";

  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points.\" << std::endl;\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";

  Export[outputFile, code, "Text"];

  Print["Generated C++ Monte Carlo code: ", outputFile];
  Print["  ", nConvergent, " convergent sectors"];
  Print["  Total: ", nIntegrands, " integrand functions"];

  Module[{codeStr, badPatterns, warnings},
    codeStr = code;
    badPatterns = {"Sin[", "Cos[", "Sqrt[", "Plus[", "Times[",
                   "Power[", "Rule[", "List["};
    warnings = Select[badPatterns, StringContainsQ[codeStr, #] &];
    If[Length[warnings] > 0,
      Print["WARNING: unresolved Mathematica symbols in C++ output: ",
            warnings]
    ];
  ];

  <|"Code" -> code, "OutputFile" -> outputFile,
    "NConvergent" -> nConvergent,
    "NTotal" -> nIntegrands,
    "Dimensions" -> integrandDims,
    "Integrator" -> integrator,
    "NeedsCuba" -> isVegas|>
];

(* --------------------------------------------------------------------------
   Private: locate a CUBA installation (header + static/shared lib).
   Returns the prefix dir (a string) or $Failed.  Logic is kept identical to
   EXAMPLES/cuba_common.wl`cubaFindPrefix on purpose (that file is standalone
   for the raw cross-check) -- do not let the two diverge.
   -------------------------------------------------------------------------- *)

findCubaPrefix[] := SelectFirst[
  {"/opt/homebrew", "/usr/local", "/usr"},
  FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
  (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &,
  $Failed];

(* --------------------------------------------------------------------------
   computeFanScaled: robust normal-fan computation.

   The normal fan of the Newton polytope is SCALE-INVARIANT, but the packaged
   fan code needs an integer interior lattice point; thin lattice simplices like
   conv{0, e_i} lack one for ambient dimension >= 4, so a direct
   ComputeDecomposition there leaks $Failed into the Polymake input and fails.
   Computing the fan from a scaled copy K*verts (which gains the interior point
   (1,...,1) once K > dim) yields the identical fan.  We try the cheap unscaled
   call first, then escalate K.  Returns {dualVertices, simplexList} or $Failed.
   (Mirrors TEST/gen_sectors.wl`computeFanRobust, which validated this against
   direct NIntegrate.)
   -------------------------------------------------------------------------- *)
computeFanScaled[verts_List] := Module[{nn, fd},
  nn = Length[First[verts]];
  Do[
    fd = Quiet[ComputeDecomposition[K*verts, "ShowProgress" -> False],
               TropicalFan::polymake];
    If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] &&
       Length[fd[[1]]] > 0 && Length[fd[[2]]] > 0,
      Return[fd, Module]],
    {K, {1, nn + 2, 2 nn + 4, 6 nn + 6}}];
  $Failed
];

(* --------------------------------------------------------------------------
   CompileCpp: Compile generated C++ code

   The optional 4th positional argument cubaPrefix links the CUBA library:
   when it is a string, -I<prefix>/include is added to the flags and
   -L<prefix>/lib -lcuba to the link line (before -lm).  When it is None
   (default) the command is byte-identical to the pre-VEGAS 3-arg call, so
   the MC path is unaffected.
   -------------------------------------------------------------------------- *)

CompileCpp[sourceFile_String, outputBinary_String,
           debug_: False, cubaPrefix_: None] :=
Module[{compiler, flags, includeFlags, linkFlags, buildCmd, cmd, result},
  compiler = "g++";

  flags = If[debug,
    {"-std=c++17", "-O2", "-fopenmp", "-DTROPICAL_MC_DEBUG",
     "-Wall", "-Wextra"},
    {"-std=c++17", "-O3", "-fopenmp", "-Wall", "-Wextra"}
  ];

  (* CUBA flags only when a prefix is supplied; otherwise both are empty and
     the command below collapses to the original byte-identical form. *)
  includeFlags = If[StringQ[cubaPrefix],
    {"-I" <> FileNameJoin[{cubaPrefix, "include"}]}, {}];
  linkFlags = If[StringQ[cubaPrefix],
    {"-L" <> FileNameJoin[{cubaPrefix, "lib"}], "-lcuba"}, {}];

  (* source must precede -lcuba so the static libcuba.a resolves its symbols *)
  buildCmd[fl_] := {compiler, Sequence @@ fl, Sequence @@ includeFlags,
    "-o", outputBinary, sourceFile, Sequence @@ linkFlags, "-lm"};

  cmd = buildCmd[flags];

  Print["Compiling: ", StringRiffle[cmd, " "]];
  result = RunProcess[cmd];

  (* If compilation fails due to -fopenmp (e.g. macOS clang), retry without it.
     CUBA include/link flags are preserved across the retry. *)
  If[result["ExitCode"] != 0 &&
     StringContainsQ[result["StandardError"], "fopenmp"],
    Print["  OpenMP not supported, retrying without -fopenmp..."];
    flags = DeleteCases[flags, "-fopenmp"];
    cmd = buildCmd[flags];
    Print["Compiling: ", StringRiffle[cmd, " "]];
    result = RunProcess[cmd];
  ];

  If[result["ExitCode"] != 0,
    Print["Compilation FAILED:"];
    Print[result["StandardError"]];
    Return[$Failed]
  ];

  If[StringLength[StringTrim[result["StandardError"]]] > 0,
    Print["Compiler warnings:"];
    Print[result["StandardError"]]
  ];

  Print["Compilation successful: ", outputBinary];
  outputBinary
];


(* ============================================================================
   MODULE 4: EvaluateTropicalMC (Driver)
   ============================================================================ *)

Options[EvaluateTropicalMC] = {
  "NSamples"       -> 1000000,
  "NThreads"       -> Automatic,
  "RunChecks"      -> True,
  "PrecisionGoal"  -> 3,
  "WorkingDirectory" -> Automatic,
  "Verbose"        -> True,
  "LiftData"       -> None,
  "SeedBase"       -> 42,
  "Integrator"     -> "MC",        (* "MC" | "VEGAS" *)
  (* VEGAS tuning -- forwarded to GenerateCppMonteCarlo when VEGAS *)
  "VegasEpsRel"    -> Automatic,    (* Automatic => 10^-PrecisionGoal, clamped *)
  "VegasSeed"      -> 0,
  (* Automatic => dimension-aware sizing (see resolveVegasSizing).  The fixed
     low-dim defaults 1000/500/1000 are FAR too small in high dimension: VEGAS
     then converges to a CONFIDENTLY WRONG value with a deceptively tight error
     bar (e.g. ~45% off at n=8).  Automatic scales the per-iteration grid
     resolution up with the sector dimension; n<=4 is left exactly at the old
     defaults.  Explicit integer values are always respected. *)
  "VegasNStart"    -> Automatic,
  "VegasNIncrease" -> Automatic,
  "VegasNBatch"    -> Automatic,
  "VegasMinEval"   -> 0
};

(* --------------------------------------------------------------------------
   resolveVegasSizing: dimension-aware VEGAS grid sizing.

   The CUBA-Vegas per-iteration sample budget (NStart, growing by NIncrease)
   must populate an adaptive grid in `dim` dimensions.  The shipped low-dim
   defaults (1000/500/1000) leave the grid catastrophically under-resolved for
   dim >= ~6: VEGAS locks onto a wrong stratification and returns a wrong value
   with a tiny (lying) error bar.  This resolves an Automatic sizing option to
   a value that scales with the maximum sector dimension, while reproducing the
   historical defaults exactly for dim <= 4 (so low-dim results are unchanged).
   -------------------------------------------------------------------------- *)
resolveVegasSizing[optVal_, maxDim_Integer, base_Integer, kind_String] :=
  If[optVal =!= Automatic,
    Round[optVal],
    (* dim<=4: historical default; dim>=5: grow ~3x per extra dimension *)
    Module[{scale = If[maxDim <= 4, 1, 3^(maxDim - 4)], v},
      v = Round[base * scale];
      Switch[kind,
        "NBatch", Min[v, 50000],   (* keep the per-call batch buffer modest *)
        _,        v
      ]
    ]
  ];

EvaluateTropicalMC[integrandSpec_Association, fanData_List,
                   kinematicPoints_List, OptionsPattern[]] :=
Module[
  {dualVertices, simplexList, n, nKP, nParams,
   allSectorData, convergentSectors, divergentSectors,
   cppFile, cppBinary, kinFile, resultFile,
   cppResult, mcResults, finalResults,
   runChecks, verbose, nSamples, nThreads,
   workDir, precGoal, seedBase,
   liftData, isLifted, liftedSpec, originalSpec, emptyDomainCount,
   integrator, isVegas, cubaPrefix, vegasEpsRel, vegasSeed,
   vegasNStart, vegasNIncrease, vegasNBatch, vegasMinEval},

  runChecks  = OptionValue["RunChecks"];
  verbose    = OptionValue["Verbose"];
  nSamples   = OptionValue["NSamples"];
  nThreads   = OptionValue["NThreads"];
  workDir    = OptionValue["WorkingDirectory"];
  precGoal   = OptionValue["PrecisionGoal"];
  seedBase   = OptionValue["SeedBase"];
  liftData   = OptionValue["LiftData"];
  isLifted   = (liftData =!= None);
  liftedSpec = integrandSpec;   (* when lifted, first arg IS the lifted spec *)
  originalSpec = If[isLifted, liftData["OriginalSpec"], integrandSpec];

  (* --- Integrator selection (default "MC" => path identical to before) --- *)
  integrator = OptionValue["Integrator"];
  If[!MatchQ[integrator, "MC" | "VEGAS"],
    Message[TropicalEval::badintegrator, integrator];
    Return[$Failed]
  ];
  isVegas = (integrator === "VEGAS");

  (* Resolve VEGAS tuning.  VegasEpsRel == Automatic -> 10^-PrecisionGoal,
     clamped to [1e-12, 1e-2]; pass a concrete number to the generator. *)
  vegasEpsRel    = OptionValue["VegasEpsRel"];
  If[vegasEpsRel === Automatic,
    vegasEpsRel = Clip[N[10^(-precGoal)], {1.*^-12, 1.*^-2}]];
  vegasSeed      = OptionValue["VegasSeed"];
  vegasNStart    = OptionValue["VegasNStart"];
  vegasNIncrease = OptionValue["VegasNIncrease"];
  vegasNBatch    = OptionValue["VegasNBatch"];
  vegasMinEval   = OptionValue["VegasMinEval"];

  (* VEGAS needs CUBA: locate it up front so we fail BEFORE codegen if absent.
     The user explicitly asked for VEGAS -- do NOT silently fall back to MC. *)
  cubaPrefix = None;
  If[isVegas,
    cubaPrefix = findCubaPrefix[];
    If[cubaPrefix === $Failed,
      Message[TropicalEval::nocuba];
      Return[$Failed]
    ]
  ];

  If[workDir === Automatic,
    workDir = DirectoryName[$InputFileName];
    If[workDir === "" || !StringQ[workDir],
      workDir = Directory[]
    ];
    workDir = FileNameJoin[{workDir, "INTERFILES"}]
  ];
  If[!DirectoryQ[workDir], Quiet[CreateDirectory[workDir]]];

  {dualVertices, simplexList} = fanData;
  n      = Length[integrandSpec["Variables"]];
  nKP    = Length[kinematicPoints];
  nParams = Length[integrandSpec["KinematicSymbols"]];

  If[verbose,
    Print["EvaluateTropicalMC: ", Length[simplexList], " sectors, ",
          nKP, " kinematic points, ", n, " variables"]
  ];

  (* --- Step 1: Validate inputs --- *)
  If[!MatchQ[integrandSpec["RegulatorSymbol"], None | _Missing],
    Message[TropicalEval::noregulator];
    Return[$Failed]
  ];

  (* Lifted-mode fan dimension assertion *)
  If[isLifted,
    Module[{nOrig, fanDim},
      nOrig  = Length[liftData["OriginalSpec"]["Variables"]];
      fanDim = If[Length[dualVertices] > 0, Length[dualVertices[[1]]], 0];
      If[fanDim != nOrig + 1,
        Message[TropicalEval::liftfandim, fanDim, nOrig + 1];
        Return[$Failed]
      ]
    ]
  ];

  If[nKP == 0,
    Print["ERROR: no kinematic points provided"];
    Return[$Failed]
  ];

  If[nParams > 0,
    If[!AllTrue[kinematicPoints,
                (ListQ[#] && Length[#] == nParams) &],
      Print["ERROR: kinematic points must have ", nParams, " parameters each"];
      Return[$Failed]
    ];
    If[!AllTrue[Flatten[kinematicPoints], (NumericQ[#] && Im[#] == 0 &&
                                          Abs[#] < 10^15) &],
      Print["ERROR: kinematic values must be real and finite"];
      Return[$Failed]
    ];
  ];

  (* --- Step 2: Process all sectors --- *)
  If[verbose, Print["Processing ", Length[simplexList], " sectors..."]];

  emptyDomainCount = 0;

  If[isLifted,
    (* Lifted mode: route through ProcessSectorLifted *)
    allSectorData = {};
    Do[
      Module[{sd},
        sd = ProcessSectorLifted[liftedSpec, dualVertices,
                                 simplexList[[s]], s, liftData,
                                 "Verbose" -> verbose];
        Which[
          (* $Failed from ProcessSectorLifted -> abort the whole call *)
          sd === $Failed,
            Return[$Failed, Module],
          (* EmptyDomain -> drop it, count it *)
          AssociationQ[sd] && KeyExistsQ[sd, "EmptyDomain"] && sd["EmptyDomain"],
            emptyDomainCount++;
            (* do NOT append to allSectorData — these contribute 0 *)
            ,
          (* Normal convergent sector *)
          True,
            AppendTo[allSectorData, sd]
        ]
      ],
      {s, Length[simplexList]}
    ];
    (* In lifted mode there are never "divergent sectors" — ProcessSectorLifted
       either errors via liftdivergent or returns a convergent sector. *)
    convergentSectors = allSectorData;
    divergentSectors  = {};   (* always empty in lifted mode *)
    If[verbose,
      Print["  ", Length[convergentSectors], " convergent sectors, ",
            emptyDomainCount, " EmptyDomain sectors dropped"]
    ];
    ,
    (* Standard (unlifted) mode: unchanged behavior *)
    allSectorData = Table[
      ProcessSector[integrandSpec, dualVertices,
                    simplexList[[s]], s, "Verbose" -> verbose],
      {s, Length[simplexList]}
    ];

    If[Count[allSectorData, _Association] != Length[simplexList],
      Print["WARNING: ", Count[allSectorData, $Failed],
            " sectors failed to process"];
    ];

    convergentSectors = Select[allSectorData,
      (AssociationQ[#] && !#["IsDivergent"]) &];
    divergentSectors  = Select[allSectorData,
      (AssociationQ[#] && #["IsDivergent"]) &];

    If[verbose,
      Print["  ", Length[convergentSectors], " convergent, ",
            Length[divergentSectors], " divergent sectors"]
    ]
  ];

  (* --- Step 3: Validation checks --- *)
  If[runChecks && Length[convergentSectors] > 0 && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[3, nKP]];
      Do[
        kinRules = Thread[
          originalSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        If[verbose,
          Print["Validating decomposition at kinematic point ", i, "..."]
        ];
        Module[{vr},
          If[isLifted,
            vr = Quiet@ValidateLiftedDecomposition[
              originalSpec, liftedSpec, fanData, liftData, kinRules, precGoal
            ],
            vr = Quiet@ValidateDecomposition[integrandSpec, fanData,
                                             kinRules, precGoal]
          ];
          If[AssociationQ[vr],
            Print["  Point ", i, ": rel error = ", vr["RelativeError"]]
          ];
        ];,
        {i, Length[testKP]}
      ];
    ]
  ];

  (* --- Step 4: Error on divergent sectors --- *)
  If[Length[divergentSectors] > 0,
    Message[TropicalEval::divergentinput,
            divergentSectors[[All, "ConeIndex"]],
            divergentSectors[[All, "NewExponents"]]];
    Return[$Failed]
  ];

  (* --- Step 5b: resolve dimension-aware VEGAS sizing (Automatic) --- *)
  If[isVegas,
    Module[{maxSecDim},
      maxSecDim = If[Length[convergentSectors] > 0,
        Max[#["Dimension"] & /@ convergentSectors], 1];
      vegasNStart    = resolveVegasSizing[vegasNStart,    maxSecDim, 1000, "NStart"];
      vegasNIncrease = resolveVegasSizing[vegasNIncrease, maxSecDim,  500, "NIncrease"];
      vegasNBatch    = resolveVegasSizing[vegasNBatch,    maxSecDim, 1000, "NBatch"];
      vegasMinEval   = Round[vegasMinEval];
      If[verbose && maxSecDim >= 5,
        Print["  VEGAS sizing (max sector dim ", maxSecDim, "): NStart=", vegasNStart,
              ", NIncrease=", vegasNIncrease, ", NBatch=", vegasNBatch]];
      (* The maxeval budget (n_samples) must allow several refinement iterations
         on top of the (now larger) per-iteration grid, or VEGAS stops before it
         has adapted.  Warn -- do NOT silently override the user's budget. *)
      If[maxSecDim >= 5 && nSamples < 20 vegasNStart,
        Message[TropicalEval::vegasbudget, maxSecDim, nSamples, vegasNStart,
                20 vegasNStart]]
    ],
    vegasNStart    = If[vegasNStart    === Automatic, 1000, Round[vegasNStart]];
    vegasNIncrease = If[vegasNIncrease === Automatic,  500, Round[vegasNIncrease]];
    vegasNBatch    = If[vegasNBatch    === Automatic, 1000, Round[vegasNBatch]];
    vegasMinEval   = Round[vegasMinEval]
  ];

  (* --- Step 6: Generate C++ code --- *)
  cppFile    = FileNameJoin[{workDir, "tropical_mc_generated.cpp"}];
  cppBinary  = FileNameJoin[{workDir, "tropical_mc"}];
  kinFile    = FileNameJoin[{workDir, "kinematic_data.txt"}];
  resultFile = FileNameJoin[{workDir, "mc_results.txt"}];

  cppResult = GenerateCppMonteCarlo[
    convergentSectors,
    {},
    integrandSpec, cppFile,
    "NSamples" -> nSamples,
    "SeedBase" -> seedBase,
    "Integrator" -> integrator,
    "VegasEpsRel" -> vegasEpsRel,
    "VegasSeed" -> vegasSeed,
    "VegasNStart" -> vegasNStart,
    "VegasNIncrease" -> vegasNIncrease,
    "VegasNBatch" -> vegasNBatch,
    "VegasMinEval" -> vegasMinEval,
    "ImagExponents" -> If[isLifted, Lookup[liftData, "ImagExponents", None], None]
  ];

  If[!AssociationQ[cppResult],
    Print["ERROR: C++ code generation failed"];
    Return[$Failed]
  ];

  (* --- Step 7: Debug compile and test --- *)
  If[runChecks,
    Module[{dbgBinary, dbgResult, dbgKinFile},
      dbgBinary  = FileNameJoin[{workDir, "tropical_mc_dbg"}];
      dbgKinFile = FileNameJoin[{workDir, "kinematic_data_dbg.txt"}];

      If[CompileCpp[cppFile, dbgBinary, True, cubaPrefix] =!= $Failed,
        Module[{testKP},
          testKP = Take[kinematicPoints, Min[5, nKP]];
          Export[dbgKinFile,
            StringRiffle[
              StringRiffle[ToString[CForm[#]] & /@ #, " "] & /@ testKP,
              "\n"
            ] <> "\n",
            "Text"
          ];

          dbgResult = RunProcess[{dbgBinary, dbgKinFile,
            FileNameJoin[{workDir, "mc_results_dbg.txt"}],
            "100000", "2"}];

          If[dbgResult["ExitCode"] == 0,
            Print["Debug run successful"];
            If[StringLength[StringTrim[dbgResult["StandardError"]]] > 0,
              Print["Debug output:\n", dbgResult["StandardError"]]
            ];,
            Print["Debug run FAILED:"];
            Print[dbgResult["StandardError"]];
          ];
        ];
      ];
    ]
  ];

  (* --- Step 8: Write kinematic data --- *)
  If[nParams > 0,
    Export[kinFile,
      StringRiffle[
        StringRiffle[
          ToString[CForm[#]] & /@ #, " "
        ] & /@ kinematicPoints,
        "\n"
      ] <> "\n",
      "Text"
    ];,
    Export[kinFile,
      StringRiffle[
        Table["0", {nKP}], "\n"
      ] <> "\n",
      "Text"
    ];
  ];

  (* --- Step 9: Release compile and run --- *)
  If[CompileCpp[cppFile, cppBinary, False, cubaPrefix] === $Failed,
    Print["ERROR: Release compilation failed"];
    Return[$Failed]
  ];

  Module[{nThreadsStr, result},
    nThreadsStr = If[nThreads === Automatic,
      ToString[$ProcessorCount],
      ToString[nThreads]
    ];

    If[verbose, Print["Running ", If[isVegas, "CUBA-Vegas", "Monte Carlo"],
                      " (", nSamples,
                      If[isVegas, " maxeval/sector", " samples"],
                      ", ", nThreadsStr, " threads)..."]];

    result = RunProcess[{cppBinary, kinFile, resultFile,
      ToString[nSamples], nThreadsStr, ToString[seedBase]}];

    If[result["ExitCode"] != 0,
      Print["ERROR: Monte Carlo execution failed:"];
      Print[result["StandardError"]];
      Return[$Failed]
    ];

    If[verbose && StringLength[StringTrim[result["StandardError"]]] > 0,
      Print[result["StandardError"]]
    ];
  ];

  (* --- Step 10: Read results --- *)
  Module[{rawResults, lines, parsed},
    rawResults = Import[resultFile, "Text"];
    If[rawResults === $Failed,
      Print["ERROR: cannot read results file"];
      Return[$Failed]
    ];

    lines = Select[StringSplit[rawResults, "\n"],
                   StringLength[StringTrim[#]] > 0 &];

    (* Use Read[StringToStream[...], Number] instead of ToExpression
       because C++ outputs scientific notation like 2.05e-05 which
       ToExpression misparses (treats 'e' as a symbol). *)
    parsed = Table[
      Read[StringToStream[#], Number] & /@ StringSplit[line],
      {line, lines}
    ];

    If[Length[parsed] != nKP,
      Print["WARNING: expected ", nKP, " result rows, got ",
            Length[parsed]];
    ];

    Module[{badRows},
      badRows = Select[Range[Length[parsed]],
        !AllTrue[parsed[[#]], NumericQ[#] && Abs[#] < 10^30 &] &
      ];
      If[Length[badRows] > 0,
        Print["WARNING: ", Length[badRows],
              " rows contain non-finite values"]
      ];
    ];

    mcResults = parsed;
  ];

  (* --- Step 11: Assemble final results --- *)
  finalResults = Table[
    If[i <= Length[mcResults] && Length[mcResults[[i]]] >= 4,
      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "Re"     -> mcResults[[i, 1]],
        "Im"     -> mcResults[[i, 2]],
        "ReErr"  -> mcResults[[i, 3]],
        "ImErr"  -> mcResults[[i, 4]]|>,
      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "Re" -> 0., "Im" -> 0., "ReErr" -> 0., "ImErr" -> 0.|>
    ],
    {i, nKP}
  ];

  (* --- Step 12: Error summary --- *)
  If[verbose,
    Module[{reErrs, imErrs, reMags},
      reErrs = #["ReErr"] & /@ finalResults;
      imErrs = #["ImErr"] & /@ finalResults;
      reMags = Abs[#["Re"]] & /@ finalResults;

      Print["\n=== Monte Carlo Error Summary ==="];
      Print["  Re errors: mean=", Mean[reErrs],
            " median=", Median[reErrs],
            " max=", Max[reErrs]];
      Print["  Im errors: mean=", Mean[imErrs],
            " median=", Median[imErrs],
            " max=", Max[imErrs]];

      Module[{badPts},
        badPts = Select[Range[nKP],
          (reMags[[#]] > 0 &&
           reErrs[[#]] / reMags[[#]] > 0.1) &
        ];
        If[Length[badPts] > 0,
          Print["  WARNING: ", Length[badPts],
                " points have Re error > 10% of result magnitude"]
        ];
      ];
    ]
  ];

  If[runChecks && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[3, nKP]];
      Print["\n=== NIntegrate Cross-Check ==="];
      Do[
        kinRules = Thread[
          originalSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        Module[{directResult, mcRe, relErr},
          If[isLifted,
            directResult = Quiet@ValidateLiftedDecomposition[
              originalSpec, liftedSpec, fanData, liftData, kinRules, precGoal
            ],
            directResult = Quiet@ValidateDecomposition[
              integrandSpec, fanData, kinRules, precGoal
            ]
          ];
          If[AssociationQ[directResult],
            mcRe = finalResults[[i, "Re"]] + I * finalResults[[i, "Im"]];
            relErr = Abs[(mcRe - directResult["DirectResult"]) /
                        directResult["DirectResult"]];
            Print["  Point ", i, ": MC = ", mcRe,
                  ", NIntegrate = ", directResult["DirectResult"],
                  ", rel err = ", relErr];
          ];
        ];,
        {i, Length[testKP]}
      ];
    ]
  ];

  <|"Results"           -> finalResults,
    "ConvergentSectors" -> Length[convergentSectors],
    "CppFile"           -> cppFile,
    "ResultFile"        -> resultFile|>
];



(* ============================================================================
   MODULE 4b: EvaluateTropicalMCLifted  (B7b)
   Convenience wrapper: auto-detect -> lift -> fan -> EvaluateTropicalMC
   ============================================================================ *)

Options[EvaluateTropicalMCLifted] = Join[
  {"LiftRules"            -> Automatic,
   "Threshold"            -> 1000,
   "AnchorRule"           -> "kStar",
   "BandEdgeGuard"        -> False,
   "FanData"              -> Automatic,
   (* "Reject": fire liftcomplexexponents and return $Failed when any B_k is complex.
      "SplitRealImag": decompose P^B = P^{Re(B)} * exp(i*Im(B)*log|P|); the real
      exponent drives lifting/sectors while the imaginary part contributes an
      oscillatory phase per MC sample.  Variance is not reduced by this splitting
      when Im(B) is large (rapid oscillation of the phase). *)
   "ComplexExponentMode"  -> "Reject"},
  Options[EvaluateTropicalMC]
];

EvaluateTropicalMCLifted[integrandSpec_Association, kinematicPoints_List,
                         opts: OptionsPattern[]] :=
Module[
  {liftRulesOpt, threshold, anchorRule, bandEdgeGuard, fanDataOpt,
   complexExpMode, hasComplexExps, imagExps, realSpec,
   rules, liftResult, liftedSpec, liftData,
   verts, fan, n,
   (* passthrough EvaluateTropicalMC opts *)
   passThroughOpts},

  liftRulesOpt  = OptionValue["LiftRules"];
  threshold     = OptionValue["Threshold"];
  anchorRule    = OptionValue["AnchorRule"];
  bandEdgeGuard = OptionValue["BandEdgeGuard"];
  fanDataOpt    = OptionValue["FanData"];

  (* Collect passthrough opts for EvaluateTropicalMC *)
  passThroughOpts = Sequence @@ FilterRules[{opts}, Options[EvaluateTropicalMC]];

  (* --- Determine lift rules --- *)
  If[liftRulesOpt === Automatic,
    (* SuggestedK (hence the anchor z0) is set by the "AnchorRule"/"BandEdgeGuard"
       options, defaulting to the kStar rule (Manual sec:anchor). *)
    rules = DetectExtremeCoefficients[integrandSpec, threshold,
              "AnchorRule" -> anchorRule, "BandEdgeGuard" -> bandEdgeGuard];
    If[rules === {},
      (* No extreme coefficients found: fall back to plain EvaluateTropicalMC *)
      Print["EvaluateTropicalMCLifted: no extreme coefficients detected \
(threshold=", threshold, "); falling back to plain EvaluateTropicalMC on \
the original spec."];
      Module[{origFan},
        origFan = fanDataOpt;
        If[origFan === Automatic,
          Module[{origVerts},
            origVerts = PolytopeVertices[
              (Times @@ integrandSpec["Polynomials"])^(-1),
              integrandSpec["Variables"]
            ];
            origFan = ComputeDecomposition[origVerts, "ShowProgress" -> False];
          ]
        ];
        Return[EvaluateTropicalMC[integrandSpec, origFan,
                                  kinematicPoints, passThroughOpts]]
      ]
    ];
    (* DetectExtremeCoefficients uses "SuggestedK" but LiftCoefficients
       expects "k" — remap here so the two functions interoperate. *)
    rules = Map[
      <|"PolyIndex"      -> #["PolyIndex"],
        "ExponentVector" -> #["ExponentVector"],
        "k"              -> #["SuggestedK"]|> &,
      rules
    ]
    ,
    rules = liftRulesOpt
  ];

  (* --- Detect and handle complex polynomial exponents (Fix 2/3) --- *)
  complexExpMode = OptionValue["ComplexExponentMode"];
  hasComplexExps = AnyTrue[N[integrandSpec["PolynomialExponents"]],
    (NumericQ[#] && Abs[Im[#]] >= 1*^-12) &];
  If[hasComplexExps && complexExpMode === "Reject",
    Message[TropicalEval::liftcomplexexponents,
            Count[N[integrandSpec["PolynomialExponents"]],
                  _?(Abs[Im[N[#]]] >= 1*^-12 &)],
            Length[integrandSpec["PolynomialExponents"]]];
    Return[$Failed]
  ];
  If[hasComplexExps && complexExpMode === "SplitRealImag",
    (* Split: real part drives lifting/sectors; imaginary part → phase per sample.
       Store imagExps in realSpec so LiftCoefficients forwards it to liftData. *)
    imagExps = Im[N[#]] & /@ integrandSpec["PolynomialExponents"];
    realSpec = Association[integrandSpec,
      "PolynomialExponents" -> (Re[N[#]] & /@ integrandSpec["PolynomialExponents"]),
      "ImagExponents"       -> imagExps]
    ,
    imagExps = ConstantArray[0., Length[integrandSpec["PolynomialExponents"]]];
    realSpec = integrandSpec
  ];

  (* --- Lift the integrand --- *)
  liftResult = LiftCoefficients[realSpec, rules];
  If[!AssociationQ[liftResult], Return[$Failed]];
  liftedSpec = liftResult["LiftedSpec"];
  (* Restore original spec (with complex exponents) in liftData so that
     ValidateLiftedDecomposition uses the correct reference integrand. *)
  liftData   = Association[liftResult["LiftData"],
    "OriginalSpec"  -> integrandSpec,
    "ImagExponents" -> imagExps
  ];
  n          = Length[integrandSpec["Variables"]];

  (* --- Build or use the (n+1)-dimensional fan --- *)
  Module[{liftedFan},
    If[fanDataOpt =!= Automatic,
      liftedFan = fanDataOpt
      ,
      (* Automatic: compute from lifted integrand.
         Wrap in Quiet so that polymake errors for degenerate (lower-
         dimensional) lifted polytopes are caught cleanly by the guard
         below rather than leaking as raw polymake messages. *)
      verts = Quiet[
        PolytopeVertices[
          (Times @@ liftedSpec["Polynomials"])^(-1),
          liftedSpec["Variables"]
        ],
        TropicalFan::polymake
      ];
      (* Use the scale-robust fan computation: a direct ComputeDecomposition
         leaks $Failed for n+1 >= 4 (thin lattice simplices lack an interior
         lattice point), which previously made the automatic path fail with a
         misleading liftdegenerate even for genuinely full-dimensional lifted
         polytopes.  computeFanScaled retries on a scaled copy (same fan). *)
      liftedFan = If[ListQ[verts], computeFanScaled[verts], $Failed];

      (* Degeneracy guard: check that fan was computed and simplices have
         the right length (n+2 rays for n+1 variables) *)
      If[!ListQ[liftedFan] || Length[liftedFan] < 2,
        Message[TropicalEval::liftdegenerate];
        Return[$Failed]
      ];
      Module[{dv, sl},
        {dv, sl} = liftedFan;
        If[Length[sl] == 0 ||
           Length[dv] == 0 ||
           !AllTrue[sl, Length[#] == n + 1 &],
          (* Simplices whose length != n+1 indicate a lower-dimensional polytope *)
          Message[TropicalEval::liftdegenerate];
          Return[$Failed]
        ]
      ]
    ];

    (* --- Call EvaluateTropicalMC with LiftData --- *)
    EvaluateTropicalMC[liftedSpec, liftedFan, kinematicPoints,
                       "LiftData" -> liftData, passThroughOpts]
  ]
];

(* --------------------------------------------------------------------------
   Package end
   -------------------------------------------------------------------------- *)

End[]

EndPackage[]
