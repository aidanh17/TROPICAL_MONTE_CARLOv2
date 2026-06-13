(* ============================================================================
   sandbox_lift_common.wl
   Standalone lifting mathematics for Phase 1 SANDBOX.
   Implements: liftSpec, deltaEliminate, sectorNIntegrate, sectorMC,
               unliftedSectorStats.

   Load AFTER Get["tropical_eval.wl"] (which loads tropical_fan.wl).
   ============================================================================ *)

(* ----------------------------------------------------------------
   liftSpec[spec, polyIdx, expVec, k]
   Replace the coefficient C of monomial expVec in polynomial polyIdx
   by z^k, keeping z0 = C^(1/k) exact.
   Augmented spec: aux variable x[n+1], MonomialExponents appended with 0.
   Returns {liftedSpec, liftData}.
   ---------------------------------------------------------------- *)
liftSpec[spec_Association, polyIdx_Integer, expVec_List, k_Integer] :=
Module[
  {polys, vars, monoExps, polyExps, n, poly, parsedPoly,
   matchPos, C0, z0, auxVar, auxIdx, newVars, newMonoExps,
   liftedPolys, liftedPoly, termToReplace, replacement,
   liftRules, liftData, liftedSpec},

  polys    = spec["Polynomials"];
  vars     = spec["Variables"];
  monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"];
  n        = Length[vars];
  auxIdx   = n + 1;
  auxVar   = x[auxIdx];

  (* Parse the target polynomial using the public TropicalEval function *)
  poly       = polys[[polyIdx]];
  parsedPoly = ParsePolynomial[poly, vars];

  (* Find the monomial matching expVec *)
  matchPos = Position[parsedPoly[[All, 2]], expVec];
  If[matchPos === {},
    Print["liftSpec: monomial with exponent ", expVec,
          " not found in polynomial ", polyIdx, ".  parsed=", parsedPoly];
    Return[$Failed]
  ];
  C0 = parsedPoly[[matchPos[[1, 1]], 1]];

  (* Anchor: z0 = C0^(1/k), exact *)
  z0 = If[IntegerQ[C0^(1/k)], C0^(1/k), Power[C0, 1/k]];

  (* Build the term to replace: C0 * prod_i x[i]^expVec[i]
     (skip factors where exponent = 0) *)
  termToReplace = Times @@ Flatten[{C0,
    MapThread[If[#2 == 0, Nothing, Power[#1, #2]]&, {vars, expVec}]
  }];

  (* Replacement monomial: auxVar^k * prod_i x[i]^expVec[i] *)
  replacement = Times @@ Flatten[{auxVar^k,
    MapThread[If[#2 == 0, Nothing, Power[#1, #2]]&, {vars, expVec}]
  }];

  liftedPoly  = poly - termToReplace + replacement;

  newVars     = Append[vars, auxVar];
  newMonoExps = Append[monoExps, 0];
  liftedPolys = ReplacePart[polys, polyIdx -> liftedPoly];

  liftedSpec = <|
    "Polynomials"         -> liftedPolys,
    "MonomialExponents"   -> newMonoExps,
    "PolynomialExponents" -> polyExps,
    "Variables"           -> newVars,
    "KinematicSymbols"    -> spec["KinematicSymbols"],
    "RegulatorSymbol"     -> spec["RegulatorSymbol"]
  |>;

  liftRules = {<|"PolyIndex"     -> polyIdx,
                 "ExponentVector"-> expVec,
                 "k"             -> k|>};

  liftData = <|
    "z0"          -> z0,
    "k"           -> k,
    "AuxIndex"    -> auxIdx,
    "AuxVariable" -> auxVar,
    "Rules"       -> liftRules,
    "OriginalSpec"-> spec
  |>;

  {liftedSpec, liftData}
];


(* ----------------------------------------------------------------
   deltaEliminate[sdAug, liftData]
   Resolve the delta(z - z0) in sector sdAug (output of ProcessSector
   on the augmented (n+1)-dim spec).
   Returns a lifted sector association,
   or <|"EmptyDomain"->True, "ConeIndex"->...|>,
   or $Failed with a printed diagnostic.

   Works on BOTH branches of ProcessSector (IsDivergent True or False);
   reads ClearedPolys, NewExponents, DetM, RayMatrix.  Per §3.4 of plan:
   IsDivergent is MEANINGLESS for the lifted path — do not branch on it.
   ---------------------------------------------------------------- *)
deltaEliminate[sdAug_Association, liftData_Association] :=
Module[
  {z0, auxIdx, coneIndex,
   a, clearedPolys, detM, mMatrix, polyExps,
   n1, n, mVec,
   candidates, bestPivot,
   tryPivot, result},

  z0        = liftData["z0"];
  auxIdx    = liftData["AuxIndex"];
  coneIndex = sdAug["ConeIndex"];

  (* Extract from sdAug — both branches of ProcessSector expose these *)
  a            = sdAug["NewExponents"];      (* length n+1 *)
  clearedPolys = sdAug["ClearedPolys"];      (* {poly1, poly2, ...} each = {{c,e},...} *)
  detM         = sdAug["DetM"];
  mMatrix      = sdAug["RayMatrix"];
  polyExps     = sdAug["PolynomialExponents"];

  n1 = Length[a];      (* n+1 *)
  n  = n1 - 1;         (* original dimension *)

  (* z-row: row auxIdx of M, i.e. mMatrix[[auxIdx]] *)
  mVec = mMatrix[[auxIdx]];    (* length n+1 *)

  (* ----------------------------------------------------------------
     Per-pivot computation (§3.2 - §3.3):
     Given pivot p, compute: monomial substitution -> re-clear -> atilde.
     Returns an association with all relevant data, or $Failed if
     a purely symbolic issue arises.
     ---------------------------------------------------------------- *)
  tryPivot[p_] := Module[
    {mp, ap, remainIdx, mOtherVec, aOtherVec,
     subPolys, reclearMin, newPolyList, atildeRaw, hasConst},

    mp        = mVec[[p]];
    ap        = a[[p]];
    remainIdx = DeleteCases[Range[n1], p];    (* indices of the n remaining vars *)
    mOtherVec = mVec[[remainIdx]];            (* m_j for j != p *)
    aOtherVec = a[[remainIdx]];               (* a_j for j != p *)

    (* §3.2 step 4: monomial substitution
       {c, e} -> {c * z0^(e_p/mp),  e_j - e_p*m_j/mp  for j in remainIdx}
       Use Function to avoid name clashes with local vars. *)
    subPolys = Table[
      Function[{cpoly},
        Map[Function[{cmono},
          Module[{cep = cmono[[2, p]]},
            {cmono[[1]] * z0^(cep / mp),
             Table[cmono[[2, remainIdx[[j]]]] - cep * mOtherVec[[j]] / mp, {j, n}]}
          ]
        ], cpoly]
      ][clearedPolys[[k]]],
      {k, Length[clearedPolys]}
    ];

    (* §3.3 Re-clear: componentwise minima *)
    reclearMin = Table[
      Table[Min[#[[2, j]] & /@ subPolys[[k]]], {j, n}],
      {k, Length[subPolys]}
    ];

    (* Subtract minima from every monomial *)
    newPolyList = Table[
      Map[Function[{cmono}, {cmono[[1]], cmono[[2]] - reclearMin[[k]]}],
          subPolys[[k]]],
      {k, Length[subPolys]}
    ];

    (* §3.2 step 5: atilde_j = a_j - a_p * m_j/m_p  for j in remainIdx *)
    atildeRaw = Table[
      aOtherVec[[j]] - ap * mOtherVec[[j]] / mp,
      {j, n}
    ];

    (* Absorb re-clear minima: atilde_j += Sum_k B_k * dtilde_{k,j} *)
    atildeRaw = atildeRaw +
      Total[Table[polyExps[[k]] * reclearMin[[k]], {k, Length[polyExps]}]];

    (* HasConstantTerm: does each Q_k have a monomial with all-zero exponents? *)
    hasConst = And @@ Table[
      AnyTrue[newPolyList[[k]], (#[[2]] === ConstantArray[0, n])&],
      {k, Length[newPolyList]}
    ];

    <|"pivot"      -> p,
      "mp"         -> mp,
      "ap"         -> ap,
      "remainIdx"  -> remainIdx,
      "mOtherVec"  -> mOtherVec,
      "atilde"     -> atildeRaw,
      "newPolys"   -> newPolyList,
      "hasConst"   -> hasConst,
      "reclearMin" -> reclearMin
    |>
  ];

  (* ---- Pivot search: §3.5 ---- *)
  candidates = {};
  Do[
    If[mVec[[p]] != 0,
      Module[{res = tryPivot[p]},
        (* Admissibility: all atilde_j real numeric (Im < 10^-12) and > 0 *)
        If[And @@ Table[
             With[{av = N[res["atilde"][[j]]]},
               NumericQ[av] && Abs[Im[av]] < 10^-12 && Re[av] > 0
             ],
             {j, n}
           ],
          AppendTo[candidates, res]
        ]
      ]
    ],
    {p, n1}
  ];

  If[candidates === {},
    (* Full diagnostics *)
    Print["deltaEliminate: cone ", coneIndex, " — no admissible pivot."];
    Print["  z-row m = ", mVec];
    Print["  a = ", a];
    Print["  Per-pivot atilde:"];
    Table[
      If[mVec[[p]] != 0,
        Module[{res = tryPivot[p]},
          Print["    p=", p, "  mp=", mVec[[p]],
                "  atilde=", N[res["atilde"]]]
        ]
      ],
      {p, n1}
    ];
    Return[$Failed]
  ];

  (* §3.5 ranking: (1) |mp|=1, (2) hasConst, (3) max min Re[atilde] *)
  candidates = SortBy[candidates,
    {-Boole[Abs[#["mp"]] == 1],
     -Boole[#["hasConst"]],
     -Min[Re[N[#["atilde"]]]]} &
  ];
  bestPivot = candidates[[1]];

  (* ---- Domain constraint classification: §3.4 ---- *)
  Module[
    {p, mp, mOtherVec, atilde, logZ0, allOtherZero,
     domainClass, flatPolys, prefactorBase, prefactor, ap},

    p           = bestPivot["pivot"];
    mp          = bestPivot["mp"];
    mOtherVec   = bestPivot["mOtherVec"];   (* length n *)
    atilde      = bestPivot["atilde"];       (* length n *)
    ap          = bestPivot["ap"];
    logZ0       = Log[z0];

    (* Check constant-root case: all m_{j != p} = 0 *)
    allOtherZero = And @@ (# == 0 & /@ mOtherVec);

    If[allOtherZero,
      (* Constant root: y_p* = z0^(1/mp) *)
      If[N[z0^(1/mp)] > 1,
        Return[<|"EmptyDomain" -> True, "ConeIndex" -> coneIndex|>]
      ];
      domainClass = None;
      ,
      (* General case *)
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
            "IndicatorCoeffs" -> Table[mOtherVec[[j]] / atilde[[j]], {j, n}]
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
            "IndicatorCoeffs" -> Table[mOtherVec[[j]] / atilde[[j]], {j, n}]
          |>
        ]
      ]
    ];

    (* ---- Flatten: §3.4 ---- *)
    flatPolys = Table[
      Map[Function[{fm}, {fm[[1]], fm[[2]] / atilde}],
          bestPivot["newPolys"][[k]]],
      {k, Length[bestPivot["newPolys"]]}
    ];

    (* §3.2 step 6 + §3.4:
       prefactor = (|detM| / |mp|) * z0^(ap/mp - 1) / Times@@atilde *)
    prefactorBase = (Abs[detM] / Abs[mp]) * z0^(ap / mp - 1);
    prefactor     = prefactorBase / (Times @@ atilde);

    <|"FlattenedPolys"      -> flatPolys,
      "ClearedPolys"        -> bestPivot["newPolys"],
      "Prefactor"           -> prefactor,
      "Dimension"           -> n,
      "PolynomialExponents" -> polyExps,
      "DomainConstraint"    -> domainClass,
      "PivotIndex"          -> p,
      "ZRow"                -> mVec,
      "ATilde"              -> atilde,
      "HasConstantTerm"     -> bestPivot["hasConst"],
      "ConeIndex"           -> coneIndex,
      "DetM"                -> detM,
      "RayMatrix"           -> mMatrix
    |>
  ]
];


(* ----------------------------------------------------------------
   sectorNIntegrate[liftedSD, kinRules, pg]
   NIntegrate the lifted sector on [0,1]^n using the log-exp form
   from tropical_eval.wl:469-496, with optional domain indicator.
   ---------------------------------------------------------------- *)
sectorNIntegrate[liftedSD_Association, kinRules_List, pg_Integer] :=
Module[
  {flatPolys, polyExps, pf, dim, yVars, domConstr,
   polyVals, integrand, logYpStar, result,
   logZ0, mp, indCoeffs},

  (* EmptyDomain sectors contribute exactly 0 *)
  If[KeyExistsQ[liftedSD, "EmptyDomain"] && TrueQ[liftedSD["EmptyDomain"]],
    Return[0]
  ];

  flatPolys = liftedSD["FlattenedPolys"];
  polyExps  = liftedSD["PolynomialExponents"] /. kinRules;
  pf        = liftedSD["Prefactor"] /. kinRules;
  dim       = liftedSD["Dimension"];
  domConstr = liftedSD["DomainConstraint"];
  yVars     = Table[Unique["yv"], {dim}];

  (* Log-exp pattern from tropical_eval.wl:469-496 *)
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
      {polyVals, polyExps}
    ];

  (* Domain indicator (§3.4) *)
  If[domConstr =!= None,
    logZ0     = domConstr["LogZ0"];
    mp        = domConstr["MP"];
    indCoeffs = domConstr["IndicatorCoeffs"] /. kinRules;
    logYpStar = (logZ0 - Total[indCoeffs * Log /@ yVars]) / mp;
    integrand = integrand * Boole[logYpStar <= 0];
  ];

  result = Quiet@NIntegrate[
    integrand,
    Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
    MaxRecursion -> 15,
    PrecisionGoal -> pg,
    Method -> "GlobalAdaptive"
  ];
  result
];


(* ----------------------------------------------------------------
   sectorMC[liftedSD, kinRules, nSamples]
   Plain uniform MC on [0,1]^n (full range, per §5.1).
   Returns <|"Mean","StdErr","FeasibleFraction","MinMag","MaxMag"|>.
   ---------------------------------------------------------------- *)
sectorMC[liftedSD_Association, kinRules_List, nSamples_Integer] :=
Module[
  {flatPolys, polyExps, pf, dim, domConstr,
   logZ0, mp, indCoeffs,
   vals, feasCount, minMag, maxMag,
   yv, logYv, polyVals, ival, logYpStar, feasible,
   mean, stderr, feasFrac},

  If[KeyExistsQ[liftedSD, "EmptyDomain"] && TrueQ[liftedSD["EmptyDomain"]],
    Return[<|"Mean" -> 0, "StdErr" -> 0, "FeasibleFraction" -> 0,
             "MinMag" -> 0, "MaxMag" -> 0, "SigmaPerSample" -> 0|>]
  ];

  flatPolys = liftedSD["FlattenedPolys"];
  polyExps  = N[liftedSD["PolynomialExponents"] /. kinRules];
  pf        = N[liftedSD["Prefactor"] /. kinRules];
  dim       = liftedSD["Dimension"];
  domConstr = liftedSD["DomainConstraint"];

  If[domConstr =!= None,
    logZ0     = N[domConstr["LogZ0"]];
    mp        = domConstr["MP"];
    indCoeffs = N[domConstr["IndicatorCoeffs"] /. kinRules];
  ];

  vals      = {};
  feasCount = 0;
  minMag    = Infinity;
  maxMag    = 0;

  Do[
    yv    = RandomReal[{0, 1}, dim];
    logYv = Log /@ yv;

    feasible = True;
    If[domConstr =!= None,
      logYpStar = (logZ0 - Total[indCoeffs * logYv]) / mp;
      If[logYpStar > 0, feasible = False]
    ];

    If[feasible,
      polyVals = Table[
        Total[
          Table[
            N[mono[[1]]] * Exp[Total[N[mono[[2]]] * logYv]],
            {mono, flatPolys[[j]]}
          ]
        ],
        {j, Length[flatPolys]}
      ];
      ival = pf * Times @@ MapThread[
               Function[{pv, be}, If[pv > 0, Exp[be * Log[pv]], 0.]],
               {polyVals, polyExps}
             ];
      AppendTo[vals, ival];
      feasCount++;
      minMag = Min[minMag, Abs[ival]];
      maxMag = Max[maxMag, Abs[ival]];
    ],
    {nSamples}
  ];

  feasFrac = N[feasCount / nSamples];

  If[feasCount == 0,
    Return[<|"Mean" -> 0, "StdErr" -> 0, "FeasibleFraction" -> 0,
             "MinMag" -> 0, "MaxMag" -> 0, "SigmaPerSample" -> 0|>]
  ];

  (* Compute per-sample sigma of g = f * Indicator(feasible).
     Infeasible draws contribute g = 0.
     mean1 = Mean[f | feasible], mean2 = Mean[f^2 | feasible]
     sigmaPerSample = Sqrt[ feasFrac * mean2 - (feasFrac * mean1)^2 ] *)
  Module[{absVals, mean1, mean2, sigmaPS},
    absVals = Abs /@ vals;   (* vals are real here; Abs for complex-safety *)
    mean1   = N[Mean[vals]];
    mean2   = N[Mean[vals^2]];
    sigmaPS = N[Sqrt[Max[0, feasFrac * mean2 - (feasFrac * mean1)^2]]];

    <|"Mean"              -> mean1 * feasFrac,
      "StdErr"            -> sigmaPS / Sqrt[nSamples],
      "FeasibleFraction"  -> feasFrac,
      "MinMag"            -> minMag,
      "MaxMag"            -> maxMag,
      "SigmaPerSample"    -> sigmaPS
    |>
  ]
];


(* ----------------------------------------------------------------
   unliftedSectorStats[spec, fanData, nSamples]
   ProcessSector on original spec per cone; per-sector magnitude
   stats + MC sigma (the "before" column).
   ---------------------------------------------------------------- *)
unliftedSectorStats[spec_Association, fanData_List, nSamples_Integer] :=
Module[
  {dualVertices, simplexList, nSectors, polyExps},

  {dualVertices, simplexList} = fanData;
  nSectors = Length[simplexList];
  polyExps = N[spec["PolynomialExponents"]];

  Table[
    Module[{sd, flatPolys, pf, dim, vals, yv, logYv, pvs, ival,
            mean, stderr},
      sd = ProcessSector[spec, dualVertices, simplexList[[s]], s];
      If[sd === $Failed || TrueQ[sd["IsDivergent"]],
        <|"Sector" -> s, "IsDivergent" -> True,
          "Min" -> Null, "Max" -> Null, "Sigma" -> Null|>
        ,
        flatPolys = sd["FlattenedPolys"];
        pf        = N[sd["Prefactor"]];
        dim       = sd["Dimension"];
        vals = Table[
          yv    = RandomReal[{0, 1}, dim];
          logYv = Log /@ yv;
          pvs = Table[
            Total[Table[
              N[mono[[1]]] * Exp[Total[N[mono[[2]]] * logYv]],
              {mono, flatPolys[[j]]}
            ]],
            {j, Length[flatPolys]}
          ];
          pf * Times @@ MapThread[
            Function[{pv, be}, If[pv > 0, Exp[be * Log[pv]], 0.]],
            {pvs, polyExps}
          ],
          {nSamples}
        ];
        mean   = N[Mean[vals]];
        stderr = N[Sqrt[Variance[vals] / nSamples]];
        <|"Sector"      -> s,
          "IsDivergent" -> False,
          "Mean"        -> mean,
          "StdErr"      -> stderr,
          "Sigma"       -> N[StandardDeviation[vals]],
          "Min"         -> Min[Abs /@ vals],
          "Max"         -> Max[Abs /@ vals]
        |>
      ]
    ],
    {s, nSectors}
  ]
];


(* ----------------------------------------------------------------
   trueSigma[sectorLikeData, kinRules, pg]
   Compute exact per-sample sigma via NIntegrate:
     sigma^2 = I2 - I1^2,  I1 = Int g,  I2 = Int g^2
   where g = integrand (with Boole indicator when DomainConstraint present).
   Works for both lifted sectors (from deltaEliminate) and unlifted
   ProcessSector outputs (DomainConstraint key absent => None).
   EmptyDomain sectors return all zeros.
   I2 divergence (non-integrable singularity) is detected and flagged.
   ---------------------------------------------------------------- *)
trueSigma[sectorLikeData_Association, kinRules_List, pg_Integer] :=
Module[
  {flatPolys, polyExps, pf, dim, yVars, domConstr,
   polyVals, integrand, logYpStar, logZ0, mp, indCoeffs,
   i1, i2, i2result, i2converged},

  (* EmptyDomain sectors contribute zero *)
  If[KeyExistsQ[sectorLikeData, "EmptyDomain"] &&
     TrueQ[sectorLikeData["EmptyDomain"]],
    Return[<|"I1" -> 0, "I2" -> 0, "Sigma" -> 0, "I2Converged" -> True|>]
  ];

  flatPolys = sectorLikeData["FlattenedPolys"];
  polyExps  = sectorLikeData["PolynomialExponents"] /. kinRules;
  pf        = sectorLikeData["Prefactor"] /. kinRules;
  dim       = sectorLikeData["Dimension"];
  domConstr = Lookup[sectorLikeData, "DomainConstraint", None];
  yVars     = Table[Unique["tsyv"], {dim}];

  (* Build log-exp integrand (same pattern as sectorNIntegrate) *)
  polyVals = Table[
    Total[
      Table[
        Module[{coeff, alphas},
          coeff  = mono[[1]] /. kinRules;
          alphas = mono[[2]] /. kinRules;
          coeff * Exp[Total[alphas * (Log /@ yVars)]]
        ],
        {mono, flatPolys[[j]]}
      ]
    ],
    {j, Length[flatPolys]}
  ];

  integrand = pf *
    Times @@ MapThread[
      Function[{pv, be}, Exp[be * Log[pv]]],
      {polyVals, polyExps}
    ];

  (* Domain indicator *)
  If[domConstr =!= None,
    logZ0     = domConstr["LogZ0"];
    mp        = domConstr["MP"];
    indCoeffs = domConstr["IndicatorCoeffs"] /. kinRules;
    logYpStar = (logZ0 - Total[indCoeffs * (Log /@ yVars)]) / mp;
    integrand = integrand * Boole[logYpStar <= 0];
  ];

  (* I1 = Integral of g *)
  i1 = Quiet@NIntegrate[
    integrand,
    Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
    MaxRecursion -> 12,
    PrecisionGoal -> pg,
    Method -> "GlobalAdaptive"
  ];

  (* I2 = Integral of g^2; may diverge *)
  i2result = Quiet@Check[
    TimeConstrained[
      NIntegrate[
        integrand^2,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 12,
        PrecisionGoal -> pg,
        Method -> "GlobalAdaptive"
      ],
      120,
      $TimedOut
    ],
    $Failed
  ];

  (* Classify I2 convergence *)
  Which[
    i2result === $Failed || i2result === $TimedOut,
      i2converged = False;
      i2 = Infinity,
    !NumericQ[i2result],
      i2converged = False;
      i2 = Infinity,
    Abs[i2result] > 10^30,
      i2converged = False;
      i2 = Infinity,
    True,
      i2converged = True;
      i2 = i2result
  ];

  <|"I1"          -> i1,
    "I2"          -> i2,
    "Sigma"       -> If[i2converged, Sqrt[Max[0, i2 - i1^2]], Infinity],
    "I2Converged" -> i2converged
  |>
];


(* ----------------------------------------------------------------
   computePivotRow[sdAug, liftData, p]
   Compute the per-pivot data for pivot index p, reusing the same
   math as deltaEliminate's internal tryPivot.  Returns an association
   with fields: pivot, mp, atilde (numeric), allPositive, hasConst.
   Returns $Failed if mVec[[p]] == 0.
   ---------------------------------------------------------------- *)
computePivotRow[sdAug_Association, liftData_Association, p_Integer] :=
Module[
  {z0, auxIdx, a, clearedPolys, mMatrix, polyExps,
   n1, n, mVec, mp, ap, remainIdx, mOtherVec, aOtherVec,
   subPolys, reclearMin, newPolyList, atildeRaw, hasConst, allPos},

  z0        = liftData["z0"];
  auxIdx    = liftData["AuxIndex"];
  a         = sdAug["NewExponents"];
  clearedPolys = sdAug["ClearedPolys"];
  mMatrix   = sdAug["RayMatrix"];
  polyExps  = sdAug["PolynomialExponents"];
  n1        = Length[a];
  n         = n1 - 1;
  mVec      = mMatrix[[auxIdx]];

  mp = mVec[[p]];
  If[mp == 0, Return[$Failed]];

  ap        = a[[p]];
  remainIdx = DeleteCases[Range[n1], p];
  mOtherVec = mVec[[remainIdx]];
  aOtherVec = a[[remainIdx]];

  subPolys = Table[
    Function[{cpoly},
      Map[Function[{cmono},
        Module[{cep = cmono[[2, p]]},
          {cmono[[1]] * z0^(cep / mp),
           Table[cmono[[2, remainIdx[[j]]]] - cep * mOtherVec[[j]] / mp, {j, n}]}
        ]
      ], cpoly]
    ][clearedPolys[[k]]],
    {k, Length[clearedPolys]}
  ];

  reclearMin = Table[
    Table[Min[#[[2, j]] & /@ subPolys[[k]]], {j, n}],
    {k, Length[subPolys]}
  ];

  newPolyList = Table[
    Map[Function[{cmono}, {cmono[[1]], cmono[[2]] - reclearMin[[k]]}],
        subPolys[[k]]],
    {k, Length[subPolys]}
  ];

  atildeRaw = Table[
    aOtherVec[[j]] - ap * mOtherVec[[j]] / mp,
    {j, n}
  ];
  atildeRaw = atildeRaw +
    Total[Table[polyExps[[k]] * reclearMin[[k]], {k, Length[polyExps]}]];

  hasConst = And @@ Table[
    AnyTrue[newPolyList[[k]], (#[[2]] === ConstantArray[0, n])&],
    {k, Length[newPolyList]}
  ];

  allPos = And @@ Table[
    With[{av = N[atildeRaw[[j]]]},
      NumericQ[av] && Abs[Im[av]] < 10^-12 && Re[av] > 0
    ],
    {j, n}
  ];

  <|"pivot"      -> p,
    "mp"         -> mp,
    "atilde"     -> N[atildeRaw],
    "allPositive" -> allPos,
    "hasConst"   -> hasConst
  |>
];


(* ----------------------------------------------------------------
   pivotTable[sdAug, liftData]
   Print a per-pivot diagnostic table for every p with m_p != 0.
   ---------------------------------------------------------------- *)
pivotTable[sdAug_Association, liftData_Association] :=
Module[
  {auxIdx, mMatrix, mVec, n1},

  auxIdx  = liftData["AuxIndex"];
  mMatrix = sdAug["RayMatrix"];
  mVec    = mMatrix[[auxIdx]];
  n1      = Length[sdAug["NewExponents"]];

  Print["  Pivot table (cone ", sdAug["ConeIndex"], "):"];
  Print["  p | m_p | atilde (numeric) | admissible? | HasConstantTerm"];
  Do[
    If[mVec[[p]] != 0,
      Module[{row = computePivotRow[sdAug, liftData, p]},
        If[row =!= $Failed,
          Print["    p=", p,
                "  m_p=", row["mp"],
                "  atilde=", row["atilde"],
                "  admissible=", row["allPositive"],
                "  HasConst=", row["hasConst"]]
        ]
      ]
    ],
    {p, n1}
  ]
];
