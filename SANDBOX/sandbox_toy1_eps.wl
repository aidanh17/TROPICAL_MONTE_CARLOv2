(* ============================================================================
   sandbox_toy1_eps.wl
   Small-coefficient toy with closed form (§5.3 of plan).
   I(eps) = Int_{[0,inf)^2} (1 + x[1] + eps^2 x[2])^{-3} dx1 dx2 = 1/(2 eps^2)
   eps = 10^-3 exact  =>  exact = 500000

   Fan construction note: the lifted polynomial 1+x[1]+x[3]^2*x[2] has
   only 3 monomials {(0,0,0),(1,0,0),(0,1,2)} forming a 2D triangle in 3D
   space — a degenerate (non-full-rank) Newton polytope. Polymake cannot
   produce proper 3-ray simplices from this. We use the EXTENDED FAN
   construction: take the original (unlifted) 2D fan and add a z-direction
   ray (0,0,-1) (for z0 < 1) or (0,0,+1) (for z0 > 1), extending each 2D
   sector to a 3D sector with z as an independent direction. This is
   mathematically equivalent to the product fan construction and is exact.

   Numbered procedure per §5.3.

   Run: cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file SANDBOX/sandbox_toy1_eps.wl
   ============================================================================ *)

useSandboxImpl = False;   (* False = use package LiftCoefficients + ProcessSectorLifted *)

(* ---- Step 1: Load packages ---- *)
Module[{pkgRoot},
  pkgRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
  SetDirectory[pkgRoot];
  Get[FileNameJoin[{pkgRoot, "tropical_eval.wl"}]];
  Get[FileNameJoin[{pkgRoot, "SANDBOX", "sandbox_lift_common.wl"}]];
];

Print["========================================"];
Print["TOY1: I(eps) = Int (1 + x1 + eps^2 x2)^{-3} dx1 dx2"];
Print["      eps=10^-3,  exact = 1/(2 eps^2) = 500000"];
Print["========================================"];

ep    = 10^-3;   (* exact rational *)
exact = 1 / (2 * ep^2);
Print["exact = ", exact];

(* ---- Step 2: Unlifted baseline ---- *)
Print["\n=== STEP 2: Unlifted baseline ==="];

spec0 = <|
  "Polynomials"         -> {1 + x[1] + 10^-6 * x[2]},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-3},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Fan from unit-coefficient proxy *)
vertsU = PolytopeVertices[(1 + x[1] + x[2])^(-1), {x[1], x[2]}];
fanU   = ComputeDecomposition[vertsU, "ShowProgress" -> False];
Print["Unlifted fan: ", Length[fanU[[2]]], " sectors"];
Print["Dual verts: ", fanU[[1]]];

(* ValidateDecomposition *)
valU = ValidateDecomposition[spec0, fanU, {}, 3];
Print["Unlifted direct NIntegrate = ", valU["DirectResult"]];
Print["Unlifted sector sum        = ", valU["SectorSum"]];
Print["Unlifted relative error    = ", valU["RelativeError"]];

(* unliftedSectorStats at 10^5 samples *)
Print["Running unliftedSectorStats at 10^5 samples..."];
ustats = unliftedSectorStats[spec0, fanU, 100000];

Print["Unlifted per-sector statistics:"];
Do[
  With[{st = ustats[[s]]},
    If[st["IsDivergent"],
      Print["  Sector ", s, ": DIVERGENT (skipped)"],
      Print["  Sector ", s, ": min|f|=", st["Min"],
            "  max|f|=", st["Max"],
            "  sigma=", st["Sigma"],
            "  mean=", st["Mean"]]
    ]
  ],
  {s, Length[ustats]}
];
uSigmaList = Map[#["Sigma"]&, Select[ustats, !TrueQ[#["IsDivergent"]]&]];
utSigma = N@Sqrt@Total[N[uSigmaList]^2];
Print["Unlifted total sigma = ", utSigma];

(* ================================================================
   Helper: buildPrincipalFan[k]
   Principled complete simplicial fan for lifted poly 1 + x1 + z^k x2.
   Support points in exponent space: {(0,0,0),(1,0,0),(0,1,k)}.
   The polytope lies in the plane k*e2 - e3 = 0; lineality direction nu=(0,k,-1).
   The 2D reduced normal fan (standard simplex) has rays:
     r1=(1,0), r2=(0,1), r3=(-1,-1)
   Lifted to R^3 with third coord 0; lineality split by +nu and -nu:
     dualVertices = {{1,0,0},{0,1,0},{-1,-1,0},{0,k,-1},{0,-k,1}}
     simplexList  = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}}
   These 6 cones cover R^3, are simplicial, and preserve tropical dominance.
   ================================================================ *)

buildPrincipalFan[k_Integer] :=
Module[{dv, sl},
  dv = {{1,0,0}, {0,1,0}, {-1,-1,0}, {0,k,-1}, {0,-k,1}};
  sl = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}};
  {dv, sl}
];

(* ================================================================
   Steps 3-7: k=2 lift
   ================================================================ *)

Print["\n=== STEP 3: Lift k=2 (z0 = 10^-3 exact) ==="];
If[useSandboxImpl,
  {liftedSpec2, liftData2} = liftSpec[spec0, 1, {0,1}, 2],
  Module[{lcRes = LiftCoefficients[spec0, {<|"PolyIndex"->1, "ExponentVector"->{0,1}, "k"->2|>}]},
    liftedSpec2 = lcRes["LiftedSpec"];
    liftData2   = lcRes["LiftData"]
  ]
];
Print["z0 = ", liftData2["z0"]];
Print["Lifted poly: ", liftedSpec2["Polynomials"]];
Print["Lifted variables: ", liftedSpec2["Variables"]];

(* Identity check *)
Module[{lp = liftedSpec2["Polynomials"][[1]], z0v = liftData2["z0"],
        auxVar = liftData2["AuxVariable"]},
  If[Expand[lp /. auxVar -> z0v] === Expand[spec0["Polynomials"][[1]]],
    Print["Identity check (k=2): PASS"],
    Print["Identity check (k=2): FAIL  got=", Expand[lp /. auxVar->z0v]]
  ]
];

Print["\n=== STEP 4: Lifted fan k=2 ==="];
{extVerts2, extSects2} = buildPrincipalFan[2];
Print["Lifted k=2 fan: ", Length[extSects2], " sectors  (unlifted had ", Length[fanU[[2]]], ")"];
Print["Principal dual verts (k=2): ", extVerts2];
Print["Principal simplices: ", extSects2];
Print["Note: principled normal fan — 6 sectors, simplicial, dominance-preserving"];

Print["\n=== STEP 5: Per-sector (k=2) ==="];
sectors2 = {};
Do[
  Module[{sdAug, lsd},
    sdAug = ProcessSector[liftedSpec2, extVerts2, extSects2[[s]], s];
    If[sdAug === $Failed, Print["  Sector ", s, ": ProcessSector FAILED"]; Return[]];
    lsd = If[useSandboxImpl,
      deltaEliminate[sdAug, liftData2],
      ProcessSectorLifted[liftedSpec2, extVerts2, extSects2[[s]], s, liftData2]
    ];
    Print["  Sector ", s, ":"];
    If[KeyExistsQ[lsd, "EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["    EmptyDomain -> True (dropped)"];
      ,
      If[lsd === $Failed,
        Print["    FAILED (deltaEliminate/$Failed)"],
        Print["    Pivot=", lsd["PivotIndex"],
              "  ZRow=", lsd["ZRow"],
              "  ATilde=", If[useSandboxImpl, lsd["ATilde"], lsd["NewExponents"]]];
        Print["    DomainConstraint=", lsd["DomainConstraint"]];
        Print["    HasConstantTerm=", lsd["HasConstantTerm"]];
        AppendTo[sectors2, lsd]
      ]
    ]
  ],
  {s, Length[extSects2]}
];
Print["Non-empty sectors k=2: ", Length[sectors2]];

Print["\n=== STEP 6: Exactness k=2 ==="];
Print["Running NIntegrate per sector (k=2, PG=3)..."];
svals2 = Table[sectorNIntegrate[sectors2[[i]], {}, 3], {i, Length[sectors2]}];
total2 = Total[svals2];
Print["k=2 sector NIntegrate results: ", svals2];
Print["k=2 total = ", total2];
Print["Exact     = ", exact];
Module[{relerr = Abs[(total2 - exact) / exact]},
  Print["k=2 relative error = ", relerr];
  If[relerr < 10^-3,
    Print["TOY1 k=2 PASS relerr=", N[relerr]],
    Print["TOY1 k=2 FAIL expected=", exact, " got=", total2, " relerr=", N[relerr]]
  ]
];
Print["Unlifted sector sum = ", valU["SectorSum"]];

Print["\n=== STEP 7: Variance comparison k=2 ==="];
Print["Running sectorMC at 10^5 per sector (k=2)..."];
mcResults2 = Table[sectorMC[sectors2[[i]], {}, 100000], {i, Length[sectors2]}];

Print["k=2 MC results per sector:"];
Do[
  With[{ff = mcResults2[[i]]["FeasibleFraction"],
        se = mcResults2[[i]]["StdErr"],
        mn = mcResults2[[i]]["MinMag"],
        mx = mcResults2[[i]]["MaxMag"]},
    Print["  Sector ", i, ": feasFrac=", ff,
          "  min|f|=", mn,
          "  max|f|=", mx,
          "  StdErr=", se]
  ],
  {i, Length[sectors2]}
];
lifted2SigmaList = Map[
  Function[{r}, r["SigmaPerSample"]],
  mcResults2
];
lt2Sigma = N@Sqrt@Total[N[lifted2SigmaList]^2];
Print["k=2 lifted total sigma   = ", lt2Sigma];
Print["Unlifted total sigma     = ", utSigma];
If[NumberQ[lt2Sigma] && lt2Sigma > 0,
  Print["Sigma ratio (unlifted/lifted) = ", utSigma / lt2Sigma]
];

(* ================================================================
   Step 8: Repeat with k=1
   ================================================================ *)

Print["\n=== STEP 8: k=1 (z0=10^-6) ==="];
If[useSandboxImpl,
  {liftedSpec1, liftData1} = liftSpec[spec0, 1, {0,1}, 1],
  Module[{lcRes = LiftCoefficients[spec0, {<|"PolyIndex"->1, "ExponentVector"->{0,1}, "k"->1|>}]},
    liftedSpec1 = lcRes["LiftedSpec"];
    liftData1   = lcRes["LiftData"]
  ]
];
Print["z0 = ", liftData1["z0"]];
Print["Lifted poly (k=1): ", liftedSpec1["Polynomials"]];

{extVerts1, extSects1} = buildPrincipalFan[1];
Print["k=1 fan: ", Length[extSects1], " sectors"];

sectors1 = {};
Do[
  Module[{sdAug, lsd},
    sdAug = ProcessSector[liftedSpec1, extVerts1, extSects1[[s]], s];
    If[sdAug === $Failed, Print["  Sector ", s, ": ProcessSector FAILED"]; Return[]];
    lsd = If[useSandboxImpl,
      deltaEliminate[sdAug, liftData1],
      ProcessSectorLifted[liftedSpec1, extVerts1, extSects1[[s]], s, liftData1]
    ];
    Print["  Sector ", s, ":"];
    If[KeyExistsQ[lsd, "EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["    EmptyDomain -> True (dropped)"];
      ,
      If[lsd === $Failed,
        Print["    FAILED (deltaEliminate/$Failed)"],
        Print["    Pivot=", lsd["PivotIndex"],
              "  ZRow=", lsd["ZRow"],
              "  ATilde=", If[useSandboxImpl, lsd["ATilde"], lsd["NewExponents"]]];
        Print["    DomainConstraint=", lsd["DomainConstraint"]];
        Print["    HasConstantTerm=", lsd["HasConstantTerm"]];
        AppendTo[sectors1, lsd]
      ]
    ]
  ],
  {s, Length[extSects1]}
];
Print["Non-empty sectors k=1: ", Length[sectors1]];

Print["Running NIntegrate per sector (k=1, PG=3)..."];
svals1 = Table[sectorNIntegrate[sectors1[[i]], {}, 3], {i, Length[sectors1]}];
total1 = Total[svals1];
Print["k=1 total = ", total1, "  exact = ", exact];
Module[{relerr = Abs[(total1 - exact) / exact]},
  If[relerr < 10^-3,
    Print["TOY1 k=1 PASS relerr=", N[relerr]],
    Print["TOY1 k=1 FAIL expected=", exact, " got=", total1, " relerr=", N[relerr]]
  ]
];

Print["Running sectorMC at 10^5 per sector (k=1)..."];
mcResults1 = Table[sectorMC[sectors1[[i]], {}, 100000], {i, Length[sectors1]}];
Do[
  With[{ff = mcResults1[[i]]["FeasibleFraction"],
        se = mcResults1[[i]]["StdErr"],
        mn = mcResults1[[i]]["MinMag"],
        mx = mcResults1[[i]]["MaxMag"]},
    Print["  Sector ", i, ": feasFrac=", ff,
          "  min|f|=", mn,
          "  max|f|=", mx,
          "  StdErr=", se]
  ],
  {i, Length[sectors1]}
];
lifted1SigmaList = Map[
  Function[{r}, r["SigmaPerSample"]],
  mcResults1
];
lt1Sigma = N@Sqrt@Total[N[lifted1SigmaList]^2];
Print["k=1 lifted total sigma = ", lt1Sigma];

(* ================================================================
   FINAL SUMMARY BLOCK
   ================================================================ *)
Print["\n"];
Print["========================================"];
Print["TOY1 FINAL SUMMARY BLOCK"];
Print["========================================"];
Print["Integral: I = Int_{[0,inf)^2} (1 + x1 + 10^-6 x2)^{-3} dx1 dx2"];
Print["Exact:    ", exact];
Print[""];
Print["--- Exactness ---"];
Print["Unlifted NIntegrate:  ", valU["DirectResult"],
      "  (sector sum: ", valU["SectorSum"], ")"];
Print["k=2 sector sum:  ", total2, "  (relerr vs exact: ",
      N[Abs[(total2-exact)/exact]], ")"];
Print["k=1 sector sum:  ", total1, "  (relerr vs exact: ",
      N[Abs[(total1-exact)/exact]], ")"];
Print[""];
Print["--- Sector counts ---"];
Print["Unlifted:         ", Length[fanU[[2]]], " sectors"];
Print["Lifted k=2:       ", Length[extSects2], " sectors total, ",
      Length[sectors2], " non-empty"];
Print["Lifted k=1:       ", Length[extSects1], " sectors total, ",
      Length[sectors1], " non-empty"];
Print[""];
Print["--- Sigma comparison ---"];
Print["Unlifted total sigma:  ", utSigma];
Print["Lifted k=2 sigma:      ", lt2Sigma];
Print["Lifted k=1 sigma:      ", lt1Sigma];
If[NumberQ[lt2Sigma] && lt2Sigma > 0,
  Print["Ratio unlifted/k=2:    ", utSigma / lt2Sigma]
];
If[NumberQ[lt1Sigma] && lt1Sigma > 0,
  Print["Ratio unlifted/k=1:    ", utSigma / lt1Sigma]
];
Print[""];
Print["--- Feasible fractions (k=2) ---"];
Do[
  Print["  Sector ", i, ": ", mcResults2[[i]]["FeasibleFraction"]],
  {i, Length[sectors2]}
];
Print["--- Feasible fractions (k=1) ---"];
Do[
  Print["  Sector ", i, ": ", mcResults1[[i]]["FeasibleFraction"]],
  {i, Length[sectors1]}
];
Print[""];
Print["--- HasConstantTerm ---"];
Print["k=2:"];
Do[
  Print["  Sector ", i, ": HasConstantTerm=", sectors2[[i]]["HasConstantTerm"]],
  {i, Length[sectors2]}
];
Print["k=1:"];
Do[
  Print["  Sector ", i, ": HasConstantTerm=", sectors1[[i]]["HasConstantTerm"]],
  {i, Length[sectors1]}
];
Print["========================================"];
Print["TOY1 COMPLETE"];
Print["========================================"];
