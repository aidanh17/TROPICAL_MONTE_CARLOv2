(* ============================================================================
   sandbox_true_sigma.wl
   Phase 1 SANDBOX — decisive variance diagnostic via exact NIntegrate sigma.

   Computes true per-sample sigma for:
     Toy1 unlifted, Toy1 lifted k=2, Toy1 lifted k=1,
     Toy2 unlifted, Toy2 lifted.
   Uses trueSigma[sd, kinRules, pg] from sandbox_lift_common.wl.

   Run: cd .../TROPICAL_MONTE_CARLOv2 &&
        wolframscript -file SANDBOX/sandbox_true_sigma.wl
   ============================================================================ *)

Module[{pkgRoot},
  pkgRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
  SetDirectory[pkgRoot];
  Get[FileNameJoin[{pkgRoot, "tropical_eval.wl"}]];
  Get[FileNameJoin[{pkgRoot, "SANDBOX", "sandbox_lift_common.wl"}]];
];

Print["============================================================"];
Print["sandbox_true_sigma.wl  — TRUE SIGMA DIAGNOSTIC"];
Print["============================================================"];

(* ============================================================
   PART 0: Helper — print a labeled sector table and compute
           total sigma for one configuration.
   ============================================================ *)

(* printSectorTable[label, tsResults, hasConstList]
   tsResults : list of trueSigma outputs (one per non-empty sector)
   hasConstList : list of HasConstantTerm values (or {} for unlifted)
   Returns total sigma (or Infinity if any I2Converged->False). *)
printSectorTable[label_String, tsResults_List, hasConstList_List] :=
Module[
  {divergentIdx, sigmaList, totalSigma},

  Print["\n--- ", label, " per-sector results ---"];

  divergentIdx = {};

  Do[
    With[{r = tsResults[[i]]},
      If[hasConstList === {},
        (* unlifted: no HasConstantTerm *)
        Print["  Sector ", i,
              ": I1=", r["I1"],
              "  I2=", If[r["I2Converged"], r["I2"], "DIVERGENT"],
              "  Sigma=", If[r["I2Converged"], r["Sigma"], "INF"],
              "  I2Converged=", r["I2Converged"]],
        (* lifted *)
        Print["  Sector ", i,
              ": I1=", r["I1"],
              "  I2=", If[r["I2Converged"], r["I2"], "DIVERGENT"],
              "  Sigma=", If[r["I2Converged"], r["Sigma"], "INF"],
              "  I2Converged=", r["I2Converged"],
              "  HasConstantTerm=", hasConstList[[i]]]
      ];
      If[!r["I2Converged"], AppendTo[divergentIdx, i]]
    ],
    {i, Length[tsResults]}
  ];

  If[divergentIdx =!= {},
    Print["  TOTALSIGMA = INFINITE (divergent I2 in sectors ", divergentIdx, ")"];
    totalSigma = Infinity,
    sigmaList  = Map[#["Sigma"]&, tsResults];
    totalSigma = Sqrt[Total[sigmaList^2]];
    Print["  TOTALSIGMA = ", totalSigma]
  ];

  totalSigma
];


(* ============================================================
   TOY1 SETUP
   ============================================================ *)
Print["\n============================================================"];
Print["TOY1: P = 1 + x[1] + 10^-6 x[2],  B={-3},  exact = 500000"];
Print["============================================================"];

spec0Toy1 = <|
  "Polynomials"         -> {1 + x[1] + 10^-6 * x[2]},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-3},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

exactToy1 = 500000;

(* Toy1 unlifted fan — use unit-coefficient proxy *)
vertsU1 = PolytopeVertices[(1 + x[1] + x[2])^(-1), {x[1], x[2]}];
fanU1   = ComputeDecomposition[vertsU1, "ShowProgress" -> False];
Print["Toy1 unlifted fan: ", Length[fanU1[[2]]], " sectors"];

(* buildPrincipalFan[k]: same as in toy1 script *)
buildPrincipalFanT1[k_Integer] :=
Module[{dv, sl},
  dv = {{1,0,0}, {0,1,0}, {-1,-1,0}, {0,k,-1}, {0,-k,1}};
  sl = {{1,2,4},{2,3,4},{3,1,4},{1,2,5},{2,3,5},{3,1,5}};
  {dv, sl}
];

(* ============================================================
   TOY1 UNLIFTED
   ============================================================ *)
Print["\n=== TOY1 UNLIFTED ==="];

(* ProcessSector each cone -> trueSigma *)
toy1UResults = {};
Do[
  Module[{sd},
    sd = ProcessSector[spec0Toy1, fanU1[[1]], fanU1[[2, s]], s];
    If[sd === $Failed || TrueQ[sd["IsDivergent"]],
      Print["  Sector ", s, ": DIVERGENT or FAILED — skipping"];
      ,
      (* Unlifted ProcessSector returns FlattenedPolys etc. directly *)
      AppendTo[toy1UResults, trueSigma[sd, {}, 3]]
    ]
  ],
  {s, Length[fanU1[[2]]]}
];

toy1UI1Sum = Total[Map[#["I1"]&, toy1UResults]];
Print["Toy1 unlifted I1 sum = ", toy1UI1Sum,
      "  (sanity gate: within 0.1% of ", exactToy1, ")"];
If[Abs[(toy1UI1Sum - exactToy1)/exactToy1] < 0.001,
  Print["Toy1 unlifted I1 sanity: PASS"],
  Print["Toy1 unlifted I1 sanity: FAIL  got=", toy1UI1Sum]
];

toy1USigmaTotal = printSectorTable["TOY1 unlifted", toy1UResults, {}];

(* ============================================================
   TOY1 LIFTED k=2
   ============================================================ *)
Print["\n=== TOY1 LIFTED k=2 ==="];

{liftedSpecT1k2, liftDataT1k2} = liftSpec[spec0Toy1, 1, {0,1}, 2];
{extVertsT1k2, extSectsT1k2}   = buildPrincipalFanT1[2];
Print["z0 = ", liftDataT1k2["z0"]];

toy1K2Sectors = {};
toy1K2SdAugs  = {};  (* keep sdAug for pivot table *)
Do[
  Module[{sdAug, lsd},
    sdAug = ProcessSector[liftedSpecT1k2, extVertsT1k2, extSectsT1k2[[s]], s];
    If[sdAug === $Failed,
      Print["  Sector ", s, ": ProcessSector FAILED"]; Return[]
    ];
    lsd = deltaEliminate[sdAug, liftDataT1k2];
    If[KeyExistsQ[lsd, "EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["  Sector ", s, ": EmptyDomain (dropped)"],
      If[lsd === $Failed,
        Print["  Sector ", s, ": deltaEliminate FAILED"],
        AppendTo[toy1K2Sectors, lsd];
        AppendTo[toy1K2SdAugs, sdAug]
      ]
    ]
  ],
  {s, Length[extSectsT1k2]}
];
Print["Non-empty sectors k=2: ", Length[toy1K2Sectors]];

toy1K2TsResults = Table[trueSigma[toy1K2Sectors[[i]], {}, 3],
                         {i, Length[toy1K2Sectors]}];
toy1K2I1Sum = Total[Map[#["I1"]&, toy1K2TsResults]];
Print["Toy1 k=2 I1 sum = ", toy1K2I1Sum,
      "  (sanity gate: within 0.1% of ", exactToy1, ")"];
If[Abs[(toy1K2I1Sum - exactToy1)/exactToy1] < 0.001,
  Print["Toy1 k=2 I1 sanity: PASS"],
  Print["Toy1 k=2 I1 sanity: FAIL  got=", toy1K2I1Sum]
];

toy1K2HasConst = Map[#["HasConstantTerm"]&, toy1K2Sectors];
toy1K2SigmaTotal = printSectorTable["TOY1 lifted k=2", toy1K2TsResults, toy1K2HasConst];

(* ============================================================
   TOY1 LIFTED k=1
   ============================================================ *)
Print["\n=== TOY1 LIFTED k=1 ==="];

{liftedSpecT1k1, liftDataT1k1} = liftSpec[spec0Toy1, 1, {0,1}, 1];
{extVertsT1k1, extSectsT1k1}   = buildPrincipalFanT1[1];
Print["z0 = ", liftDataT1k1["z0"]];

toy1K1Sectors = {};
toy1K1SdAugs  = {};
Do[
  Module[{sdAug, lsd},
    sdAug = ProcessSector[liftedSpecT1k1, extVertsT1k1, extSectsT1k1[[s]], s];
    If[sdAug === $Failed,
      Print["  Sector ", s, ": ProcessSector FAILED"]; Return[]
    ];
    lsd = deltaEliminate[sdAug, liftDataT1k1];
    If[KeyExistsQ[lsd, "EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["  Sector ", s, ": EmptyDomain (dropped)"],
      If[lsd === $Failed,
        Print["  Sector ", s, ": deltaEliminate FAILED"],
        AppendTo[toy1K1Sectors, lsd];
        AppendTo[toy1K1SdAugs, sdAug]
      ]
    ]
  ],
  {s, Length[extSectsT1k1]}
];
Print["Non-empty sectors k=1: ", Length[toy1K1Sectors]];

toy1K1TsResults = Table[trueSigma[toy1K1Sectors[[i]], {}, 3],
                         {i, Length[toy1K1Sectors]}];
toy1K1I1Sum = Total[Map[#["I1"]&, toy1K1TsResults]];
Print["Toy1 k=1 I1 sum = ", toy1K1I1Sum,
      "  (sanity gate: within 0.1% of ", exactToy1, ")"];
If[Abs[(toy1K1I1Sum - exactToy1)/exactToy1] < 0.001,
  Print["Toy1 k=1 I1 sanity: PASS"],
  Print["Toy1 k=1 I1 sanity: FAIL  got=", toy1K1I1Sum]
];

toy1K1HasConst = Map[#["HasConstantTerm"]&, toy1K1Sectors];
toy1K1SigmaTotal = printSectorTable["TOY1 lifted k=1", toy1K1TsResults, toy1K1HasConst];


(* ============================================================
   TOY2 SETUP
   ============================================================ *)
Print["\n============================================================"];
Print["TOY2: P = 1 + 10^6 x[1]^2 + x[2]^2 + x[1]x[2]^2,  B={-2}"];
Print["============================================================"];

spec0Toy2 = <|
  "Polynomials"         -> {1 + 10^6 * x[1]^2 + x[2]^2 + x[1] * x[2]^2},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Exact reference: ValidateDecomposition gives ~0.000785 *)
(* We will compute I1 sum and check against ~0.000785 *)
exactToy2Approx = 0.000785;

(* ============================================================
   TOY2 UNLIFTED
   ============================================================ *)
Print["\n=== TOY2 UNLIFTED ==="];

(* Mirror sandbox_toy2_largecoeff.wl EXACTLY for the unlifted fan *)
vertsU2 = PolytopeVertices[
  (1 + x[1]^2 + x[2]^2 + x[1] * x[2]^2)^(-1),
  {x[1], x[2]}
];
fanU2 = ComputeDecomposition[vertsU2, "ShowProgress" -> False];
Print["Toy2 unlifted fan: ", Length[fanU2[[2]]], " sectors"];

toy2UResults = {};
Do[
  Module[{sd},
    sd = ProcessSector[spec0Toy2, fanU2[[1]], fanU2[[2, s]], s];
    If[sd === $Failed || TrueQ[sd["IsDivergent"]],
      Print["  Sector ", s, ": DIVERGENT or FAILED — skipping"],
      AppendTo[toy2UResults, trueSigma[sd, {}, 3]]
    ]
  ],
  {s, Length[fanU2[[2]]]}
];

toy2UI1Sum = Total[Map[#["I1"]&, toy2UResults]];
Print["Toy2 unlifted I1 sum = ", toy2UI1Sum,
      "  (sanity gate: within 0.1% of ~", exactToy2Approx, ")"];
toy2UGate = Abs[(toy2UI1Sum - exactToy2Approx)/exactToy2Approx] < 0.001;
If[toy2UGate,
  Print["Toy2 unlifted I1 sanity: PASS"],
  Print["Toy2 unlifted I1 sanity: FAIL  got=", toy2UI1Sum,
        "  (note: gate uses approx 0.000785 — checking abs diff < 1e-6)"];
  (* Fallback: absolute tolerance *)
  If[Abs[toy2UI1Sum - exactToy2Approx] < 1*^-6,
    Print["Toy2 unlifted I1 sanity (abs tol): PASS"]]
];

toy2USigmaTotal = printSectorTable["TOY2 unlifted", toy2UResults, {}];

(* ============================================================
   TOY2 LIFTED k=2
   ============================================================ *)
Print["\n=== TOY2 LIFTED k=2 ==="];

(* Mirror sandbox_toy2_largecoeff.wl EXACTLY *)
{liftedSpecT2, liftDataT2} = liftSpec[spec0Toy2, 1, {2, 0}, 2];
Print["z0 = ", liftDataT2["z0"]];

vertsL2 = PolytopeVertices[
  (1 + x[3]^2 * x[1]^2 + x[2]^2 + x[1] * x[2]^2)^(-1),
  {x[1], x[2], x[3]}
];
fanL2 = ComputeDecomposition[vertsL2, "ShowProgress" -> False];
Print["Toy2 lifted fan: ", Length[fanL2[[2]]], " sectors"];

toy2LiftedSectors = {};
toy2LiftedSdAugs  = {};
Do[
  Module[{sdAug, lsd},
    sdAug = ProcessSector[liftedSpecT2, fanL2[[1]], fanL2[[2, s]], s];
    lsd   = deltaEliminate[sdAug, liftDataT2];
    If[KeyExistsQ[lsd, "EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["  Sector ", s, ": EmptyDomain (dropped)"],
      If[lsd === $Failed,
        Print["  Sector ", s, ": deltaEliminate FAILED"],
        AppendTo[toy2LiftedSectors, lsd];
        AppendTo[toy2LiftedSdAugs, sdAug]
      ]
    ]
  ],
  {s, Length[fanL2[[2]]]}
];
Print["Non-empty sectors toy2 lifted: ", Length[toy2LiftedSectors]];

toy2LTsResults = Table[trueSigma[toy2LiftedSectors[[i]], {}, 3],
                        {i, Length[toy2LiftedSectors]}];
toy2LI1Sum = Total[Map[#["I1"]&, toy2LTsResults]];
Print["Toy2 lifted I1 sum = ", toy2LI1Sum,
      "  (sanity gate: within 0.1% of ~", exactToy2Approx, ")"];
toy2LGate = Abs[(toy2LI1Sum - exactToy2Approx)/exactToy2Approx] < 0.001;
If[toy2LGate,
  Print["Toy2 lifted I1 sanity: PASS"],
  Print["Toy2 lifted I1 sanity: FAIL  got=", toy2LI1Sum,
        "  (checking abs diff < 1e-6)"];
  If[Abs[toy2LI1Sum - exactToy2Approx] < 1*^-6,
    Print["Toy2 lifted I1 sanity (abs tol): PASS"]]
];

toy2LHasConst = Map[#["HasConstantTerm"]&, toy2LiftedSectors];
toy2LSigmaTotal = printSectorTable["TOY2 lifted", toy2LTsResults, toy2LHasConst];


(* ============================================================
   FINAL SUMMARY
   ============================================================ *)
Print["\n"];
Print["=== TRUE SIGMA SUMMARY ==="];

(* Helper to format sigma for summary *)
fmtSigma[s_] := If[s === Infinity || !NumberQ[s], "INF", ToString[N[s]]];

Print["TOY1 unlifted : ", fmtSigma[toy1USigmaTotal]];
Print["TOY1 lifted k=2: ", fmtSigma[toy1K2SigmaTotal]];
Print["TOY1 lifted k=1: ", fmtSigma[toy1K1SigmaTotal]];

If[NumberQ[toy1USigmaTotal] && NumberQ[toy1K2SigmaTotal] && toy1K2SigmaTotal > 0,
  Print["TOY1 ratio unlifted/lifted(k=2): ",
        fmtSigma[toy1USigmaTotal / toy1K2SigmaTotal]],
  Print["TOY1 ratio unlifted/lifted(k=2): N/A"]
];
Print["TOY2 unlifted : ", fmtSigma[toy2USigmaTotal]];
Print["TOY2 lifted : ", fmtSigma[toy2LSigmaTotal]];
If[NumberQ[toy2USigmaTotal] && NumberQ[toy2LSigmaTotal] && toy2LSigmaTotal > 0,
  Print["TOY2 ratio unlifted/lifted: ",
        fmtSigma[toy2USigmaTotal / toy2LSigmaTotal]],
  Print["TOY2 ratio unlifted/lifted: N/A"]
];

Print["\n=== I1 SANITY SUMMARY ==="];
Print["TOY1 unlifted I1 sum = ", toy1UI1Sum,
      "  (exact=", exactToy1, ")"];
Print["TOY1 k=2 I1 sum = ", toy1K2I1Sum,
      "  (exact=", exactToy1, ")"];
Print["TOY1 k=1 I1 sum = ", toy1K1I1Sum,
      "  (exact=", exactToy1, ")"];
Print["TOY2 unlifted I1 sum = ", toy2UI1Sum,
      "  (approx=", exactToy2Approx, ")"];
Print["TOY2 lifted I1 sum = ", toy2LI1Sum,
      "  (approx=", exactToy2Approx, ")"];

(* ============================================================
   PART B: Per-pivot tables for lifted sectors
   ============================================================ *)
Print["\n=== PART B: PER-PIVOT TABLES ==="];

Print["\n--- TOY1 lifted k=2 per-pivot tables ---"];
Do[
  Module[{lsd = toy1K2Sectors[[i]], sdAug = toy1K2SdAugs[[i]]},
    Print["Sector ", i, " (ConeIndex=", lsd["ConeIndex"],
          ",  HasConstantTerm=", lsd["HasConstantTerm"], "):"];
    pivotTable[sdAug, liftDataT1k2]
  ],
  {i, Length[toy1K2Sectors]}
];

Print["\n--- TOY1 lifted k=1 per-pivot tables ---"];
Do[
  Module[{lsd = toy1K1Sectors[[i]], sdAug = toy1K1SdAugs[[i]]},
    Print["Sector ", i, " (ConeIndex=", lsd["ConeIndex"],
          ",  HasConstantTerm=", lsd["HasConstantTerm"], "):"];
    pivotTable[sdAug, liftDataT1k1]
  ],
  {i, Length[toy1K1Sectors]}
];

Print["\n--- TOY2 lifted per-pivot tables (HasConstantTerm=False sectors) ---"];
Module[{noConstIdx},
  noConstIdx = Select[Range[Length[toy2LiftedSectors]],
                      !toy2LHasConst[[#]]&];
  If[noConstIdx === {},
    Print["All toy2 lifted sectors have HasConstantTerm -> True — no table needed."],
    Do[
      Module[{i = idx, lsd = toy2LiftedSectors[[idx]],
              sdAug = toy2LiftedSdAugs[[idx]]},
        Print["Sector ", i, " (ConeIndex=", lsd["ConeIndex"],
              ",  HasConstantTerm=", lsd["HasConstantTerm"], "):"];
        pivotTable[sdAug, liftDataT2]
      ],
      {idx, noConstIdx}
    ]
  ]
];

Print["\n============================================================"];
Print["sandbox_true_sigma.wl COMPLETE"];
Print["============================================================"];
