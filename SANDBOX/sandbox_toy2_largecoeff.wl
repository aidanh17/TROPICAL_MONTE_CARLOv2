(* ============================================================================
   sandbox_toy2_largecoeff.wl
   Test 6 Case A lifted (§5.4 of plan).
   P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2},  A={0,0}
   Lift rule: PolyIndex 1, ExponentVector {2,0}, k=2 => z0=10^3

   Run: cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file SANDBOX/sandbox_toy2_largecoeff.wl
   ============================================================================ *)

useSandboxImpl = False;   (* False = use package LiftCoefficients + ProcessSectorLifted *)

Module[{pkgRoot},
  pkgRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
  SetDirectory[pkgRoot];
  Get[FileNameJoin[{pkgRoot, "tropical_eval.wl"}]];
  Get[FileNameJoin[{pkgRoot, "SANDBOX", "sandbox_lift_common.wl"}]];
];

Print["========================================"];
Print["TOY2: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2}"];
Print["      Lift: 10^6 x1^2 -> z^2 x1^2,  z0=10^3"];
Print["========================================"];

(* ---- Original spec ---- *)
spec0 = <|
  "Polynomials"         -> {1 + 10^6 * x[1]^2 + x[2]^2 + x[1] * x[2]^2},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1], x[2]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* ---- Unlifted baseline ---- *)
Print["\n=== Unlifted baseline ==="];

vertsU = PolytopeVertices[
  (1 + x[1]^2 + x[2]^2 + x[1] * x[2]^2)^(-1),
  {x[1], x[2]}
];
fanU = ComputeDecomposition[vertsU, "ShowProgress" -> False];
Print["Unlifted fan: ", Length[fanU[[2]]], " sectors"];

(* ValidateDecomposition to get the reference value *)
Print["Running ValidateDecomposition (unlifted)..."];
valU = ValidateDecomposition[spec0, fanU, {}, 4];
Print["Unlifted direct NIntegrate = ", valU["DirectResult"]];
Print["Unlifted sector sum        = ", valU["SectorSum"]];
Print["Unlifted relative error    = ", valU["RelativeError"]];
refValue = valU["DirectResult"];

(* Unlifted stats *)
Print["Running unliftedSectorStats at 10^5 samples..."];
ustats = unliftedSectorStats[spec0, fanU, 100000];
Print["Unlifted per-sector stats:"];
Do[
  With[{st = ustats[[s]]},
    If[st["IsDivergent"],
      Print["  Sector ", s, ": DIVERGENT"],
      Print["  Sector ", s, ": min|f|=", st["Min"],
            "  max|f|=", st["Max"],
            "  sigma=", st["Sigma"]]
    ]
  ],
  {s, Length[ustats]}
];
uSigmaList2 = Map[#["Sigma"]&, Select[ustats, !TrueQ[#["IsDivergent"]]&]];
utSigma2 = N@Sqrt@Total[N[uSigmaList2]^2];
Print["Unlifted total sigma = ", utSigma2];

(* ---- Lift k=2 ---- *)
Print["\n=== Lift: 10^6 x1^2 -> z^2 x1^2,  k=2,  z0=10^3 ==="];

If[useSandboxImpl,
  {liftedSpec, liftData} = liftSpec[spec0, 1, {2, 0}, 2],
  Module[{lcRes = LiftCoefficients[spec0, {<|"PolyIndex"->1, "ExponentVector"->{2,0}, "k"->2|>}]},
    liftedSpec = lcRes["LiftedSpec"];
    liftData   = lcRes["LiftData"]
  ]
];
Print["z0 = ", liftData["z0"]];
Print["Lifted poly: ", liftedSpec["Polynomials"]];

(* Identity check *)
Module[{lp, z0v, auxVar},
  lp     = liftedSpec["Polynomials"][[1]];
  z0v    = liftData["z0"];
  auxVar = liftData["AuxVariable"];
  If[Expand[lp /. auxVar -> z0v] === Expand[spec0["Polynomials"][[1]]],
    Print["Identity check: PASS"],
    Print["Identity check: FAIL  got=", Expand[lp /. auxVar->z0v]]
  ]
];

(* Lifted fan: support of 1 + z^2 x1^2 + x2^2 + x1 x2^2 *)
vertsL = PolytopeVertices[
  (1 + x[3]^2 * x[1]^2 + x[2]^2 + x[1] * x[2]^2)^(-1),
  {x[1], x[2], x[3]}
];
fanL = ComputeDecomposition[vertsL, "ShowProgress" -> False];
Print["Lifted fan: ", Length[fanL[[2]]], " sectors  (unlifted had ", Length[fanU[[2]]], ")"];

(* Per-sector processing *)
Print["\n=== Per-sector delta elimination ==="];
liftedSectors = {};
emptyCount = 0;
noConstTermSectors = {};

Do[
  Module[{sdAug, lsd, atildeKey},
    sdAug = ProcessSector[liftedSpec, fanL[[1]], fanL[[2, s]], s];
    lsd = If[useSandboxImpl,
      deltaEliminate[sdAug, liftData],
      ProcessSectorLifted[liftedSpec, fanL[[1]], fanL[[2, s]], s, liftData]
    ];
    atildeKey = If[useSandboxImpl, "ATilde", "NewExponents"];
    If[KeyExistsQ[lsd, "EmptyDomain"] && lsd["EmptyDomain"],
      emptyCount++;
      Print["  Sector ", s, ": EmptyDomain (dropped)"];
      ,
      Print["  Sector ", s, ": Pivot=", lsd["PivotIndex"],
            "  ATilde=", N[lsd[atildeKey]],
            "  DomainConstraint=", If[lsd["DomainConstraint"] === None, "None", "present"],
            "  HasConstantTerm=", lsd["HasConstantTerm"]];
      AppendTo[liftedSectors, lsd];
      If[!lsd["HasConstantTerm"],
        AppendTo[noConstTermSectors, {s, lsd}]
      ]
    ]
  ],
  {s, Length[fanL[[2]]]}
];
Print["Non-empty sectors: ", Length[liftedSectors], "  Empty: ", emptyCount];

(* HasConstantTerm report *)
If[noConstTermSectors === {},
  Print["All non-empty sectors have HasConstantTerm -> True"],
  Print["WARNING: ", Length[noConstTermSectors],
        " sector(s) have HasConstantTerm -> False:"];
  Do[
    Print["  Sector ", noConstTermSectors[[i, 1]]],
    {i, Length[noConstTermSectors]}
  ]
];

(* ---- Exactness ---- *)
Print["\n=== Exactness ==="];
Print["Running NIntegrate per sector (PG=3)..."];
svals = Table[sectorNIntegrate[liftedSectors[[i]], {}, 3],
              {i, Length[liftedSectors]}];
totalLifted = Total[svals];
Print["Sector NIntegrate values: ", svals];
Print["Lifted total  = ", totalLifted];
Print["Direct NInt   = ", refValue];
Module[{relerr},
  relerr = Abs[(totalLifted - refValue) / refValue];
  Print["Relative error = ", N[relerr]];
  If[relerr < 10^-3,
    Print["TOY2 PASS relerr=", N[relerr]],
    Print["TOY2 FAIL expected=", refValue, " got=", totalLifted, " relerr=", N[relerr]]
  ]
];

(* ---- Variance comparison ---- *)
Print["\n=== Variance comparison ==="];
Print["Running sectorMC at 10^5 per sector..."];
mcResults = Table[sectorMC[liftedSectors[[i]], {}, 100000],
                  {i, Length[liftedSectors]}];

Print["Lifted per-sector MC:"];
Do[
  With[{ff  = mcResults[[i]]["FeasibleFraction"],
        se  = mcResults[[i]]["StdErr"],
        mn  = mcResults[[i]]["MinMag"],
        mx  = mcResults[[i]]["MaxMag"],
        hct = liftedSectors[[i]]["HasConstantTerm"]},
    Print["  Sector ", i, ": feasFrac=", ff,
          "  min|f|=", mn, "  max|f|=", mx,
          "  StdErr=", se, "  HasConstantTerm=", hct]
  ],
  {i, Length[liftedSectors]}
];
liftedSigmaList = Map[
  Function[{r}, r["SigmaPerSample"]],
  mcResults
];
ltSigma2 = N@Sqrt@Total[N[liftedSigmaList]^2];

(* Report on HasConstantTerm=False sectors vs others *)
If[noConstTermSectors =!= {},
  Print["\n--- HasConstantTerm=False sector magnitude analysis ---"];
  noConstIdx  = noConstTermSectors[[All, 1]];
  Print["Sectors with HasConstantTerm=False: ", noConstIdx];
  Do[
    Module[{sIdx = noConstTermSectors[[i, 1]],
            (* find position in liftedSectors *)
            lPos},
      lPos = Position[liftedSectors, noConstTermSectors[[i, 2]]][[1, 1]];
      Print["  Sector ", sIdx, " (pos ", lPos, "): ",
            "  min|f|=", mcResults[[lPos]]["MinMag"],
            "  max|f|=", mcResults[[lPos]]["MaxMag"],
            "  feasFrac=", mcResults[[lPos]]["FeasibleFraction"]]
    ],
    {i, Length[noConstTermSectors]}
  ];
  Print["Sectors with HasConstantTerm=True: magnitude spreads:"];
  Do[
    If[liftedSectors[[i]]["HasConstantTerm"],
      Print["  Sector ", i, ": min|f|=", mcResults[[i]]["MinMag"],
            "  max|f|=", mcResults[[i]]["MaxMag"]]
    ],
    {i, Length[liftedSectors]}
  ]
];

(* ================================================================
   FINAL SUMMARY BLOCK
   ================================================================ *)
Print["\n"];
Print["========================================"];
Print["TOY2 FINAL SUMMARY BLOCK"];
Print["========================================"];
Print["Integral: P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2,  B={-2}"];
Print["Reference (direct NIntegrate): ", refValue];
Print[""];
Print["--- Exactness ---"];
Print["Unlifted sector sum:  ", valU["SectorSum"],
      "  (relerr: ", valU["RelativeError"], ")"];
Print["Lifted sector sum:    ", totalLifted,
      "  (relerr vs NInt: ", N[Abs[(totalLifted-refValue)/refValue]], ")"];
Print[""];
Print["--- Sector counts ---"];
Print["Unlifted: ", Length[fanU[[2]]], " sectors"];
Print["Lifted:   ", Length[fanL[[2]]], " sectors total,  ",
      Length[liftedSectors], " non-empty,  ", emptyCount, " empty"];
Print[""];
Print["--- Sigma comparison ---"];
Print["Unlifted total sigma:  ", utSigma2];
Print["Lifted total sigma:    ", ltSigma2];
If[NumberQ[ltSigma2] && ltSigma2 > 0,
  Print["Ratio (unlifted/lifted): ", utSigma2 / ltSigma2]
];
Print[""];
Print["--- HasConstantTerm=False sectors ---"];
If[noConstTermSectors === {},
  Print["None"],
  Print["Count: ", Length[noConstTermSectors]];
  Do[
    Module[{sIdx = noConstTermSectors[[i, 1]], lPos},
      lPos = Position[liftedSectors, noConstTermSectors[[i, 2]]][[1, 1]];
      Print["  Sector ", sIdx, ": min|f|=", mcResults[[lPos]]["MinMag"],
            "  max|f|=", mcResults[[lPos]]["MaxMag"],
            "  feasFrac=", mcResults[[lPos]]["FeasibleFraction"]]
    ],
    {i, Length[noConstTermSectors]}
  ]
];
Print["--- Feasible fractions ---"];
Do[
  Print["  Sector ", i, ": ", mcResults[[i]]["FeasibleFraction"]],
  {i, Length[liftedSectors]}
];
Print["========================================"];
Print["TOY2 COMPLETE"];
Print["========================================"];
