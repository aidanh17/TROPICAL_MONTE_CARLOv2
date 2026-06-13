(* ============================================================================
   sandbox_toy0_1d.wl
   Hand-checkable 1D example (§3.6 of plan).
   I = Int_0^inf (1 + 10^6 x)^(-2) dx = 10^-6

   Tests k=1 and k=2 lifts.
   The 1D lifted polynomial has a 1D Newton polytope in 2D space;
   Polymake cannot produce proper 2D simplicial sectors for it.
   We therefore use manually-constructed complete unimodular triangulations
   of R^2 (the dual/normal fan space) that respect the normal-fan division.

   For k=1: polytope support {(0,0),(1,1)}, dividing hyperplane n1+n2=0.
   For k=2: polytope support {(0,0),(1,2)}, dividing hyperplane n1+2n2=0.

   The §3.6 hand-check: the sector M={{1,1},{0,1}} (rays {(-1,0),(-1,-1)}) is
   EmptyDomain for z0=10^6 > 1 — explicitly verified.

   Run: cd .../TROPICAL_MONTE_CARLOv2 && wolframscript -file SANDBOX/sandbox_toy0_1d.wl
   ============================================================================ *)

useSandboxImpl = False;   (* False = use package LiftCoefficients + ProcessSectorLifted *)

Module[{pkgRoot},
  pkgRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
  SetDirectory[pkgRoot];
  Get[FileNameJoin[{pkgRoot, "tropical_eval.wl"}]];
  Get[FileNameJoin[{pkgRoot, "SANDBOX", "sandbox_lift_common.wl"}]];
];

Print["========================================"];
Print["TOY0: 1D example  I = Int_0^inf (1+10^6 x)^{-2} dx = 10^-6"];
Print["========================================"];

exactAnswer = 10^-6;

spec0 = <|
  "Polynomials"         -> {1 + 10^6 * x[1]},
  "MonomialExponents"   -> {0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* Unlifted 1D sanity check via ValidateDecomposition *)
Print["\n--- Unlifted 1D sanity check ---"];
verts1d = PolytopeVertices[(1 + x[1])^(-1), {x[1]}];
fan1d   = ComputeDecomposition[verts1d, "ShowProgress" -> False];
val1d   = ValidateDecomposition[spec0, fan1d, {}, 5];
Print["Direct NIntegrate: ", val1d["DirectResult"]];
Print["Sector sum:        ", val1d["SectorSum"]];
Print["Relative error:    ", val1d["RelativeError"]];

(* ================================================================
   k = 1:  z0 = 10^6,  lifted P_L = 1 + z * x[1]
   Support {(0,0),(1,1)}.  Complete unimodular triangulation of R^2:
   7 sectors with ray pairs from rays {1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}.
   ================================================================ *)

Print["\n"];
Print["========================================"];
Print["k=1: z0 = 10^6"];
Print["========================================"];

If[useSandboxImpl,
  {liftedSpec1, liftData1} = liftSpec[spec0, 1, {1}, 1],
  Module[{lcRes = LiftCoefficients[spec0, {<|"PolyIndex"->1, "ExponentVector"->{1}, "k"->1|>}]},
    liftedSpec1 = lcRes["LiftedSpec"];
    liftData1   = lcRes["LiftData"]
  ]
];
Print["z0 = ", liftData1["z0"]];
Print["Lifted poly: ", liftedSpec1["Polynomials"]];

(* Identity check *)
Module[{lp = liftedSpec1["Polynomials"][[1]], z0v = liftData1["z0"],
        auxVar = liftData1["AuxVariable"]},
  If[Expand[lp /. auxVar -> z0v] === Expand[spec0["Polynomials"][[1]]],
    Print["Identity check (k=1): PASS"],
    Print["Identity check (k=1): FAIL  got=", Expand[lp /. auxVar->z0v]]
  ]
];

(* Complete unimodular triangulation of R^2 for support {(0,0),(1,1)}.
   Rays arranged cyclically around origin; consecutive pairs form sectors.
   All pairs have |det| = 1 (unimodular). *)
raysK1    = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
sectsK1   = Table[{i, If[i<Length[raysK1],i+1,1]}, {i,Length[raysK1]}];
(* Verify unimodular: *)
Assert[And @@ Table[Abs[Det[raysK1[[sectsK1[[s]]]]]] == 1,
                    {s, Length[sectsK1]}], "k=1 sectors not unimodular"];

Print["\n--- §3.6 hand-check: sector M={{1,1},{0,1}} must be EmptyDomain for z0=10^6 ---"];
(* This sector: selectedRays = {{-1,0},{-1,-1}} = rays 4,5 = sectK1[[4]] = {4,5} *)
Module[{sdChk, lsdChk},
  sdChk = ProcessSector[liftedSpec1, raysK1, {4,5}, -1];
  lsdChk = If[useSandboxImpl,
    deltaEliminate[sdChk, liftData1],
    ProcessSectorLifted[liftedSpec1, raysK1, {4,5}, -1, liftData1]
  ];
  Print["  M = ", sdChk["RayMatrix"]];
  Print["  z-row m = ", sdChk["RayMatrix"][[liftData1["AuxIndex"]]]];
  Print["  a = ", sdChk["NewExponents"]];
  If[KeyExistsQ[lsdChk,"EmptyDomain"] && TrueQ[lsdChk["EmptyDomain"]],
    Print["  §3.6 check: EmptyDomain = True — PASS (z0=", liftData1["z0"], " > 1, constant-root sector)"],
    Print["  §3.6 check: expected EmptyDomain=True but got: ", lsdChk]
  ]
];

Print["\n--- k=1 all sectors ---"];
total1 = 0;
Do[
  Module[{sdAug, lsd, sval},
    sdAug = ProcessSector[liftedSpec1, raysK1, sectsK1[[s]], s];
    lsd = If[useSandboxImpl,
      deltaEliminate[sdAug, liftData1],
      ProcessSectorLifted[liftedSpec1, raysK1, sectsK1[[s]], s, liftData1]
    ];
    Print["  Sector ", s, " rays=", raysK1[[sectsK1[[s]]]],
          "  M=", sdAug["RayMatrix"],
          "  z-row=", sdAug["RayMatrix"][[liftData1["AuxIndex"]]],
          "  a=", sdAug["NewExponents"]];
    If[KeyExistsQ[lsd,"EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["    -> EmptyDomain (dropped)"];
      ,
      If[lsd === $Failed,
        Print["    -> FAILED (deltaEliminate/$Failed)"],
        Print["    Pivot=", lsd["PivotIndex"],
              "  ATilde=", If[useSandboxImpl, lsd["ATilde"], lsd["NewExponents"]],
              "  Prefactor=", N[lsd["Prefactor"]]];
        Print["    ClearedPolys (Qtilde)=", lsd["ClearedPolys"]];
        Print["    DomainConstraint=", lsd["DomainConstraint"]];
        Print["    HasConstantTerm=", lsd["HasConstantTerm"]];
        sval = sectorNIntegrate[lsd, {}, 5];
        Print["    NIntegrate = ", sval];
        total1 += sval;
      ]
    ]
  ],
  {s, Length[sectsK1]}
];

Print["\nk=1 total sector sum = ", total1];
Print["Exact = ", exactAnswer];
Module[{relerr = Abs[(total1 - exactAnswer) / exactAnswer]},
  If[relerr < 10^-3,
    Print["TOY0 k=1 PASS relerr=", relerr],
    Print["TOY0 k=1 FAIL expected=", exactAnswer, " got=", total1, " relerr=", relerr]
  ]
];

(* ================================================================
   k = 2:  z0 = 10^3,  lifted P_L = 1 + z^2 * x[1]
   Support {(0,0),(1,2)}.  Complete unimodular triangulation of R^2
   with additional ray {-1,2} to respect the normal fan.
   ================================================================ *)

Print["\n"];
Print["========================================"];
Print["k=2: z0 = 10^3"];
Print["========================================"];

If[useSandboxImpl,
  {liftedSpec2, liftData2} = liftSpec[spec0, 1, {1}, 2],
  Module[{lcRes = LiftCoefficients[spec0, {<|"PolyIndex"->1, "ExponentVector"->{1}, "k"->2|>}]},
    liftedSpec2 = lcRes["LiftedSpec"];
    liftData2   = lcRes["LiftData"]
  ]
];
Print["z0 = ", liftData2["z0"]];
Print["Lifted poly: ", liftedSpec2["Polynomials"]];

Module[{lp = liftedSpec2["Polynomials"][[1]], z0v = liftData2["z0"],
        auxVar = liftData2["AuxVariable"]},
  If[Expand[lp /. auxVar -> z0v] === Expand[spec0["Polynomials"][[1]]],
    Print["Identity check (k=2): PASS"],
    Print["Identity check (k=2): FAIL  got=", Expand[lp /. auxVar->z0v]]
  ]
];

(* For k=2 support {(0,0),(1,2)}, divider ray (1,2)/(−1,−2).
   Unimodular triangulation adding ray (-1,2):
   rays: {1,0},{0,1},{-1,2},{-1,1},{-1,0},{-1,-1},{-1,-2},{0,-1},{1,-2},{1,-1}
   Consecutive unimodular pairs: *)
raysK2All = {{1,0},{0,1},{-1,2},{-1,1},{-1,0},{-1,-1},{-1,-2},{0,-1},{1,-2},{1,-1}};
(* Keep only consecutive pairs with |det|=1 *)
sectsK2 = {};
Do[
  Module[{r1 = raysK2All[[i]], r2 = raysK2All[[If[i < Length[raysK2All], i+1, 1]]]},
    If[Abs[Det[{r1,r2}]] == 1,
      AppendTo[sectsK2, {i, If[i < Length[raysK2All], i+1, 1]}]
    ]
  ],
  {i, Length[raysK2All]}
];
Print["k=2 sectors (unimodular): ", Length[sectsK2]];

Print["\n--- k=2 all sectors ---"];
total2 = 0;
Do[
  Module[{sdAug, lsd, sval, ridx = sectsK2[[s]]},
    sdAug = ProcessSector[liftedSpec2, raysK2All, ridx, s];
    lsd = If[useSandboxImpl,
      deltaEliminate[sdAug, liftData2],
      ProcessSectorLifted[liftedSpec2, raysK2All, ridx, s, liftData2]
    ];
    Print["  Sector ", s, " rays=", raysK2All[[ridx]],
          "  M=", sdAug["RayMatrix"],
          "  z-row=", sdAug["RayMatrix"][[liftData2["AuxIndex"]]],
          "  a=", sdAug["NewExponents"]];
    If[KeyExistsQ[lsd,"EmptyDomain"] && TrueQ[lsd["EmptyDomain"]],
      Print["    -> EmptyDomain (dropped)"];
      ,
      If[lsd === $Failed,
        Print["    -> FAILED (deltaEliminate/$Failed)"],
        Print["    Pivot=", lsd["PivotIndex"],
              "  ATilde=", If[useSandboxImpl, lsd["ATilde"], lsd["NewExponents"]],
              "  Prefactor=", N[lsd["Prefactor"]]];
        Print["    ClearedPolys (Qtilde)=", lsd["ClearedPolys"]];
        Print["    DomainConstraint=", lsd["DomainConstraint"]];
        Print["    HasConstantTerm=", lsd["HasConstantTerm"]];
        sval = sectorNIntegrate[lsd, {}, 5];
        Print["    NIntegrate = ", sval];
        total2 += sval;
      ]
    ]
  ],
  {s, Length[sectsK2]}
];

Print["\nk=2 total sector sum = ", total2];
Print["Exact = ", exactAnswer];
Module[{relerr = Abs[(total2 - exactAnswer) / exactAnswer]},
  If[relerr < 10^-3,
    Print["TOY0 k=2 PASS relerr=", relerr],
    Print["TOY0 k=2 FAIL expected=", exactAnswer, " got=", total2, " relerr=", relerr]
  ]
];

Print["\n========================================"];
Print["TOY0 COMPLETE"];
Print["========================================"];
