(* ============================================================================
   z0_sweep_battery.wl  --  the Tier 1-5 battery (driver script).

   Builds the magnitude-laddered case specs of z0_numerical_testing_plan.md 5,
   runs sweepK / evalRuleSet / geometryCost, and writes AUXT/z0_sweep_results.csv
   (one row per case,k,seed; plus geometry-cost rows for Tier 5).  The CSV is
   re-written after every case so a crash still leaves a valid partial file.

   Env controls:
     Z0_TIERS  = "1,2,3,4,5"   which tiers to run (default all)
     Z0_FAST   = "1"            quick shake-out: N=2e5, 2 seeds, short k-lists

   Run: cd .../TROPICAL_MONTE_CARLOv2 &&
        wolframscript -file SANDBOX/z0_sweep_battery.wl
   ============================================================================ *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "z0_sweep_common.wl"}]];

(* ---- mode toggles ---- *)
$fast = (Environment["Z0_FAST"] === "1");
If[$fast, $Z0NSamples = 200000; $Z0Seeds = {42, 101}];

$tiers = Module[{e = Environment["Z0_TIERS"]},
  If[StringQ[e], ToExpression /@ StringSplit[e, ","], {1, 2, 3, 4, 5}]];

$csvFile = FileNameJoin[{DirectoryName[$InputFileName], "..", "AUXT",
  "z0_sweep_results.csv"}];

Print["============================================================"];
Print["z0_sweep_battery: tiers=", $tiers, "  N=", $Z0NSamples,
      "  seeds=", $Z0Seeds, "  fast=", $fast];
Print["============================================================"];

$allRows = {};
flush[] := writeRows[$allRows, $csvFile];

(* ---- shared builders ---- *)
specT1[c_] := <|
  "Polynomials"         -> {1 + c x[1]^2 + x[2]^2 + x[1] x[2]^2},
  "MonomialExponents"   -> {0, 0}, "PolynomialExponents" -> {-2},
  "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

caseT1[id_, tier_, c_, kl_] := <|
  "CaseId" -> id, "Tier" -> tier, "Spec" -> specT1[c],
  "ProxyPoly" -> 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2,
  "PrimaryRule" -> <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}|>,
  "Magnitude" -> c, "Alpha" -> {2, 0},
  "FanProvider" -> Function[{k, ls}, tryAutoFan[ls]], "KList" -> kl|>;

runSweepCase[cs_] := Module[{base, rows},
  Print["\n### ", cs["CaseId"], " (", cs["Tier"], ")  |C|=", N[cs["Magnitude"]],
        "  k*=", kStarFor[cs["Magnitude"], $Z0Threshold]];
  base = unliftedBaseline[cs];
  If[base === $Failed,
    Print["  unlifted baseline FAILED -> recording skip rows"];
    rows = {withConfig[makeSkipRow[cs["CaseId"], Length[cs["Spec"]["Variables"]],
      cs["Magnitude"], cs["Alpha"], 0,
      kStarFor[cs["Magnitude"], $Z0Threshold], "unlifted baseline failed"],
      cs["Tier"], "single"]},
    Print["  unlifted: sigmaSample=", base["sigmaSampleMean"],
          " sigmaTrue=", base["sigmaTrue"], " sigmaTrueConv=", base["sigmaTrueConv"],
          " truth=", base["truth"]["value"], " trusted=", base["truth"]["trusted"]];
    rows = sweepK[cs, base]
  ];
  $allRows = Join[$allRows, rows]; flush[];
];


(* ============================================================================
   TIER 1  --  single extreme coefficient, non-degenerate (H1, H3)
   ============================================================================ *)
If[MemberQ[$tiers, 1],
  Print["\n========== TIER 1 =========="];
  Module[{cases},
    cases = {
      caseT1["T1a", "T1", 10^3, If[$fast, {1, 2}, Range[1, 4]]],
      caseT1["T1b", "T1", 10^4, If[$fast, {1, 2}, Range[1, 4]]],
      caseT1["T1c", "T1", 10^6, If[$fast, {2}, Range[1, 5]]],
      caseT1["T1d", "T1", 10^7, If[$fast, {2, 3}, Range[1, 5]]],
      caseT1["T1e", "T1", 10^9, If[$fast, {3}, Range[1, 6]]]};
    runSweepCase /@ cases;
  ]
];


(* ============================================================================
   TIER 2  --  small coefficient (symmetry of |log|C|| form)
   ============================================================================ *)
If[MemberQ[$tiers, 2],
  Print["\n========== TIER 2 =========="];
  Module[{cases},
    cases = {
      caseT1["T2a", "T2", 10^-4, If[$fast, {1, 2}, Range[1, 4]]],
      caseT1["T2b", "T2", 10^-6, If[$fast, {2}, Range[1, 4]]],
      caseT1["T2c", "T2", 10^-7, If[$fast, {3}, Range[1, 5]]]};
    runSweepCase /@ cases;
  ]
];


(* ============================================================================
   TIER 3  --  gate-dominated / degenerate (H2): Case C, explicit fans
   ============================================================================ *)
If[MemberQ[$tiers, 3],
  Print["\n========== TIER 3 =========="];
  Module[{caseC},
    (* lineality pair {k,0,-3}/{-k,0,3} on the {3,1} lift; 6 shared simplices
       (benchmark_results.md Case C). *)
    caseC = <|
      "CaseId" -> "T3_CaseC", "Tier" -> "T3",
      "Spec" -> <|"Polynomials" -> {1 + 10^8 x[1]^3 x[2] + x[2]^3},
        "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-3},
        "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {},
        "RegulatorSymbol" -> None|>,
      "ProxyPoly" -> 1 + x[1]^3 x[2] + x[2]^3,
      "PrimaryRule" -> <|"PolyIndex" -> 1, "ExponentVector" -> {3, 1}|>,
      "Magnitude" -> 10^8, "Alpha" -> {3, 1},
      "FanProvider" -> Function[{k, ls},
        {{{1, 0, 0}, {0, 1, 0}, {-1, -3, 0}, {k, 0, -3}, {-k, 0, 3}},
         {{1, 2, 4}, {2, 3, 4}, {3, 1, 4}, {1, 2, 5}, {2, 3, 5}, {3, 1, 5}}}],
      "KList" -> If[$fast, {2, 4}, {2, 3, 4}]|>;
    runSweepCase[caseC];
  ]
];


(* ============================================================================
   TIER 4  --  multi-rule / shared anchor (companion 5): Case B
   single-rule k* vs two-rule k=1, scored on sigmaTrue AND relErr.
   ============================================================================ *)
If[MemberQ[$tiers, 4],
  Print["\n========== TIER 4 =========="];
  Module[{specB, proxyB, csB, base, rowsSingle, rows2, kStarB},
    specB = <|"Polynomials" -> {10^-4 + 10^4 x[1]^2 + 10^-4 x[2]^2 +
        10^4 x[1] x[2]^2 + x[1]^2 x[2]},
      "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
      "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {},
      "RegulatorSymbol" -> None|>;
    proxyB = 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2 + x[1]^2 x[2];
    csB = <|"CaseId" -> "T4_CaseB", "Tier" -> "T4", "Spec" -> specB,
      "ProxyPoly" -> proxyB,
      "PrimaryRule" -> <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}|>,
      "Magnitude" -> 10^4, "Alpha" -> {2, 0},
      "FanProvider" -> Function[{k, ls}, tryAutoFan[ls]],
      "KList" -> If[$fast, {1, 2}, {1, 2, 3}]|>;
    kStarB = kStarFor[10^4, $Z0Threshold];

    Print["\n### T4_CaseB  |C_primary|=10^4  k*=", kStarB];
    base = unliftedBaseline[csB];
    If[base === $Failed,
      Print["  Case B unlifted baseline FAILED"];
      $allRows = Join[$allRows, {withConfig[makeSkipRow["T4_CaseB", 2, 10^4,
        {2, 0}, 0, kStarB, "unlifted baseline failed"], "T4", "single"]}]; flush[],

      Print["  unlifted: sigmaSample=", base["sigmaSampleMean"],
            " sigmaTrue=", base["sigmaTrue"], " truth=", base["truth"]["value"],
            " trusted=", base["truth"]["trusted"]];
      (* single-rule sweep on the primary {2,0} *)
      rowsSingle = sweepK[csB, base];
      (* two-rule k=1: primary {2,0} + secondary {1,2}, shared anchor *)
      rows2 = evalRuleSet["T4_CaseB", "T4", "2rule", specB,
        {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 1|>,
         <|"PolyIndex" -> 1, "ExponentVector" -> {1, 2}, "k" -> 1|>},
        1, kStarB, 10^4, {2, 0}, Function[{k, ls}, tryAutoFan[ls]], base];
      $allRows = Join[$allRows, rowsSingle, rows2]; flush[];
    ];
  ]
];


(* ============================================================================
   TIER 5  --  dimension scaling (geometry cost of larger k); geometry-only.
   ============================================================================ *)
If[MemberQ[$tiers, 5],
  Print["\n========== TIER 5 =========="];
  Module[{spec3, spec4, kl, rows},
    kl = If[$fast, {1, 2, 3}, Range[1, 5]];
    spec3 = <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3]},
      "MonomialExponents" -> {0, 0, 0}, "PolynomialExponents" -> {-2},
      "Variables" -> {x[1], x[2], x[3]}, "KinematicSymbols" -> {},
      "RegulatorSymbol" -> None|>;
    spec4 = <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 +
        x[1] x[2] x[3] x[4]},
      "MonomialExponents" -> {0, 0, 0, 0}, "PolynomialExponents" -> {-2},
      "Variables" -> {x[1], x[2], x[3], x[4]}, "KinematicSymbols" -> {},
      "RegulatorSymbol" -> None|>;

    Print["\n### T5_n3 geometry cost"];
    rows = geometryCost["T5_n3", "T5", spec3,
      <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0, 0}|>, kl, 10^6, {2, 0, 0},
      Function[{k, ls}, tryAutoFan[ls]]];
    $allRows = Join[$allRows, rows]; flush[];

    Print["\n### T5_n4 geometry cost"];
    rows = geometryCost["T5_n4", "T5", spec4,
      <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0, 0, 0}|>, kl, 10^6, {2, 0, 0, 0},
      Function[{k, ls}, tryAutoFan[ls]]];
    $allRows = Join[$allRows, rows]; flush[];
  ]
];

flush[];
Print["\n============================================================"];
Print["z0_sweep_battery COMPLETE -- ", Length[$allRows], " rows -> ", $csvFile];
Print["============================================================"];
