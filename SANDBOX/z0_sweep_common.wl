(* ============================================================================
   z0_sweep_common.wl  --  core engine for the z0 / k-sweep test harness.

   Implements the experiment of AUXT/z0_numerical_testing_plan.md under the rule
   of AUXT/anchor_selection_procedure.md, following the build guide AUXT/plan.md.

   Provides (all top-level, so z0_sweep_battery.wl can Get this file):
     kStarFor, kListFor          - anchor-rule prediction + sweep range
     tryAutoFan, degenerateQ     - lifted-fan construction + degeneracy detect
     runMCSeeds                  - multi-seed MC (compile once, re-seed cheaply)
     gateScan                    - per-cone gate metrics + sector data for trueSigma
     aggTrueSigma                - decisive true-sigma aggregation (plan 3.2)
     unliftedBaseline            - k-independent reference (sigma denominator)
     hardenedOracle              - truth value with trust flag (plan 5)
     sweepK                      - the core loop (plan 3.1)
     writeRows                   - CSV writer with the plan 3.4 schema

   Load: done automatically below (tropical_eval.wl + sandbox_lift_common.wl).
   ============================================================================ *)

Module[{pkgRoot},
  pkgRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
  SetDirectory[pkgRoot];
  Get[FileNameJoin[{pkgRoot, "tropical_eval.wl"}]];
  Get[FileNameJoin[{pkgRoot, "SANDBOX", "sandbox_lift_common.wl"}]];   (* trueSigma *)
];

(* ----------------------------------------------------------------------------
   Global configuration (the controls of plan 6).
   sigma convention everywhere: sigma = stderr * Sqrt[NSamples]  (bench:73,135)
   ---------------------------------------------------------------------------- *)
$Z0NSamples  = 1000000;            (* fixed N across all cells so VarRed compares *)
$Z0Threshold = 1000;               (* detector band tau *)
$Z0Seeds     = {42, 101, 202, 303, 404};   (* >= 5 distinct seeds (plan 6) *)
$Z0NThreads  = $ProcessorCount;
$Z0TruePG    = 4;                  (* PrecisionGoal for trueSigma NIntegrate *)


(* ----------------------------------------------------------------------------
   kStarFor / kListFor  --  the anchor rule and the sweep window.
   k* = max(1, ceil(|log10|C|| / log10 tau))   (anchor_selection_procedure 3)
   ---------------------------------------------------------------------------- *)
kStarFor[magnitude_, tau_:1000] :=
  Max[1, Ceiling[Abs[Log[10, N[magnitude]]] / Log[10, N[tau]]]];

(* default sweep 1 .. kStar+2 so the H3 plateau is visible *)
kListFor[kStar_Integer] := Range[1, kStar + 2];


(* ----------------------------------------------------------------------------
   tryAutoFan[liftedSpec] -> {rays, simplices}  or  $Failed
   Mirrors the EvaluateTropicalMCLifted auto-fan path + degeneracy guard
   (tropical_eval.wl:2340-2370): a lower-dimensional lifted Newton polytope
   yields simplices whose length != n+1, which we reject.
   ---------------------------------------------------------------------------- *)
tryAutoFan[liftedSpec_Association] := Module[
  {nv, verts, fan, dv, sl},
  nv = Length[liftedSpec["Variables"]];          (* = n+1 *)
  verts = Quiet[
    PolytopeVertices[(Times @@ liftedSpec["Polynomials"])^(-1),
                     liftedSpec["Variables"]],
    TropicalFan::polymake];
  If[!ListQ[verts], Return[$Failed]];
  fan = Quiet[ComputeDecomposition[verts, "ShowProgress" -> False],
              TropicalFan::polymake];
  If[!ListQ[fan] || Length[fan] < 2, Return[$Failed]];
  {dv, sl} = fan;
  If[Length[sl] == 0 || Length[dv] == 0 || !AllTrue[sl, Length[#] == nv &],
    Return[$Failed]];
  fan
];

degenerateQ[liftedSpec_Association] := (tryAutoFan[liftedSpec] === $Failed);


(* ----------------------------------------------------------------------------
   runMCSeeds[firstRes, nSamples, seeds, nThreads] -> list of per-seed assocs
   firstRes is the EvaluateTropicalMC[Lifted] result of seed seeds[[1]] (already
   compiled+run).  Subsequent seeds re-run the SAME compiled binary with the
   runtime seed override (argv[5]) -- no recompile.  This makes seed-robustness
   real and cheap (plan 6, "extend the driver").
   Each row: <|"seed","value","err","sigmaSample"|>.
   ---------------------------------------------------------------------------- *)
runMCSeeds[firstRes_, nSamples_Integer, seeds_List, nThreads_] := Module[
  {workDir, binary, kinFile, threadsStr, rows, sqrtN},
  If[!AssociationQ[firstRes] || !KeyExistsQ[firstRes, "Results"] ||
     !KeyExistsQ[firstRes, "ResultFile"],
    Return[$Failed]];
  sqrtN      = Sqrt[nSamples];
  workDir    = DirectoryName[firstRes["ResultFile"]];
  binary     = FileNameJoin[{workDir, "tropical_mc"}];
  kinFile    = FileNameJoin[{workDir, "kinematic_data.txt"}];
  threadsStr = ToString[nThreads];

  rows = {<|"seed" -> seeds[[1]],
            "value" -> firstRes["Results"][[1]]["Re"],
            "err"   -> firstRes["Results"][[1]]["ReErr"],
            "sigmaSample" -> firstRes["Results"][[1]]["ReErr"] * sqrtN|>};

  Do[
    Module[{rf, proc, raw, lines, parsed, val, err},
      rf = FileNameJoin[{workDir, "mc_results_seed_" <> ToString[sd] <> ".txt"}];
      proc = RunProcess[{binary, kinFile, rf,
        ToString[nSamples], threadsStr, ToString[sd]}];
      If[AssociationQ[proc] && proc["ExitCode"] == 0 && FileExistsQ[rf],
        raw   = Import[rf, "Text"];
        lines = Select[StringSplit[raw, "\n"], StringLength[StringTrim[#]] > 0 &];
        If[lines =!= {},
          parsed = Read[StringToStream[#], Number] & /@ StringSplit[lines[[1]]];
          val = parsed[[1]]; err = parsed[[3]];
          AppendTo[rows, <|"seed" -> sd, "value" -> val, "err" -> err,
                           "sigmaSample" -> err * sqrtN|>],
          AppendTo[rows, <|"seed" -> sd, "value" -> $Failed, "err" -> $Failed,
                           "sigmaSample" -> $Failed|>]
        ],
        AppendTo[rows, <|"seed" -> sd, "value" -> $Failed, "err" -> $Failed,
                         "sigmaSample" -> $Failed|>]
      ]
    ],
    {sd, Rest[seeds]}
  ];
  rows
];


(* ----------------------------------------------------------------------------
   gateScan[liftedSpec, fan, liftData] -> association of gate metrics
   Runs ProcessSectorLifted over every cone; collects per-cone EmptyDomain,
   HasConstantTerm, min Re atilde (NewExponents), |mp| (ZRow[[PivotIndex]]).
   Derives gateOK = exists surviving cone with HasConstantTerm and all Re atilde>0
   (surviving cones already have positive margins by construction in PSL).
   Returns also the surviving sector associations for trueSigma.
   ---------------------------------------------------------------------------- *)
gateScan[liftedSpec_Association, fan_List, liftData_Association] := Module[
  {dv, sl, sectors, survivors, hasConstList, margins, absMps, detMs,
   nDropped, nFailed, gateOK},
  {dv, sl} = fan;
  sectors = Table[
    Quiet @ ProcessSectorLifted[liftedSpec, dv, sl[[s]], s, liftData],
    {s, Length[sl]}
  ];

  (* surviving = AssociationQ, not EmptyDomain, not $Failed *)
  survivors = Select[sectors,
    AssociationQ[#] && !TrueQ[Lookup[#, "EmptyDomain", False]] &];
  nDropped = Count[sectors,
    a_ /; AssociationQ[a] && TrueQ[Lookup[a, "EmptyDomain", False]]];
  nFailed  = Count[sectors, $Failed];

  hasConstList = Lookup[#, "HasConstantTerm", Missing[]] & /@ survivors;
  margins = Table[
    Min[Re[N[Lookup[s, "NewExponents", {Infinity}]]]], {s, survivors}];
  absMps = Table[
    Module[{zr = Lookup[s, "ZRow", {}], pi = Lookup[s, "PivotIndex", 0]},
      If[ListQ[zr] && IntegerQ[pi] && 1 <= pi <= Length[zr],
        Abs[zr[[pi]]], Missing[]]],
    {s, survivors}];

  gateOK = MemberQ[hasConstList, True];
  detMs = Cases[Abs[N[Lookup[#, "DetM", Missing[]]]] & /@ survivors, _?NumericQ];

  <|"nSectors"     -> Length[sl],
    "nSurvivors"   -> Length[survivors],
    "nDropped"     -> nDropped,
    "nFailed"      -> nFailed,
    "hasConstList" -> hasConstList,
    "hasConstAny"  -> MemberQ[hasConstList, True],
    "gateOK"       -> gateOK,
    "minMargin"    -> If[margins === {}, Missing[], Min[margins]],
    "absMp"        -> With[{v = Cases[absMps, _?NumericQ]},
                        If[v === {}, Missing[], Min[v]]],
    "maxDetM"      -> If[detMs === {}, Missing[], Max[detMs]],
    "sumDetM"      -> If[detMs === {}, Missing[], Total[detMs]],
    "survivors"    -> survivors|>
];


(* ----------------------------------------------------------------------------
   geometryCost[caseId, tier, spec, prim, kList, mag, alpha, fanProvider]
   Geometry-only sweep (NO MC, NO trueSigma): per k, lift -> fan -> gateScan,
   recording cone count, #EmptyDomain, |det M|, gate status.  This is the
   Tier-5 cost curve (testing plan Tier 5) and is robust even where the integral
   is marginal.  Returns one row per k.
   ---------------------------------------------------------------------------- *)
geometryCost[caseId_, tier_, spec_, prim_, kList_List, mag_, alpha_,
             fanProvider_] := Module[{n, kStar},
  n = Length[spec["Variables"]];
  kStar = kStarFor[mag, $Z0Threshold];
  Table[
    Module[{rule, liftRes, liftedSpec, liftData, auxVar, z0, z0Num, fan, gate},
      rule = Append[prim, "k" -> k];
      liftRes = Quiet @ LiftCoefficients[spec, {rule}];
      Which[
        !AssociationQ[liftRes],
          withConfig[makeSkipRow[caseId, n, mag, alpha, k, kStar,
            "LiftCoefficients failed"], tier, "geom"],
        True,
          liftedSpec = liftRes["LiftedSpec"]; liftData = liftRes["LiftData"];
          auxVar = liftData["AuxVariable"]; z0 = liftData["z0"]; z0Num = N[z0];
          fan = fanProvider[k, liftedSpec];
          If[fan === $Failed,
            withConfig[makeSkipRow[caseId, n, mag, alpha, k, kStar,
              "no fan (degenerate)", z0Num], tier, "geom"],
            gate = gateScan[liftedSpec, fan, liftData];
            Join[
              makeSkipRow[caseId, n, mag, alpha, k, kStar,
                "geometry-only (no MC)", z0Num],
              <|"tier" -> tier, "config" -> "geom", "exact" -> True,
                "gateOK" -> gate["gateOK"], "hasConstAny" -> gate["hasConstAny"],
                "minMargin" -> gate["minMargin"], "absMp" -> gate["absMp"],
                "nDropped" -> gate["nDropped"], "nSectors" -> gate["nSectors"],
                "maxDetM" -> gate["maxDetM"], "sumDetM" -> gate["sumDetM"]|>]
          ]
      ]
    ],
    {k, kList}]
];


(* ----------------------------------------------------------------------------
   aggTrueSigma[survivors, kinRules, pg] -> rich association.
   Var(total) = Sum Var(sector) (sectors independent).  A single non-convergent
   I2 => infinite TOTAL true variance => heavy-tailed (the Case B/C trap; plan
   3.2).  We ALSO report the convergent-sector aggregate and the dominant-cone
   sigma/mu, because H3 is stated for "the dominant surviving cone" (testing
   plan H3) and phase1 reported Case A's honest reduction by excluding a tiny
   (0.2% mass) divergent sector.  The report applies the strict rule for the
   verdict and the convergent view for the H3 trend.
   Per-sector records carry I1, Sigma, I2Converged, HasConstantTerm.
   ---------------------------------------------------------------------------- *)
aggTrueSigma[survivors_List, kinRules_List:{}, pg_Integer:4] := Module[
  {per, conv, divg, sumI1All, sumI1Conv, sigStrict, sigConv, dom,
   domConv, domSOM},
  If[survivors === {},
    Return[<|"sigmaTrue" -> Indeterminate, "heavyTail" -> True,
             "sigmaTrueConv" -> Indeterminate, "convMassFrac" -> Indeterminate,
             "sumI1" -> Indeterminate, "nDivergent" -> 0, "nConv" -> 0,
             "domSigOverMu" -> Indeterminate, "domConverged" -> False,
             "perSector" -> {}|>]];

  per = Table[
    Module[{t = Quiet @ trueSigma[s, kinRules, pg]},
      If[AssociationQ[t],
        <|"I1" -> Lookup[t, "I1", 0], "Sigma" -> Lookup[t, "Sigma", 0],
          "I2Converged" -> TrueQ[t["I2Converged"]],
          "HasConstantTerm" -> Lookup[s, "HasConstantTerm", Missing[]]|>,
        Nothing]],
    {s, survivors}];

  conv = Select[per, #["I2Converged"] &];
  divg = Select[per, !#["I2Converged"] &];
  sumI1All  = Total[Abs[#["I1"]] & /@ per];
  sumI1Conv = Total[Abs[#["I1"]] & /@ conv];
  sigConv   = Sqrt[Total[(#["Sigma"] & /@ conv)^2]];
  sigStrict = If[divg === {}, sigConv, Infinity];

  dom = If[per === {}, None, First @ MaximalBy[per, Abs[#["I1"]] &]];
  domConv = dom =!= None && dom["I2Converged"];
  domSOM  = If[domConv && Abs[dom["I1"]] > 0,
    dom["Sigma"] / Abs[dom["I1"]], If[dom === None, Indeterminate, Infinity]];

  <|"sigmaTrue" -> sigStrict, "heavyTail" -> (divg =!= {}),
    "sigmaTrueConv" -> sigConv,
    "convMassFrac" -> If[sumI1All > 0, sumI1Conv / sumI1All, Indeterminate],
    "sumI1" -> Total[#["I1"] & /@ per],
    "nDivergent" -> Length[divg], "nConv" -> Length[conv],
    "domSigOverMu" -> domSOM, "domConverged" -> domConv,
    "perSector" -> per|>
];


(* ----------------------------------------------------------------------------
   unliftedSigmaTrue[spec, fan] -> aggTrueSigma over ProcessSector cones.
   Unlifted ProcessSector exposes FlattenedPolys/Prefactor/Dimension directly
   (DomainConstraint absent => None), exactly what trueSigma consumes.
   ---------------------------------------------------------------------------- *)
unliftedSigmaTrue[spec_Association, fan_List] := Module[
  {dv, sl, sectors},
  {dv, sl} = fan;
  sectors = Table[
    Module[{sd = Quiet @ ProcessSector[spec, dv, sl[[s]], s]},
      If[sd === $Failed || TrueQ[Lookup[sd, "IsDivergent", False]], Nothing, sd]],
    {s, Length[sl]}];
  aggTrueSigma[sectors, {}, $Z0TruePG]
];


(* ----------------------------------------------------------------------------
   hardenedOracle[spec, crossCheck] -> <|"value","trusted"|>
   Hardened direct NIntegrate of the original integrand, cross-checked against
   an independent estimate (the unlifted tropical sector I1-sum).  Accept only
   on < 1% agreement (plan 5 step 2-3); else trusted=False, rank by sigmaTrue.
   ---------------------------------------------------------------------------- *)
hardenedOracle[spec_Association, crossCheck_:Indeterminate] := Module[
  {vars, integrand, direct, rel, trusted},
  vars = spec["Variables"];
  integrand = (Times @@ MapThread[Power,
                {spec["Polynomials"], spec["PolynomialExponents"]}]) *
              (Times @@ MapThread[Power, {vars, spec["MonomialExponents"]}]);
  direct = Quiet @ NIntegrate[
    integrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 30, PrecisionGoal -> 6, Method -> "GlobalAdaptive"];

  trusted = NumericQ[direct];
  If[trusted && NumericQ[crossCheck] && Abs[crossCheck] > 0,
    rel = Abs[(direct - crossCheck) / crossCheck];
    trusted = (rel < 0.01)
  ];
  <|"value" -> If[NumericQ[direct], direct, Indeterminate],
    "trusted" -> TrueQ[trusted],
    "crossCheck" -> crossCheck|>
];


(* ----------------------------------------------------------------------------
   unliftedBaseline[caseSpec] -> association
   The k-independent reference shared by every VarRed.  Computes sigma_sample
   (multi-seed MC on the unit-coefficient proxy fan) AND sigma_true (decisive),
   plus the truth oracle for this case.
   caseSpec keys used: "Spec", "ProxyPoly", "CaseId".
   ---------------------------------------------------------------------------- *)
unliftedBaseline[caseSpec_Association] := Module[
  {spec, proxy, vars, verts, fan, res1, seedRows, ts, oracle, truth},
  spec  = caseSpec["Spec"];
  proxy = caseSpec["ProxyPoly"];
  vars  = spec["Variables"];

  verts = Quiet @ PolytopeVertices[proxy^(-1), vars];
  fan   = Quiet @ ComputeDecomposition[verts, "ShowProgress" -> False];
  If[!ListQ[fan] || Length[fan] < 2, Return[$Failed]];

  res1 = Quiet @ EvaluateTropicalMC[spec, fan, {{}},
    "NSamples" -> $Z0NSamples, "RunChecks" -> False, "Verbose" -> False,
    "SeedBase" -> $Z0Seeds[[1]]];
  If[!AssociationQ[res1] || !KeyExistsQ[res1, "Results"], Return[$Failed]];

  seedRows = runMCSeeds[res1, $Z0NSamples, $Z0Seeds, $Z0NThreads];
  ts       = unliftedSigmaTrue[spec, fan];

  (* truth: closed form if supplied, else hardened oracle cross-checked vs I1 sum *)
  oracle = hardenedOracle[spec, ts["sumI1"]];
  truth  = If[KeyExistsQ[caseSpec, "Truth"] && NumericQ[caseSpec["Truth"]],
    <|"value" -> caseSpec["Truth"], "trusted" -> True, "closedForm" -> True|>,
    Append[oracle, "closedForm" -> False]];

  <|"fan" -> fan, "nSectors" -> Length[fan[[2]]],
    "seedRows" -> seedRows,
    "sigmaSampleMean" -> Mean[Cases[Lookup[#, "sigmaSample"] & /@ seedRows, _?NumericQ]],
    "valueMean" -> Mean[Cases[Lookup[#, "value"] & /@ seedRows, _?NumericQ]],
    "sigmaTrue" -> ts["sigmaTrue"], "heavyTail" -> ts["heavyTail"],
    "sigmaTrueConv" -> ts["sigmaTrueConv"], "convMassFrac" -> ts["convMassFrac"],
    "domSigOverMu" -> ts["domSigOverMu"],
    "sumI1" -> ts["sumI1"],
    "truth" -> truth|>
];


(* ----------------------------------------------------------------------------
   evalRuleSet  --  evaluate ONE lift configuration (a list of rules) end to end:
   lift + exactness gate, fan, gate scan, true sigma, multi-seed MC, var-reduction
   and relErr.  Returns the per-seed rows.  Used by sweepK (single primary rule
   per k) and by the Tier-4 two-rule comparison.
   ---------------------------------------------------------------------------- *)
evalRuleSet[caseId_, tier_, config_, spec_, rules_List, kLabel_, kStar_,
            mag_, alpha_, fanProvider_, baseline_] := Module[
  {n, truth, sigTrueU, sigTrueUConv, rows,
   liftRes, liftedSpec, liftData, auxVar, z0, z0Num, exact,
   fan, gate, agg, seedRows, varRedTrue, varRedTrueConv, timing, baseRow},

  n = Length[spec["Variables"]];
  truth = baseline["truth"];
  sigTrueU = baseline["sigmaTrue"];
  sigTrueUConv = baseline["sigmaTrueConv"];
  rows = {};

  (* --- lift + exactness precondition (plan 3.1 step 2; plan 6) --- *)
  liftRes = Quiet @ LiftCoefficients[spec, rules];
  If[!AssociationQ[liftRes],
    Return[{withConfig[makeSkipRow[caseId, n, mag, alpha, kLabel, kStar,
      "LiftCoefficients failed"], tier, config]}]];
  liftedSpec = liftRes["LiftedSpec"];
  liftData   = liftRes["LiftData"];
  auxVar     = liftData["AuxVariable"];
  z0         = liftData["z0"];
  z0Num      = N[z0];

  exact = AllTrue[
    Range[Length[spec["Polynomials"]]],
    Quiet @ TrueQ[Simplify[
      (liftedSpec["Polynomials"][[#]] /. auxVar -> z0) -
      spec["Polynomials"][[#]]] === 0] &];
  If[!exact,
    Return[{withConfig[makeSkipRow[caseId, n, mag, alpha, kLabel, kStar,
      "liftidentity round-trip failed (non-exact)", z0Num], tier, config]}]];

  (* --- fan (auto, with degeneracy detect; or explicit) --- *)
  fan = fanProvider[kLabel, liftedSpec];
  If[fan === $Failed,
    Return[{withConfig[makeSkipRow[caseId, n, mag, alpha, kLabel, kStar,
      "no fan (degenerate, no explicit fan for this k)", z0Num], tier, config]}]];

  (* --- gate pass + true sigma (decisive) --- *)
  gate = gateScan[liftedSpec, fan, liftData];
  agg  = aggTrueSigma[gate["survivors"], {}, $Z0TruePG];

  (* --- MC over seeds --- *)
  {timing, seedRows} = AbsoluteTiming[
    Module[{res1},
      res1 = Quiet @ EvaluateTropicalMCLifted[spec, {{}},
        "LiftRules" -> rules, "FanData" -> fan,
        "NSamples" -> $Z0NSamples, "RunChecks" -> False, "Verbose" -> False,
        "SeedBase" -> $Z0Seeds[[1]]];
      If[AssociationQ[res1] && KeyExistsQ[res1, "Results"],
        runMCSeeds[res1, $Z0NSamples, $Z0Seeds, $Z0NThreads], $Failed]]];
  If[seedRows === $Failed,
    Return[{withConfig[makeSkipRow[caseId, n, mag, alpha, kLabel, kStar,
      "MC driver failed", z0Num], tier, config]}]];

  varRedTrue = Which[
    agg["sigmaTrue"] === Infinity, 0,
    !NumericQ[sigTrueU] || !NumericQ[agg["sigmaTrue"]] ||
      agg["sigmaTrue"] <= 0, Indeterminate,
    True, (sigTrueU / agg["sigmaTrue"])^2];
  varRedTrueConv = If[
    NumericQ[sigTrueUConv] && NumericQ[agg["sigmaTrueConv"]] &&
      agg["sigmaTrueConv"] > 0,
    (sigTrueUConv / agg["sigmaTrueConv"])^2, Indeterminate];

  baseRow = <|
    "caseId" -> caseId, "tier" -> tier, "config" -> config,
    "n" -> n, "magnitude" -> N[mag], "alpha" -> alpha,
    "k" -> kLabel, "kStar" -> kStar, "z0" -> z0, "z0Num" -> z0Num,
    "sigmaTrue" -> agg["sigmaTrue"], "heavyTail" -> agg["heavyTail"],
    "sigmaTrueConv" -> agg["sigmaTrueConv"],
    "convMassFrac" -> agg["convMassFrac"], "nDivergent" -> agg["nDivergent"],
    "domSigOverMu" -> agg["domSigOverMu"],
    "gateOK" -> gate["gateOK"], "hasConstAny" -> gate["hasConstAny"],
    "minMargin" -> gate["minMargin"], "absMp" -> gate["absMp"],
    "nDropped" -> gate["nDropped"], "nSectors" -> gate["nSectors"],
    "maxDetM" -> gate["maxDetM"], "sumDetM" -> gate["sumDetM"],
    "varRedTrue" -> varRedTrue, "varRedTrueConv" -> varRedTrueConv,
    "wallClock" -> timing, "exact" -> True, "skipReason" -> ""|>;

  Do[
    Module[{relErr},
      relErr = If[NumericQ[truth["value"]] && truth["trusted"] &&
                  NumericQ[sr["value"]] && Abs[truth["value"]] > 0,
        Abs[(sr["value"] - truth["value"]) / truth["value"]], Indeterminate];
      AppendTo[rows, Join[baseRow, <|
        "seed" -> sr["seed"], "value" -> sr["value"], "err" -> sr["err"],
        "sigmaSample" -> sr["sigmaSample"], "relErr" -> relErr|>]]],
    {sr, seedRows}];
  rows
];

(* tag a skip row with its tier + config (makeSkipRow leaves tier="") *)
withConfig[row_Association, tier_, config_] :=
  Join[row, <|"tier" -> tier, "config" -> config|>];


(* ----------------------------------------------------------------------------
   sweepK[caseSpec, baseline] -> list of row-assocs (one per (k, seed))
   The core loop (plan 3.1).  caseSpec keys:
     "CaseId","Spec","ProxyPoly","PrimaryRule" (no k),"Magnitude","Alpha",
     "FanProvider" (function (k, liftedSpec) -> fan|$Failed),
     optional "KList","Truth","Tier".
   baseline = unliftedBaseline[caseSpec] (shared denominator + truth).
   ---------------------------------------------------------------------------- *)
sweepK[caseSpec_Association, baseline_Association] := Module[
  {spec, prim, mag, alpha, kStar, kList, caseId, tier, fanProvider},
  spec   = caseSpec["Spec"];
  prim   = caseSpec["PrimaryRule"];
  mag    = caseSpec["Magnitude"];
  alpha  = caseSpec["Alpha"];
  caseId = caseSpec["CaseId"];
  tier   = Lookup[caseSpec, "Tier", ""];
  fanProvider = caseSpec["FanProvider"];
  kStar  = kStarFor[mag, $Z0Threshold];
  kList  = Lookup[caseSpec, "KList", kListFor[kStar]];

  Join @@ Table[
    evalRuleSet[caseId, tier, "single", spec, {Append[prim, "k" -> k]},
      k, kStar, mag, alpha, fanProvider, baseline],
    {k, kList}]
];


(* ----------------------------------------------------------------------------
   makeSkipRow  --  a fully-populated row for a skipped cell (plan 6: every row
   has exact=True or a skipReason).  Seed/MC/gate fields are NA.
   ---------------------------------------------------------------------------- *)
makeSkipRow[caseId_, n_, mag_, alpha_, k_, kStar_, reason_String,
            z0Num_:Indeterminate] := <|
  "caseId" -> caseId, "tier" -> "", "config" -> "single", "n" -> n,
  "magnitude" -> N[mag], "alpha" -> alpha, "k" -> k, "kStar" -> kStar,
  "z0" -> Indeterminate, "z0Num" -> z0Num, "seed" -> "",
  "value" -> Indeterminate, "err" -> Indeterminate,
  "sigmaSample" -> Indeterminate, "relErr" -> Indeterminate,
  "sigmaTrue" -> Indeterminate, "heavyTail" -> Indeterminate,
  "sigmaTrueConv" -> Indeterminate, "convMassFrac" -> Indeterminate,
  "nDivergent" -> Indeterminate, "domSigOverMu" -> Indeterminate,
  "gateOK" -> Indeterminate, "hasConstAny" -> Indeterminate,
  "minMargin" -> Indeterminate, "absMp" -> Indeterminate,
  "nDropped" -> Indeterminate, "nSectors" -> Indeterminate,
  "maxDetM" -> Indeterminate, "sumDetM" -> Indeterminate,
  "varRedTrue" -> Indeterminate, "varRedTrueConv" -> Indeterminate,
  "wallClock" -> Indeterminate, "exact" -> False, "skipReason" -> reason|>;


(* ----------------------------------------------------------------------------
   CSV schema + writer (plan 3.4).
   ---------------------------------------------------------------------------- *)
$Z0Columns = {
  "caseId", "tier", "config", "n", "magnitude", "alpha", "k", "kStar",
  "z0", "z0Num",
  "seed", "value", "err", "sigmaSample", "relErr",
  "sigmaTrue", "heavyTail", "sigmaTrueConv", "convMassFrac", "nDivergent",
  "domSigOverMu", "gateOK",
  "hasConstAny", "minMargin", "absMp", "nDropped", "nSectors",
  "maxDetM", "sumDetM",
  "varRedTrue", "varRedTrueConv", "wallClock", "exact", "skipReason"};

(* Flat, single-line, CSV-safe rendering.  ToString[N[x]] renders large/small
   reals in 2-D scientific form WITH EMBEDDED NEWLINES (1.x10^6), which breaks
   CSV rows; InputForm forces a flat, ReadList/ToExpression-parseable literal
   (1.*^6) instead -- it round-trips, unlike CForm's ambiguous "1.e6". *)
csvCell[x_] := Which[
  x === Infinity,          "inf",
  x === -Infinity,         "-inf",
  x === Indeterminate,     "NA",
  MatchQ[x, _Missing],     "NA",
  x === "",                "",
  x === True,              "True",
  x === False,             "False",
  StringQ[x],              x,
  IntegerQ[x],             ToString[x],
  NumericQ[x],             StringReplace[ToString[N[x], InputForm], "\n" -> " "],
  ListQ[x],                StringReplace[ToString[x, InputForm], "\n" -> " "],
  True,                    StringReplace[ToString[x], "\n" -> " "]
];

rowToList[row_Association] := csvCell[Lookup[row, #, ""]] & /@ $Z0Columns;

writeRows[rows_List, file_String] := Module[{data},
  data = Prepend[rowToList /@ rows, $Z0Columns];
  Export[file, data, "CSV"];
  Print["Wrote ", Length[rows], " rows to ", file];
  file
];

Print["z0_sweep_common.wl loaded: NSamples=", $Z0NSamples,
      "  seeds=", $Z0Seeds, "  tau=", $Z0Threshold];
