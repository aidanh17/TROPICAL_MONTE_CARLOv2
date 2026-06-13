(* Probe: does lifting + VEGAS work on SMALL-coefficient integrands, and does
   the domain-indicator-cutoff bias shrink with VEGAS budget?  Compares
   lifted-MC vs lifted-VEGAS (several maxeval) vs NIntegrate reference. *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];
wd = FileNameJoin[{Directory[], "TEST", "INTERFILES"}];
If[!DirectoryQ[wd], CreateDirectory[wd]];

fmtR[x_, p_] := If[x == 0, "0", ToString[CForm[SetPrecision[x, p]]]];

probe[label_, spec_, liftRules_, ref_] := Module[{rMC, vbudgets, rVE},
  Print["===== ", label, "   (NIntegrate ref = ", fmtR[ref, 8], ") ====="];
  rMC = Quiet@EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> liftRules,
    "Integrator" -> "MC", "NSamples" -> 2000000,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wd];
  If[AssociationQ[rMC],
    Module[{r = rMC["Results"][[1]]},
      Print["  lifted-MC  (2e6):     ", fmtR[r["Re"], 7], " +/- ", fmtR[r["ReErr"], 2],
            "   relErr = ", fmtR[Abs[r["Re"]/ref - 1], 3]]],
    Print["  lifted-MC FAILED: ", rMC]];
  Do[
    rVE = Quiet@EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> liftRules,
      "Integrator" -> "VEGAS", "NSamples" -> vb, "VegasEpsRel" -> 1.*^-9,
      "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> wd];
    If[AssociationQ[rVE],
      Module[{r = rVE["Results"][[1]]},
        Print["  lifted-VEGAS (", ToString[vb], "): ", fmtR[r["Re"], 7], " +/- ", fmtR[r["ReErr"], 2],
              "   relErr = ", fmtR[Abs[r["Re"]/ref - 1], 3]]],
      Print["  lifted-VEGAS (", ToString[vb], ") FAILED: ", rVE]],
    {vb, {300000, 1000000, 3000000}}];
  Print[];
];

(* Case 1: 10^-4 coefficient (Example 20).  Auto-k=2 => z0 = 10^-2. *)
probe["1e-4 coeff: (1 + 1e-4 x1^2 + x2^2 + x1 x2^2)^-2",
  <|"Polynomials" -> {1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
  Quiet@NIntegrate[(1 + 10^-4 x1^2 + x2^2 + x1 x2^2)^(-2),
    {x1, 0, Infinity}, {x2, 0, Infinity}, MaxRecursion -> 40,
    PrecisionGoal -> 8, WorkingPrecision -> 30]];

(* Case 2: 10^-6 coefficient.  Auto-k=2 => z0 = 10^-3. *)
probe["1e-6 coeff: (1 + 1e-6 x1^2 + x2^2 + x1 x2^2)^-2",
  <|"Polynomials" -> {1 + 10^-6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
  Quiet@NIntegrate[(1 + 10^-6 x1^2 + x2^2 + x1 x2^2)^(-2),
    {x1, 0, Infinity}, {x2, 0, Infinity}, MaxRecursion -> 40,
    PrecisionGoal -> 8, WorkingPrecision -> 30]];

(* Case 3: LARGE coeff 10^6 (the L case) for contrast. *)
probe["1e6 coeff: (1 + 1e6 x1^2 + x2^2 + x1 x2^2)^-2",
  <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
  Quiet@NIntegrate[(1 + 10^6 x1^2 + x2^2 + x1 x2^2)^(-2),
    {x1, 0, Infinity}, {x2, 0, Infinity}, MaxRecursion -> 40,
    PrecisionGoal -> 8, WorkingPrecision -> 30]];
