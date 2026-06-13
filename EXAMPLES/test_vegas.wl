(* ============================================================================
   test_vegas.wl — MC vs VEGAS parity demonstration

   Runs the SAME tropical specs through the pipeline BOTH ways:

       Integrator -> "MC"     the shipped plain Monte-Carlo sampler
       Integrator -> "VEGAS"  CUBA-Vegas adaptive importance sampling

   and asserts they return the SAME value (the decomposition produces the
   value; the sampler is a second-order optimization on top), with VEGAS
   typically 1-3 orders more accurate at a smaller evaluation budget.  Each
   case prints a side-by-side table and checks (cf. VEGAS_PLAN.md §5, R5):

     (a) |MC - VEGAS|  <=  K * sqrt(MCerr^2 + VEGASerr^2),  K = 5
         (the two independent estimates agree within combined error), and
     (b) MC and VEGAS each within a budget-appropriate relative tolerance of
         the analytic / NIntegrate reference.

   Cases:
     V1  A_simplex      P = 1 + sum x_i, A=0, B=-(n+2); exact 1/(n+1)!  (n=2,4)
     V2  Af_frac (n=4)  A_i=-1/2, B=-(n+1); exact pi^(n/2) Gamma(n/2+1)/n!
                        — the sqrt(.) cube-edge case where MC is weakest
     V3  complex B      P = 1 + lam x1^2 + x2^2 + x1 x2^2, B=-(2+I/2), scan lam;
                        reference via NIntegrate (Example 17). Exercises ncomp=2.
     V4  coeff batch    P = c0 + sum c_i x_i (n=4), B=-6; exact 1/(120 c0^2 prod c_i)
                        over a scan of random coefficient sets.
     L1  lifted, LARGE  auxiliary-variable lift + VEGAS on (1+1e6 x1^2+...)^-2;
                        MC edges VEGAS here (cutoff bias dominates for VEGAS).
     L2  lifted, small  auxiliary lift + VEGAS on (1+1e-4 x1^2+...)^-2; lifting +
                        VEGAS BEATS lifting + MC (heavy tail; ~0.3% vs ~1.2%).
     L3  lifted, small  auxiliary lift + VEGAS on (1+1e-6 x1^2+...)^-2; VEGAS wins
                        big (~3% vs plain MC ~13% with a lying error bar).

   Requires: tropical_fan.wl, tropical_eval.wl, Polymake, g++, and the CUBA
   library (`brew install cuba`).  If CUBA is absent the script reports that
   and exits cleanly.

   Run:  wolframscript -file EXAMPLES/test_vegas.wl
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

interfilesDir = FileNameJoin[{Directory[], "INTERFILES"}];
If[!DirectoryQ[interfilesDir], CreateDirectory[interfilesDir]];

(* CUBA presence check (same dirs as the package's findCubaPrefix). *)
cubaPrefix = SelectFirst[
  {"/opt/homebrew", "/usr/local", "/usr"},
  FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
  (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &, $Failed];
If[cubaPrefix === $Failed,
  Print["CUBA not found (looked under /opt/homebrew, /usr/local, /usr)."];
  Print["Install with `brew install cuba` or from https://feynarts.de/cuba/ ."];
  Print["test_vegas.wl requires CUBA for the Integrator -> \"VEGAS\" path; exiting."];
  Quit[]];

Print["tropical_eval.wl loaded; CUBA at ", cubaPrefix];
Print[];

$vegasRows = {};

(* ASCII-clean numeric formatter.  ScientificForm does not stringify in a text
   Print, and InputForm leaks precision backticks; CForm[SetPrecision[..]] gives
   clean C-style output (e.g. 0.382021, 2.69e-6).  Re/Im are formatted
   separately so complex values read as "a - b*I"; Chop drops the ~1e-18
   imaginary noise on real integrands. *)
fmtReal[x_, p_] := If[x == 0, "0", ToString[CForm[SetPrecision[x, p]]]];
fmt[v_, p_: 4] := With[{vv = Chop[N[v], 1.*^-14]},
  If[Im[vv] == 0,
    fmtReal[Re[vv], p],
    fmtReal[Re[vv], p] <> If[Im[vv] < 0, " - ", " + "] <>
      fmtReal[Abs[Im[vv]], p] <> "*I"]];

(* --------------------------------------------------------------------------
   computeFanRobust: simplex / linear-support fans occasionally need the dual
   vertices scaled by an integer before polymake returns a complete simplicial
   fan (cf. TEST/gen_kin.wl).  Try a few scalings.
   -------------------------------------------------------------------------- *)
computeFanRobust[verts_] := Module[{n = Length[First[verts]], fd},
  Do[
    fd = Quiet@ComputeDecomposition[K*verts, "ShowProgress" -> False];
    If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] &&
       Length[fd[[1]]] > 0, Return[fd, Module]],
    {K, {1, n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

(* --------------------------------------------------------------------------
   compareMCvsVEGAS — run a spec both ways, print a table, assert R5.

   refList: list of reference (analytic/NIntegrate) complex values, one per kp,
            or None.  Options tune budgets / tolerances / which gates are hard.
   -------------------------------------------------------------------------- *)

Options[compareMCvsVEGAS] = {
  "RefList"      -> None,
  "RefLabel"     -> "exact",
  "McSamples"    -> 1000000,
  "VegasMaxeval" -> 200000,
  "VegasEpsRel"  -> 1.*^-9,   (* small => VEGAS uses the full maxeval budget *)
  "McTol"        -> 2.*^-2,   (* rel tol: MC    vs reference *)
  "VegasTol"     -> 1.*^-2,   (* rel tol: VEGAS vs reference *)
  "SigmaK"       -> 5,        (* combined-error band for |MC - VEGAS| *)
  "GateMC"       -> True,     (* hard-gate MC-vs-ref (False => report only) *)
  "GateSigma"    -> True      (* hard-gate the combined band (False => report) *)
};

compareMCvsVEGAS[label_, spec_, fanData_, kinPoints_, OptionsPattern[]] :=
Module[{refList, refLabel, ns, ve, veEps, mcTol, veTol, sigmaK, gateMC, gateSig,
        rMC, rVE, nkp, casePass = True, notes = {}, anyKpFail = False},

  refList  = OptionValue["RefList"];
  refLabel = OptionValue["RefLabel"];
  ns       = OptionValue["McSamples"];
  ve       = OptionValue["VegasMaxeval"];
  veEps    = OptionValue["VegasEpsRel"];
  mcTol    = OptionValue["McTol"];
  veTol    = OptionValue["VegasTol"];
  sigmaK   = OptionValue["SigmaK"];
  gateMC   = OptionValue["GateMC"];
  gateSig  = OptionValue["GateSigma"];
  nkp      = Length[kinPoints];

  Print["=== ", label, "  (MC ", ns, " samples vs VEGAS ", ve,
        " maxeval/sector) ==="];

  rMC = EvaluateTropicalMC[spec, fanData, kinPoints,
    "Integrator" -> "MC", "NSamples" -> ns,
    "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> interfilesDir];
  rVE = EvaluateTropicalMC[spec, fanData, kinPoints,
    "Integrator" -> "VEGAS", "NSamples" -> ve, "VegasEpsRel" -> veEps,
    "RunChecks" -> False, "Verbose" -> False,
    "WorkingDirectory" -> interfilesDir];

  If[!AssociationQ[rMC] || !AssociationQ[rVE],
    Print["  FAIL: ", If[!AssociationQ[rMC], "MC", "VEGAS"],
          " run returned ", If[!AssociationQ[rMC], rMC, rVE]];
    AppendTo[$vegasRows, <|"Label" -> label, "Pass" -> False,
      "Notes" -> "run failed"|>];
    Return[False]];

  Do[
    Module[{rm, rv, mcVal, veVal, mcErr, veErr, diff, comb, sig,
            ref, mcDev, veDev, kpFail = False, tags = {}},
      rm = rMC["Results"][[i]];  rv = rVE["Results"][[i]];
      mcVal = rm["Re"] + I rm["Im"];  veVal = rv["Re"] + I rv["Im"];
      mcErr = Sqrt[rm["ReErr"]^2 + rm["ImErr"]^2];
      veErr = Sqrt[rv["ReErr"]^2 + rv["ImErr"]^2];
      diff  = Abs[mcVal - veVal];
      comb  = Max[Sqrt[mcErr^2 + veErr^2], 1.*^-300];
      sig   = diff / comb;

      Print["  kp ", i, If[nkp > 1 && kinPoints[[i]] =!= {},
                           " = " <> ToString[kinPoints[[i]]], ""], ":"];
      Print["    MC    = ", fmt[mcVal, 6], "  +/- ", fmt[mcErr, 3]];
      Print["    VEGAS = ", fmt[veVal, 6], "  +/- ", fmt[veErr, 3]];
      Print["    |MC-VEGAS| = ", fmt[diff, 3],
            "  = ", fmt[sig, 3], " sigma (gate ", sigmaK, ")"];

      If[gateSig && !TrueQ[N@sig <= sigmaK],
        kpFail = True; AppendTo[tags, "sigma>K"]];

      If[refList =!= None,
        ref   = refList[[i]];
        mcDev = Abs[(mcVal - ref)/ref];
        veDev = Abs[(veVal - ref)/ref];
        Print["    ", refLabel, " = ", fmt[ref, 6]];
        Print["    relErr: MC = ", fmt[mcDev, 3],
              " (gate ", mcTol, If[gateMC, "", ", report"], ")",
              "   VEGAS = ", fmt[veDev, 3],
              " (gate ", veTol, ")"];
        If[gateMC && !TrueQ[N@mcDev <= mcTol],
          kpFail = True; AppendTo[tags, "MC>tol"]];
        If[!TrueQ[N@veDev <= veTol],
          kpFail = True; AppendTo[tags, "VEGAS>tol"]];
      ];

      If[kpFail, anyKpFail = True;
        Print["    -> kp ", i, " FAIL: ", StringRiffle[tags, ", "]]];
    ],
    {i, nkp}];

  casePass = !anyKpFail;
  Print["  ", If[casePass, "PASS", "FAIL"], "  (", nkp, " kinematic point",
        If[nkp == 1, "", "s"], ")"];
  Print[];
  AppendTo[$vegasRows, <|"Label" -> label, "Pass" -> casePass,
    "NKP" -> nkp, "Notes" -> ""|>];
  casePass
];

(* ==========================================================================
   V1 — A_simplex, param-free.  P = 1 + sum x_i, A=0, B=-(n+2); exact 1/(n+1)!
   ========================================================================== *)

Do[
  Module[{vars, P, spec, verts, fan, exact},
    vars = Table[x[i], {i, n}];
    P    = 1 + Total[vars];
    spec = <|"Polynomials" -> {P}, "MonomialExponents" -> ConstantArray[0, n],
      "PolynomialExponents" -> {-(n + 2)}, "Variables" -> vars,
      "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
    verts = PolytopeVertices[P^(-1), vars];
    fan   = computeFanRobust[verts];
    exact = 1/(n + 1)!;
    compareMCvsVEGAS["V1 A_simplex n=" <> ToString[n], spec, fan, {{}},
      "RefList" -> {N[exact]}, "McSamples" -> 1000000, "VegasMaxeval" -> 200000];
  ],
  {n, {2, 4}}];

(* ==========================================================================
   V2 — Af_frac, n=4.  A_i=-1/2, B=-(n+1); exact pi^(n/2) Gamma(n/2+1)/n!
   The sqrt(.) cube-edge case: MC is weakest, VEGAS shines.
   ========================================================================== *)

Module[{n = 4, vars, P, spec, verts, fan, exact},
  vars = Table[x[i], {i, n}];
  P    = 1 + Total[vars];
  spec = <|"Polynomials" -> {P}, "MonomialExponents" -> ConstantArray[-1/2, n],
    "PolynomialExponents" -> {-(n + 1)}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[P^(-1), vars];
  fan   = computeFanRobust[verts];
  exact = Pi^(n/2) Gamma[n/2 + 1] / n!;
  (* MC is heavy on the sqrt edge: larger budget + looser MC tol, but keep the
     VEGAS gate strict (it nails the edge) and report the combined band. *)
  compareMCvsVEGAS["V2 Af_frac n=4 (sqrt edge)", spec, fan, {{}},
    "RefList" -> {N[exact]}, "McSamples" -> 3000000, "VegasMaxeval" -> 300000,
    "McTol" -> 3.*^-2, "VegasTol" -> 1.*^-2, "SigmaK" -> 6];
];

(* ==========================================================================
   V3 — complex exponent, kinematic scan (Example 17).
   P = 1 + lam x1^2 + x2^2 + x1 x2^2, B=-(2+I/2); reference via NIntegrate.
   ========================================================================== *)

Module[{lam, A, P, vars, spec, verts, fan, lamValues, kinPoints, refList},
  lam = Symbol["lam"];  A = 2 + I/2;
  P = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {P}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A}, "Variables" -> vars,
    "KinematicSymbols" -> {lam}, "RegulatorSymbol" -> None|>;
  (* fan is lam-independent: compute from the unit-coefficient proxy *)
  verts = PolytopeVertices[(P /. lam -> 1)^(-1), vars];
  fan   = computeFanRobust[verts];
  lamValues = {0.5, 1.0, 2.0, 4.0, 8.0};
  kinPoints = List /@ lamValues;
  refList = Table[
    Quiet@NIntegrate[
      (1 + lv t1^2 + t2^2 + t1 t2^2)^(-A),
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6],
    {lv, lamValues}];
  compareMCvsVEGAS["V3 complex B, lam scan", spec, fan, kinPoints,
    "RefList" -> refList, "RefLabel" -> "NIntegrate",
    "McSamples" -> 500000, "VegasMaxeval" -> 150000,
    "McTol" -> 2.*^-2, "VegasTol" -> 1.*^-2];
];

(* ==========================================================================
   V4 — coefficient batch scan.  P = c0 + sum c_i x_i (n=4), B=-6;
   exact = 1/(120 c0^2 prod_{i=1}^4 c_i).  Random coefficient sets.
   ========================================================================== *)

Module[{n = 4, vars, csyms, P, spec, verts, fan, nSets, kinPoints, refList},
  vars  = Table[x[i], {i, n}];
  csyms = Table[Symbol["c" <> ToString[i]], {i, 0, n}];
  P     = csyms[[1]] + Sum[csyms[[i + 1]] x[i], {i, n}];
  spec  = <|"Polynomials" -> {P}, "MonomialExponents" -> ConstantArray[0, n],
    "PolynomialExponents" -> {-6}, "Variables" -> vars,
    "KinematicSymbols" -> csyms, "RegulatorSymbol" -> None|>;
  (* fan is coefficient-independent: unit-coefficient proxy = simplex *)
  verts = PolytopeVertices[(1 + Total[vars])^(-1), vars];
  fan   = computeFanRobust[verts];
  nSets = 24;
  SeedRandom[12345];
  kinPoints = Table[
    Prepend[Table[Exp[RandomReal[{-0.7, 0.7}]], {n}], 1.0],  (* c0=1 *)
    {nSets}];
  refList = Table[
    Module[{c = kinPoints[[s]]},
      1.0/(120. * c[[1]]^2 * Times @@ c[[2 ;;]])],
    {s, nSets}];
  compareMCvsVEGAS["V4 coeff batch (n=4, " <> ToString[nSets] <> " sets)",
    spec, fan, kinPoints, "RefList" -> refList,
    "McSamples" -> 300000, "VegasMaxeval" -> 100000,
    "McTol" -> 2.*^-2, "VegasTol" -> 1.*^-2];
];

(* ==========================================================================
   LIFTED + VEGAS: the auxiliary-variable method (small / large coefficients)
   combined with VEGAS INSTEAD OF MC.  (VEGAS_PLAN.md §1.3, extended.)

   When a polynomial coefficient is extreme (very small or very large), the
   cleared polynomial Q develops a heavy tail that uniform MC samples poorly;
   the auxiliary-variable lift (Sections sec:lifting/sec:delta) tames it -- and
   the lift can be sampled with EITHER integrator.  The lifted sectors carry a
   HARD domain-indicator cutoff (a step discontinuity emitted inside
   integrand_conv_k), so VEGAS keeps a small residual bias that SHRINKS with
   budget; both samplers' error bars are optimistic in this heavy-tail regime.
   The regime story (measured below):
     - SMALL coefficients (1e-4, 1e-6): lifting + VEGAS BEATS lifting + MC,
       decisively at 1e-6 where plain MC is catastrophically off (~13%) with a
       lying error bar, while VEGAS lands ~3% and converging.
     - LARGE coefficient (1e6): MC edges VEGAS (the cutoff bias dominates for
       VEGAS while MC's tail is mild here).
   Each case runs EvaluateTropicalMCLifted BOTH ways at EQUAL budget, checks the
   VEGAS codegen was reached, that VEGAS is within tolerance of the NIntegrate
   reference, and that the regime's expected winner is the more accurate one.
   ========================================================================== *)

Options[compareLiftedMCvsVEGAS] = {
  "McSamples"      -> 2000000,
  "VegasMaxeval"   -> 2000000,
  "VegasEpsRel"    -> 1.*^-9,   (* small => VEGAS uses the full maxeval budget *)
  "VegasTol"       -> 5.*^-2,   (* loose: heavy-tail + cutoff regime *)
  "ExpectedWinner" -> "VEGAS"   (* "VEGAS" | "MC": who should be more accurate *)
};

(* NIntegrate has HoldAll: the integration variables must be passed as plain
   symbols (not vars[[i]] Part expressions, which it cannot localize). *)
niRef[poly_, v1_, v2_] := Module[{r},
  r = Quiet@NIntegrate[poly, Evaluate@{v1, 0, Infinity}, Evaluate@{v2, 0, Infinity},
        MaxRecursion -> 40, PrecisionGoal -> 8, WorkingPrecision -> 30];
  If[NumericQ[r], N[r], $Failed]];

compareLiftedMCvsVEGAS[label_, spec_, liftRules_, ref_, OptionsPattern[]] :=
Module[{ns, ve, veEps, veTol, winner, rMC, rVE, srcOK, mcVal, veVal, mcErr,
        veErr, mcDev, veDev, betterIsVE, pass = True, notes = {}},
  ns     = OptionValue["McSamples"];
  ve     = OptionValue["VegasMaxeval"];
  veEps  = OptionValue["VegasEpsRel"];
  veTol  = OptionValue["VegasTol"];
  winner = OptionValue["ExpectedWinner"];

  Print["=== ", label, "   (NIntegrate ref = ", fmt[ref, 7], ") ==="];
  Print["    lifted, MC ", ns, " samples vs VEGAS ", ve, " maxeval/sector"];

  rMC = Quiet@EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> liftRules,
    "Integrator" -> "MC", "NSamples" -> ns,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> interfilesDir];
  rVE = Quiet@EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> liftRules,
    "Integrator" -> "VEGAS", "NSamples" -> ve, "VegasEpsRel" -> veEps,
    "RunChecks" -> False, "Verbose" -> False, "WorkingDirectory" -> interfilesDir];

  srcOK = AssociationQ[rVE] && KeyExistsQ[rVE, "CppFile"] &&
          StringContainsQ[Import[rVE["CppFile"], "Text"], "Vegas("];

  If[AssociationQ[rMC] && AssociationQ[rVE],
    mcVal = rMC["Results"][[1]]["Re"];  veVal = rVE["Results"][[1]]["Re"];
    mcErr = rMC["Results"][[1]]["ReErr"];  veErr = rVE["Results"][[1]]["ReErr"];
    mcDev = Abs[mcVal/ref - 1];  veDev = Abs[veVal/ref - 1];
    betterIsVE = veDev < mcDev;
    Print["    lifted-MC    = ", fmt[mcVal, 7], " +/- ", fmt[mcErr, 2],
          "   relErr = ", fmt[mcDev, 3]];
    Print["    lifted-VEGAS = ", fmt[veVal, 7], " +/- ", fmt[veErr, 2],
          "   relErr = ", fmt[veDev, 3], "  (VEGAS codegen reached: ", srcOK, ")"];
    Print["    more accurate: ", If[betterIsVE, "VEGAS", "MC"],
          "   (expected ", winner, ")",
          "   [error bars optimistic in this regime; cutoff bias shrinks w/ budget]"];
    (* Hard gates *)
    If[!srcOK, pass = False; AppendTo[notes, "VEGAS codegen not reached"]];
    If[!AllTrue[{veVal, veErr}, NumericQ[#] && Abs[#] < 10^30 &],
      pass = False; AppendTo[notes, "VEGAS output non-finite"]];
    If[!TrueQ[veDev <= veTol], pass = False;
      AppendTo[notes, "VEGAS > " <> fmt[veTol, 2] <> " from reference"]];
    If[winner === "VEGAS" && !betterIsVE, pass = False;
      AppendTo[notes, "expected VEGAS more accurate than MC"]];
    If[winner === "MC" && betterIsVE, pass = False;
      AppendTo[notes, "expected MC more accurate than VEGAS"]];
    ,
    pass = False; AppendTo[notes, "lifted run failed"]];

  Print["  ", If[pass, "PASS", "FAIL"],
        If[notes =!= {}, "  (" <> StringRiffle[notes, "; "] <> ")", ""]];
  Print[];
  AppendTo[$vegasRows, <|"Label" -> label, "Pass" -> pass, "NKP" -> 1,
    "Notes" -> If[pass && AssociationQ[rVE], "winner: " <>
      If[TrueQ[veDev < mcDev], "VEGAS", "MC"], StringRiffle[notes, "; "]]|>];
  pass
];

(* L1 -- LARGE coeff 1e6 (k=2 => z0=1e2): MC edges VEGAS (cutoff bias dominates) *)
compareLiftedMCvsVEGAS["L1 lifted, LARGE coeff 1e6 (k=2)",
  <|"Polynomials" -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
  niRef[(1 + 10^6 x1^2 + x2^2 + x1 x2^2)^(-2), x1, x2],
  "ExpectedWinner" -> "MC", "VegasTol" -> 5.*^-2];

(* L2 -- small coeff 1e-4 (k=2 => z0=1e-2): lifting + VEGAS beats lifting + MC *)
compareLiftedMCvsVEGAS["L2 lifted, small coeff 1e-4 (k=2)",
  <|"Polynomials" -> {1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
  niRef[(1 + 10^-4 x1^2 + x2^2 + x1 x2^2)^(-2), x1, x2],
  "ExpectedWinner" -> "VEGAS", "VegasTol" -> 2.*^-2];

(* L3 -- small coeff 1e-6 (k=2 => z0=1e-3): VEGAS wins big; plain MC ~13% off *)
compareLiftedMCvsVEGAS["L3 lifted, small coeff 1e-6 (k=2)",
  <|"Polynomials" -> {1 + 10^-6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
    "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}, "k" -> 2|>},
  niRef[(1 + 10^-6 x1^2 + x2^2 + x1 x2^2)^(-2), x1, x2],
  "ExpectedWinner" -> "VEGAS", "VegasTol" -> 5.*^-2];

(* ==========================================================================
   Summary
   ========================================================================== *)

Module[{nPass, nFail},
  nPass = Count[$vegasRows, r_ /; TrueQ[r["Pass"]]];
  nFail = Length[$vegasRows] - nPass;
  Print["================================================================"];
  Print["  MC vs VEGAS parity summary (", Length[$vegasRows], " cases)"];
  Print["================================================================"];
  Do[
    Print["  ", StringPadRight[row["Label"], 36],
      If[TrueQ[row["Pass"]], "PASS", "FAIL"],
      If[row["Notes"] =!= "" && StringQ[row["Notes"]],
        "   [" <> row["Notes"] <> "]", ""]],
    {row, $vegasRows}];
  Print["----------------------------------------------------------------"];
  Print["  ", nPass, " PASSED, ", nFail, " FAILED"];
  Print["================================================================"];
];
