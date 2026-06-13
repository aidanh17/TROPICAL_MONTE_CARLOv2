(* Head-to-head: CUBA on the RAW (un-decomposed) integrand vs. CUBA on the
   tropically decomposed + flattened sectors.  Answers: is the tropical
   preprocessing still pulling its weight even with a good adaptive sampler?

   Uses the package's own raw-integrand CUBA generator (cuba_common.wl), which
   integrates the ORIGINAL integrand on [0,inf)^n via x_i = t_i/(1-t_i). *)

SetDirectory[DirectoryName[$InputFileName]];          (* TEST/ *)
$SkipPolymakeLoad = True;                              (* fans not needed here *)
Get[FileNameJoin[{ParentDirectory[], "tropical_eval.wl"}]];
Get[FileNameJoin[{ParentDirectory[], "EXAMPLES", "cuba_common.wl"}]];

mkspec[polys_, monoExps_, polyExps_, vars_] := <|
  "Polynomials" -> polys, "MonomialExponents" -> monoExps,
  "PolynomialExponents" -> polyExps, "Variables" -> vars,
  "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

vl[n_] := Table[x[k], {k, n}];

cases = {
  <|"tag" -> "A_simplex_n6", "spec" -> mkspec[{1 + Total[vl[6]]},
      ConstantArray[0, 6], {-8}, vl[6]], "exact" -> 1/7!|>,
  <|"tag" -> "Af_fracsimplex_n4", "spec" -> mkspec[{1 + Total[vl[4]]},
      ConstantArray[-1/2, 4], {-5}, vl[4]], "exact" -> Pi^2 Gamma[3]/4!|>,
  <|"tag" -> "B_quad_n6", "spec" -> mkspec[{1 + Total[vl[6]^2]},
      ConstantArray[0, 6], {-7}, vl[6]], "exact" -> 2^(-6) Pi^3 Gamma[4]/6!|>
};

maxEval = 2000000;   (* generous: more than the tropical runs used *)

Print["RAW-INTEGRAND CUBA (compactified x=t/(1-t)), MaxEval = ", maxEval];
Print[StringRepeat["-", 78]];
Do[
  Module[{tag, spec, exact, res, rel},
    tag = c["tag"]; spec = c["spec"]; exact = N[c["exact"], 16];
    Print["CASE ", tag, "   exact = ", exact];
    res = runCubaCheck[spec, "raw_" <> tag,
      "RunVegas" -> True, "VegasEpsRel" -> 10^-12,
      "CuhreEpsRel" -> 10^-12, "MaxEval" -> maxEval];
    If[res === $Failed, Print["  [skipped]"],
      Do[
        If[KeyExistsQ[res, m],
          rel = Abs[res[m]["Value"] - exact]/Abs[exact];
          Print["  ", m, ": value=", N[Re[res[m]["Value"]], 8],
                "  rel.err=", ScientificForm[N[rel], 3],
                "  neval=", res[m]["NEval"], "  fail=", res[m]["Fail"]]],
        {m, {"VEGAS", "CUHRE"}}]
    ];
    Print[""];
  ],
  {c, cases}];
