(* ============================================================================
   cuba_common.wl

   Shared helpers for CUBA cross-checks.  Loaded by
   tropical_eval_examples3.wl (Examples 21-23) and test_cuba.wl.

   generateCubaSource[spec, srcFile, opts] writes a self-contained C++
   program that integrates the ORIGINAL integrand of an IntegrandSpec
   (numeric coefficients/exponents, no kinematic symbols) over [0,inf)^n:

     - compactification x_i = t_i/(1-t_i), Jacobian prod_i (1-t_i)^-2
     - integrand exp( sum_i A_i log x_i + sum_j B_j log P_j(x) ) with
       std::complex arithmetic; two CUBA components (Re, Im)
     - t clamped to [1e-12, 1-1e-12] to avoid log(0) at the boundary
     - 1D integrands are padded with a dummy dimension (Cuhre requires
       ndim >= 2; the integrand simply ignores t_2)
     - cubacores(0,0): single process, deterministic, macOS-safe

   runCubaCheck[spec, tag, opts] generates, compiles (g++ -lcuba), runs,
   and parses the output into <|"VEGAS" -> <|...|>, "CUHRE" -> <|...|>|>
   with keys Value (complex), Error (complex), NEval, Fail.
   Returns $Failed (with a message) if CUBA is not installed.

   Options: "RunVegas" -> True|False (default True),
            "VegasEpsRel"  (default 5*10^-4),
            "CuhreEpsRel"  (default 10^-6),
            "MaxEval"      (default 2*10^6).

   CUBA: https://feynarts.de/cuba/   (macOS: `brew install cuba`)
   Auto-detected under /opt/homebrew, /usr/local, /usr.
   ============================================================================ *)

cubaFindPrefix[] := SelectFirst[
  {"/opt/homebrew", "/usr/local", "/usr"},
  FileExistsQ[FileNameJoin[{#, "include", "cuba.h"}]] &&
  (FileExistsQ[FileNameJoin[{#, "lib", "libcuba.a"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.dylib"}]] ||
   FileExistsQ[FileNameJoin[{#, "lib", "libcuba.so"}]]) &,
  $Failed];

cNum[r_] := ToString[CForm[N[r, 17]]];

cubaComplexStr[z_] := "cx(" <> cNum[Re[z]] <> ", " <> cNum[Im[z]] <> ")";

cubaPolyStr[poly_, vars_] := Module[{parsed},
  parsed = ParsePolynomial[poly, vars];
  StringRiffle[
    Table[
      StringRiffle[
        Join[
          {cNum[mono[[1]]]},
          MapIndexed[
            If[#1 == 0, Nothing,
              If[#1 == 1,
                "x[" <> ToString[#2[[1]] - 1] <> "]",
                "std::pow(x[" <> ToString[#2[[1]] - 1] <> "], " <>
                  ToString[#1] <> ")"]] &,
            mono[[2]]]
        ], " * "],
      {mono, parsed}],
    " + "]
];

generateCubaSource[spec_, srcFile_, opts___Rule] := Module[
  {polys, monoExps, polyExps, vars, n, ndim, polyDefs, logTerms, src,
   runVegas, vegasEpsRel, cuhreEpsRel, maxEval},

  {runVegas, vegasEpsRel, cuhreEpsRel, maxEval} =
    {"RunVegas", "VegasEpsRel", "CuhreEpsRel", "MaxEval"} /.
    {opts} /. {"RunVegas" -> True, "VegasEpsRel" -> 5*10^-4,
               "CuhreEpsRel" -> 10^-6, "MaxEval" -> 2000000};

  polys    = spec["Polynomials"];
  monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"];
  vars     = spec["Variables"];
  n        = Length[vars];
  ndim     = Max[n, 2];   (* Cuhre requires ndim >= 2: pad 1D integrands *)

  polyDefs = StringRiffle[
    Table["  const double P" <> ToString[j] <> " = " <>
          cubaPolyStr[polys[[j]], vars] <> ";",
      {j, Length[polys]}], "\n"];

  (* log-integrand: B_j log P_j + A_i log x_i (skip A_i == 0 terms so the
     boundary x_i -> 0 never produces 0 * log(0) = NaN) *)
  logTerms = StringRiffle[
    Join[
      Table[cubaComplexStr[polyExps[[j]]] <> " * std::log(P" <>
            ToString[j] <> ")", {j, Length[polys]}],
      Table[
        If[TrueQ[monoExps[[i]] == 0], Nothing,
          cubaComplexStr[monoExps[[i]]] <> " * std::log(x[" <>
            ToString[i - 1] <> "])"],
        {i, n}]
    ], "\n             + "];

  src = "// Auto-generated CUBA cross-check (direct integrand, no tropical
// decomposition).  Compactification x_i = t_i/(1-t_i).
#include <cmath>
#include <complex>
#include <cstdio>
extern \"C\" {
#include <cuba.h>
}
using cx = std::complex<double>;

static int Integrand(const int *ndim, const cubareal tt[],
                     const int *ncomp, cubareal ff[], void *userdata) {
  (void)ndim; (void)ncomp; (void)userdata;
  double x[" <> ToString[n] <> "];
  double jac = 1.0;
  for (int i = 0; i < " <> ToString[n] <> "; ++i) {
    double t = tt[i];
    if (t < 1e-12) t = 1e-12;
    if (t > 1.0 - 1e-12) t = 1.0 - 1e-12;
    const double u = 1.0 - t;
    x[i] = t / u;
    jac /= (u * u);
  }
" <> polyDefs <> "
  const cx logI = " <> logTerms <> ";
  const cx val = std::exp(logI) * jac;
  ff[0] = val.real();
  ff[1] = val.imag();
  return 0;
}

int main() {
  const int zero = 0;
  cubacores(&zero, &zero);   // single process: deterministic, macOS-safe

  const int ndim = " <> ToString[ndim] <> ", ncomp = 2;
  int neval, fail, nregions;
  cubareal integral[2], error[2], prob[2];
" <>
  If[TrueQ[runVegas],
"
  Vegas(ndim, ncomp, Integrand, nullptr, 1,
        " <> cNum[vegasEpsRel] <> ", 1e-12, 0, 1,
        0, " <> ToString[maxEval] <> ", 20000, 10000, 1000,
        0, nullptr, nullptr,
        &neval, &fail, integral, error, prob);
  std::printf(\"VEGAS %.12e %.12e %.3e %.3e %d %d\\n\",
              integral[0], integral[1], error[0], error[1], neval, fail);
", ""] <> "
  Cuhre(ndim, ncomp, Integrand, nullptr, 1,
        " <> cNum[cuhreEpsRel] <> ", 1e-12, 0, 0, " <> ToString[maxEval] <> ",
        0, nullptr, nullptr,
        &nregions, &neval, &fail, integral, error, prob);
  std::printf(\"CUHRE %.12e %.12e %.3e %.3e %d %d\\n\",
              integral[0], integral[1], error[0], error[1], neval, fail);

  return 0;
}
";
  Export[srcFile, src, "Text"];
  srcFile
];

runCubaCheck[spec_, tag_, opts___Rule] := Module[
  {prefix, srcFile, binFile, cc, run, lines, out = <||>},
  prefix = cubaFindPrefix[];
  If[prefix === $Failed,
    Print["  [CUBA not found -- skipping.  Install with `brew install cuba`",
          " or from https://feynarts.de/cuba/]"];
    Return[$Failed]];
  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  srcFile = FileNameJoin[{Directory[], "INTERFILES",
                          "cuba_check_" <> tag <> ".cpp"}];
  binFile = FileNameJoin[{Directory[], "INTERFILES", "cuba_check_" <> tag}];
  generateCubaSource[spec, srcFile, opts];
  cc = RunProcess[{"g++", "-O2", "-std=c++17",
    "-I" <> FileNameJoin[{prefix, "include"}],
    srcFile,
    "-L" <> FileNameJoin[{prefix, "lib"}],
    "-lcuba", "-lm", "-o", binFile}];
  If[cc["ExitCode"] != 0,
    Print["  [CUBA compile failed:]\n", cc["StandardError"]];
    Return[$Failed]];
  run = RunProcess[{binFile}];
  If[run["ExitCode"] != 0,
    Print["  [CUBA run failed:]\n", run["StandardError"]];
    Return[$Failed]];
  lines = StringSplit[run["StandardOutput"], "\n"];
  Do[
    Module[{fields = StringSplit[line]},
      If[Length[fields] >= 7,
        out[fields[[1]]] = <|
          "Value" -> (Read[StringToStream[fields[[2]]], Number] +
                      I Read[StringToStream[fields[[3]]], Number]),
          "Error" -> (Read[StringToStream[fields[[4]]], Number] +
                      I Read[StringToStream[fields[[5]]], Number]),
          "NEval" -> ToExpression[fields[[6]]],
          "Fail"  -> ToExpression[fields[[7]]]|>]],
    {line, lines}];
  out
];
