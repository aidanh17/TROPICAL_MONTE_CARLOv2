(* ============================================================================
   TEST/gen_kin.wl
   Kinematic-parametrized variant for the BATCHING test.

   Spec:  P = c0 + c1 x1 + c2 x2 + c3 x3 + c4 x4   (n=4),  A_i=0,  B=-6.
   KinematicSymbols = {c0,c1,c2,c3,c4}  -> the polynomial COEFFICIENTS are the
   batch parameters.  The Newton polytope (hence the fan and the flattened
   sector structure) is independent of the coefficient values, so ONE
   decomposition serves all 7200 coefficient sets -- exactly the property the
   batching exploits.

   Exact value per coefficient set:
       I = 1 / ( (n+1)! * c0^2 * prod_{i=1}^n c_i ).

   Emits TEST/INTERFILES/sectors_kin.hpp containing
     (a) the standard pipeline integrand  integrand_conv_N(y, params)  (params =
         the kinematic coefficients), main() stripped  -- the "naive" evaluator;
     (b) shared-basis monomial tables  KS_*  (basis exponents + per-param
         coefficient rows + prefactor + B), so the batch harness can compute the
         coefficient-INDEPENDENT monomial basis once per sample and reuse it for
         all 7200 coefficient sets.
   ============================================================================ *)

SetDirectory[DirectoryName[$InputFileName]];
pkgDir = ParentDirectory[];
Get[FileNameJoin[{pkgDir, "tropical_eval.wl"}]];
outDir = FileNameJoin[{Directory[], "INTERFILES"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

computeFanRobust[verts_] := Module[{n, fd},
  n = Length[First[verts]];
  Do[
    fd = Quiet@ComputeDecomposition[K*verts, "ShowProgress" -> False];
    If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] &&
       Length[fd[[1]]] > 0, Return[fd, Module]],
    {K, {n + 2, 2 n + 4, 6 n + 6}}];
  $Failed];

nEnv = Environment["KIN_N"];
n = If[StringQ[nEnv], ToExpression[nEnv], 4];   (* dimension; default 4 *)
vars = Table[x[k], {k, n}];
csyms = Table[Symbol["c" <> ToString[i]], {i, 0, n}];   (* params[0..n] *)
P = csyms[[1]] + Sum[csyms[[i + 1]] x[i], {i, n}];
Bexp = -(n + 2);

spec = <|
  "Polynomials" -> {P}, "MonomialExponents" -> ConstantArray[0, n],
  "PolynomialExponents" -> {Bexp}, "Variables" -> vars,
  "KinematicSymbols" -> csyms, "RegulatorSymbol" -> None|>;

Print["P = ", P, ",  B = ", Bexp, ",  n = ", n];
verts = PolytopeVertices[P^(-1), vars];
fan = computeFanRobust[verts];
{dv, sl} = fan;
Print["fan: ", Length[dv], " rays, ", Length[sl], " cones"];

allSD = Table[ProcessSector[spec, dv, sl[[s]], s, "Verbose" -> False],
  {s, Length[sl]}];
conv = Select[allSD, (AssociationQ[#] && !#["IsDivergent"]) &];
Print["convergent sectors: ", Length[conv], " / ", Length[sl]];
If[Length[conv] != Length[sl], Print["WARNING: some sectors divergent/failed"]];

(* (a) standard integrand header (strip main) *)
tmp = FileNameJoin[{outDir, "tmp_kin.cpp"}];
GenerateCppMonteCarlo[conv, {}, spec, tmp];
prefix = First[StringSplit[Import[tmp, "Text"], "int main("]];
DeleteFile[tmp];

(* (b) shared-basis monomial tables.
   Each sector here has exactly ONE polynomial (P) raised to B=-6.
   For each monomial of the flattened poly: {coeff (linear in csyms), alphas}. *)
nParams = Length[csyms];
cf[x_] := ToString[CForm[N[x, 17]]];

nmonoList = {}; prefacList = {}; alphaFlat = {}; cconstFlat = {}; crowFlat = {};
Do[
  Module[{sd, fp, pref, monos},
    sd = conv[[s]];
    pref = N[sd["Prefactor"], 17];     (* coefficient-independent number *)
    fp = sd["FlattenedPolys"][[1]];    (* single polynomial *)
    AppendTo[prefacList, pref];
    AppendTo[nmonoList, Length[fp]];
    Do[
      Module[{coeff, alphas, cconst, crow},
        coeff  = mono[[1]];
        alphas = mono[[2]];            (* length n, flattened exponents *)
        cconst = coeff /. Thread[csyms -> 0];
        crow   = Table[Coefficient[coeff, csyms[[p]]], {p, nParams}];
        alphaFlat  = Join[alphaFlat, N[alphas, 17]];
        AppendTo[cconstFlat, N[cconst, 17]];
        crowFlat   = Join[crowFlat, N[crow, 17]];
      ],
      {mono, fp}];
  ],
  {s, Length[conv]}];

(* offsets into the flattened monomial arrays *)
moff = Prepend[Accumulate[nmonoList], 0];

cArr[name_, vals_, type_: "double"] :=
  "static const " <> type <> " " <> name <> "[] = {" <>
  StringRiffle[
    If[type == "int", ToString /@ vals, cf /@ vals], ", "] <> "};\n";

tables = StringJoin[
  "\n// ---- shared-basis monomial tables (gen_kin.wl) ----\n",
  "static const int KS_NSEC = " <> ToString[Length[conv]] <> ";\n",
  "static const int KS_DIM = " <> ToString[n] <> ";\n",
  "static const int KS_NPARAMS = " <> ToString[nParams] <> ";\n",
  "static const double KS_B = " <> cf[Bexp] <> ";\n",
  cArr["KS_NMONO", nmonoList, "int"],
  cArr["KS_MOFF", moff, "int"],
  cArr["KS_PREFAC", prefacList, "double"],
  cArr["KS_ALPHA", alphaFlat, "double"],   (* totmono * KS_DIM, row-major *)
  cArr["KS_CCONST", cconstFlat, "double"], (* totmono *)
  cArr["KS_CROW", crowFlat, "double"],     (* totmono * KS_NPARAMS, row-major *)
  "static const int CASE_DIM = " <> ToString[n] <> ";\n",
  "static const int KS_NKINPARAMS = " <> ToString[nParams] <> ";\n"
];

outName = "sectors_kin_n" <> ToString[n] <> ".hpp";
Export[FileNameJoin[{outDir, outName}], prefix <> tables, "Text"];
Print["wrote ", outName, " : ", Length[conv], " sectors, ",
      Total[nmonoList], " total monomials, nParams=", nParams];
Print["exact(per kp) = 1/((n+1)! c0^2 prod c_i),  (n+1)! = ", (n + 1)!];
