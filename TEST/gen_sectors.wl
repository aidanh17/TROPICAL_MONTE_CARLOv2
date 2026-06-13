(* ============================================================================
   TEST/gen_sectors.wl

   Generates, for each test integrand, a self-contained C++ header that
   contains the EXACT flattened per-sector integrand functions produced by the
   tropical pipeline (GenerateCppMonteCarlo), with main() stripped off, plus
   the exact reference value baked in.  bench.cpp #includes one of these and
   integrates the same flattened sectors with MC / QMC / CUBA.

   The whole point: every integration method sees the *identical*
   pre-processed (tropically decomposed + flattened) integrand.  Differences
   are purely the sampler, not the integrand.

   Run:  wolframscript -file TEST/gen_sectors.wl
   Output: TEST/INTERFILES/sectors_<name>.hpp  +  TEST/INTERFILES/manifest.csv
   ============================================================================ *)

SetDirectory[DirectoryName[$InputFileName]];   (* TEST/ *)
pkgDir = ParentDirectory[];
Get[FileNameJoin[{pkgDir, "tropical_eval.wl"}]];

outDir = FileNameJoin[{Directory[], "INTERFILES"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

manifest = {};   (* {name, dim, nSectors, refRe, refIm, family} rows *)

(* The tropical fan is the NORMAL fan of the Newton polytope, which is
   scale-invariant.  The packaged fan code (translateToOriginInteger) needs an
   integer interior point, which thin lattice simplices like conv{0,e_i} lack
   for n>=4 -- it then leaks $Failed into the Polymake input.  Computing the
   fan from a scaled copy conv{0,K e_i} (which has an interior lattice point)
   yields the identical fan.  This was validated against direct NIntegrate. *)
computeFanRobust[verts_] := Module[{n, fd},
  n = Length[First[verts]];   (* ambient dimension *)
  (* Scale up front so conv{0, K e_i} has the interior lattice point (1,..,1):
     need sum < K, i.e. K > n.  K = n+2 works for simplices, scaled simplices,
     and the hypercube alike (scaling only enlarges, never removes, interior
     points), and avoids the slow failing Polymake roundtrips at K=1. *)
  Do[
    fd = Quiet@ComputeDecomposition[K*verts, "ShowProgress" -> False];
    If[ListQ[fd] && Length[fd] == 2 && FreeQ[fd, $Failed] &&
       Length[fd[[1]]] > 0,
      Print["  (fan from x", K, "-scaled polytope; normal fan is ",
            "scale-invariant)"];
      Return[fd, Module]],
    {K, {n + 2, 2 n + 4, 6 n + 6}}];
  $Failed
];

(* ---------------------------------------------------------------------------
   buildCase: process one numeric, parameter-free integrand spec into a
   stripped sector header.  exactValue is the analytic reference (complex ok).
   --------------------------------------------------------------------------- *)
buildCase[name_String, family_String, polys_List, monoExps_List,
          polyExps_List, vars_List, exactValue_] := Module[
  {spec, verts, fanData, dv, sl, allSD, conv, divg, tmpcpp, full, prefix,
   n, nSec, refRe, refIm, hdr, vr},

  n = Length[vars];
  Print["=============================================================="];
  Print["CASE ", name, "  (", family, ", n=", n, ")"];

  spec = <|
    "Polynomials"         -> polys,
    "MonomialExponents"   -> monoExps,
    "PolynomialExponents" -> polyExps,
    "Variables"           -> vars,
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None
  |>;

  (* Fan from the Newton polytope of the product of the polynomials.
     PolytopeVertices Expands, so product forms like Product[1+x_i] work. *)
  verts   = PolytopeVertices[Times @@ (polys^(-1)), vars];
  fanData = computeFanRobust[verts];
  If[fanData === $Failed,
    Print["  !! fan computation FAILED for ", name]; Return[$Failed]];
  {dv, sl} = fanData;
  Print["  fan: ", Length[dv], " rays, ", Length[sl], " simplicial cones"];

  (* Process every cone *)
  allSD = Table[
    ProcessSector[spec, dv, sl[[s]], s, "Verbose" -> False],
    {s, Length[sl]}];
  conv = Select[allSD, (AssociationQ[#] && !#["IsDivergent"]) &];
  divg = Select[allSD, (AssociationQ[#] &&  #["IsDivergent"]) &];
  Print["  sectors: ", Length[conv], " convergent, ", Length[divg], " divergent"];
  If[Length[divg] > 0,
    Print["  !! DIVERGENT SECTORS PRESENT -- spec rejected by v2 pipeline. ",
          "Adjust exponents."];
    Return[$Failed]];

  (* Independent sanity check against direct NIntegrate (small n only) *)
  If[n <= 4,
    vr = ValidateDecomposition[spec, fanData, {}, 4];
    If[AssociationQ[vr],
      Print["  ValidateDecomposition rel.err vs NIntegrate: ", vr["RelativeError"]]]
  ];

  (* Generate the full MC C++ (we only keep everything above main()) *)
  tmpcpp = FileNameJoin[{outDir, "tmp_" <> name <> ".cpp"}];
  GenerateCppMonteCarlo[conv, {}, spec, tmpcpp];
  full   = Import[tmpcpp, "Text"];
  prefix = First[StringSplit[full, "int main("]];   (* everything before main *)
  DeleteFile[tmpcpp];

  nSec  = Length[conv];
  refRe = N[Re[exactValue], 17];
  refIm = N[Im[exactValue], 17];

  hdr = prefix <> "\n" <>
    "// ---- benchmark metadata (appended by gen_sectors.wl) ----\n" <>
    "static const char* CASE_NAME = \"" <> name <> "\";\n" <>
    "static const char* CASE_FAMILY = \"" <> family <> "\";\n" <>
    "static const int   CASE_DIM = " <> ToString[n] <> ";\n" <>
    "static const int   CASE_NSEC = " <> ToString[nSec] <> ";\n" <>
    "static const double REF_RE = " <> ToString[CForm[refRe]] <> ";\n" <>
    "static const double REF_IM = " <> ToString[CForm[refIm]] <> ";\n";

  Export[FileNameJoin[{outDir, "sectors_" <> name <> ".hpp"}], hdr, "Text"];
  AppendTo[manifest, {name, n, nSec, refRe, refIm, family}];
  Print["  -> wrote sectors_", name, ".hpp   exact = ",
        refRe, If[refIm != 0, " + " <> ToString[refIm] <> " i", ""]];
];

(* ---------------------------------------------------------------------------
   Test integrand families  (all convergent, numeric, parameter-free,
   with analytic reference values)
   --------------------------------------------------------------------------- *)
X[k_] := x[k];
vlist[n_] := Table[x[k], {k, n}];

(* FAMILY A -- simplex / Dirichlet:  P = 1 + sum x_i,  A_i = 0,  B = -(n+2)
   exact = Gamma(2)/Gamma(n+2) = 1/(n+1)!                                     *)
doA[n_] := buildCase["A_simplex_n" <> ToString[n], "simplex",
  {1 + Total[vlist[n]]}, ConstantArray[0, n], {-(n + 2)}, vlist[n],
  1/(n + 1)!];

(* FAMILY Af -- fractional simplex: A_i = -1/2, B = -(n+1)
   exact = pi^(n/2) Gamma(n/2+1)/n!   (genuine x^{-1/2} singularities ->
   exercises the flattening transform)                                        *)
doAf[n_] := buildCase["Af_fracsimplex_n" <> ToString[n], "fracsimplex",
  {1 + Total[vlist[n]]}, ConstantArray[-1/2, n], {-(n + 1)}, vlist[n],
  Pi^(n/2) Gamma[n/2 + 1]/n!];

(* FAMILY B -- quadratic:  P = 1 + sum x_i^2,  A_i = 0,  B = -(n+1)
   exact = 2^-n pi^(n/2) Gamma(n/2+1)/n!   (richer monomial maps, factor-2 rays) *)
doB[n_] := buildCase["B_quad_n" <> ToString[n], "quadratic",
  {1 + Total[vlist[n]^2]}, ConstantArray[0, n], {-(n + 1)}, vlist[n],
  2^(-n) Pi^(n/2) Gamma[n/2 + 1]/n!];

(* FAMILY D -- product / many sectors:  P = prod (1+x_i), A_i = 0, B = -2
   Newton polytope = unit hypercube -> 2^n simplicial cones.
   exact = (1/(-B-1))^n = 1                                                    *)
doD[n_] := buildCase["D_product_n" <> ToString[n], "product",
  {Expand[Product[1 + x[k], {k, n}]]}, ConstantArray[0, n], {-2}, vlist[n],
  1];

(* ----- run the matrix ----- *)
Do[doA[n],  {n, {2, 3, 4, 5, 6, 7, 8}}];
Do[doAf[n], {n, {2, 4, 6}}];
Do[doB[n],  {n, {4, 6}}];
Do[doD[n],  {n, {3, 4, 5}}];

(* ----- write manifest ----- *)
Export[FileNameJoin[{outDir, "manifest.csv"}],
  Prepend[manifest, {"name", "dim", "nsectors", "ref_re", "ref_im", "family"}],
  "CSV"];

Print["\n==============================================================="];
Print["Done. ", Length[manifest], " cases written to ", outDir];
Print[Grid[Prepend[manifest, {"name","n","nSec","refRe","refIm","family"}],
  Frame -> All]];
