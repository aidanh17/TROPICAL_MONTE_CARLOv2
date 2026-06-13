(* Faithful with/without-lifting MC comparison for Toy 0:
     I = Int_0^inf (1 + 10^6 x)^-2 dx = 10^-6.
   Uses the package's own ProcessSector / ProcessSectorLifted so the
   flattened integrands are exactly what the C++ MC would sample.
   Replicates the C++ MC: uniform y in [0,1]^dim, domain indicator zeroes
   infeasible points. Fixed seed for reproducibility. *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Nsamp = 1000000;
exact = 1.*^-6;

(* ---- Build a vectorized flattened-integrand sampler from SectorData ---- *)
(* Returns {mean, sampleSigma, stderr} for one sector via uniform MC.       *)
sectorMC[sd_, n_] := Module[
  {dim, flat, pExps, pf, dc, ys, vals, logYpStar, feasible,
   logZ0, mp, ic},
  dim   = sd["Dimension"];
  flat  = sd["FlattenedPolys"];
  pExps = N[sd["PolynomialExponents"]];
  pf    = N[sd["Prefactor"]];
  dc    = Lookup[sd, "DomainConstraint", None];

  ys = RandomReal[{0, 1}, {n, dim}];

  (* integrand value at each sample row *)
  vals = Map[
    Function[y,
      Module[{polyVals},
        polyVals = Table[
          Total[Table[
            mono[[1]] * Exp[Total[mono[[2]] * Log[y]]],
            {mono, flat[[j]]}]],
          {j, Length[flat]}];
        pf * Times @@ MapThread[Exp[#2 * Log[#1]] &, {polyVals, pExps}]
      ]],
    ys];
  vals = Re[vals];

  (* domain indicator: zero out infeasible samples (matches generated C++) *)
  If[dc =!= None,
    logZ0 = N[dc["LogZ0"]]; mp = N[dc["MP"]]; ic = N[dc["IndicatorCoeffs"]];
    feasible = Map[Function[y, ((logZ0 - Total[ic * Log[y]])/mp) <= 0], ys];
    vals = MapThread[If[#2, #1, 0.] &, {vals, feasible}];
  ];

  Module[{mean, var},
    mean = Mean[vals];
    var  = Variance[vals];          (* per-sample variance over the cube *)
    <|"Mean" -> mean, "Sigma" -> Sqrt[var], "StdErr" -> Sqrt[var/n]|>
  ]
];

Print["================================================================"];
Print["  Toy 0:  I = Int_0^inf (1 + 10^6 x)^-2 dx   (exact = 1e-6)"];
Print["  N = ", Nsamp, " uniform samples per sector, seed 42"];
Print["================================================================"];

(* ============================ UNLIFTED ============================ *)
SeedRandom[42];
specU = <|
  "Polynomials"         -> {1 + 10^6 x[1]},
  "MonomialExponents"   -> {0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> {x[1]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;
(* 1D fan: two rays +1 and -1, each a maximal cone *)
fanU = {{{1}, {-1}}, {{1}, {2}}};

Print[""];
Print["---------------- WITHOUT lifting (plain tropical) ----------------"];
sdU = Table[ProcessSector[specU, fanU[[1]], {fanU[[2, s, 1]]}, s], {s, 2}];
resU = Table[sectorMC[sdU[[s]], Nsamp], {s, 2}];
Do[
  Print["  Sector ", s, "  prefactor=", N[sdU[[s]]["Prefactor"]],
        "  Qtilde=", sdU[[s]]["FlattenedPolys"][[1]]];
  Print["     mean=", resU[[s]]["Mean"],
        "  per-sample sigma=", resU[[s]]["Sigma"],
        "  sigma/mean=", resU[[s]]["Sigma"]/Abs[resU[[s]]["Mean"]]];,
  {s, 2}];
totU    = Total[#["Mean"] & /@ resU];
errU    = Sqrt[Total[(#["StdErr"]^2) & /@ resU]];
Print["  TOTAL (unlifted) = ", totU, " +/- ", errU,
      "   (rel.dev from exact = ", Abs[totU - exact]/exact, ")"];

(* ============================ LIFTED k=1 ============================ *)
SeedRandom[42];
lc = LiftCoefficients[specU,
  {<|"PolyIndex" -> 1, "ExponentVector" -> {1}, "k" -> 1|>}];
liftedSpec = lc["LiftedSpec"];
liftData   = lc["LiftData"];

raysK1  = {{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
sectsK1 = Table[{i, If[i < Length[raysK1], i+1, 1]}, {i, Length[raysK1]}];

Print[""];
Print["---------------- WITH lifting (k=1, z0=10^6) ----------------"];
Print["  Lifted poly: ", liftedSpec["Polynomials"][[1]], "   (coeff now O(1))"];
sdL = Table[
  ProcessSectorLifted[liftedSpec, raysK1, sectsK1[[s]], s, liftData],
  {s, Length[sectsK1]}];

liveTot = 0.; liveErr2 = 0.; nDrop = 0;
Do[
  If[AssociationQ[sdL[[s]]] && TrueQ[sdL[[s]]["EmptyDomain"]],
    nDrop++;
    Print["  Sector ", s, "  -> EmptyDomain (dropped, contributes 0)"];,
    Module[{r = sectorMC[sdL[[s]], Nsamp]},
      liveTot += r["Mean"]; liveErr2 += r["StdErr"]^2;
      Print["  Sector ", s, "  prefactor=", N[sdL[[s]]["Prefactor"]],
            "  Qtilde=", sdL[[s]]["FlattenedPolys"][[1]],
            "  dom=", If[sdL[[s]]["DomainConstraint"]===None,"None","constrained"]];
      Print["     mean=", r["Mean"],
            "  per-sample sigma=", r["Sigma"],
            "  sigma/mean=", r["Sigma"]/Abs[r["Mean"]]];
    ]
  ],
  {s, Length[sectsK1]}];
errL = Sqrt[liveErr2];
Print["  Dropped sectors: ", nDrop];
Print["  TOTAL (lifted)   = ", liveTot, " +/- ", errL,
      "   (rel.dev from exact = ", Abs[liveTot - exact]/exact, ")"];

Print[""];
Print["================================================================"];
Print["  COMPARISON at N = ", Nsamp];
Print["    unlifted stderr = ", errU];
Print["    lifted   stderr = ", errL];
Print["    error reduction factor (unlifted/lifted) = ", errU/errL];
Print["    variance reduction factor                = ", (errU/errL)^2];
Print["================================================================"];
