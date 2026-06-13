(* Verify the contour-rotation picture of complex-exponent flattening.
   Example:  I = Int_0^inf (1+x)^B dx,   B = -2 + I b,   exact = 1/(1 - I b).
   The large-x sector (ray +1) has complex effective exponent a = 1 - I b.
   Show: (i) a, prefactor 1/a, flattened poly from ProcessSector;
         (ii) pre-flatten integrand g(y)=y^{a-1}(1+y)^B oscillates (phase winds);
         (iii) flattened integrand f(y')=(1/a)(1+(y')^{1/a})^B is smooth;
         (iv) variance of Re/Im drops; values agree; total matches exact. *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

bval = 6;
Bc   = -2 + I bval;
exact = 1/(1 - I bval);

spec = <|
  "Polynomials"         -> {1 + x[1]},
  "MonomialExponents"   -> {0},
  "PolynomialExponents" -> {Bc},
  "Variables"           -> {x[1]},
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;
fan = {{{1}, {-1}}, {{1}, {2}}};   (* 1D: rays +1 and -1 *)

Print["================================================================"];
Print["  I = Int_0^inf (1+x)^(", Bc, ") dx     exact = ", N[exact]];
Print["================================================================"];

sd = Table[ProcessSector[spec, fan[[1]], {fan[[2, s, 1]]}, s], {s, 2}];

Do[
  Print["--- Sector ", s, "  ray=", fan[[1, fan[[2,s,1]] ]] , " ---"];
  Print["   a (NewExponents) = ", sd[[s]]["NewExponents"]];
  Print["   prefactor        = ", N[sd[[s]]["Prefactor"]]];
  Print["   FlattenedPolys   = ", sd[[s]]["FlattenedPolys"]];
  Print["   PolyExps (B)     = ", sd[[s]]["PolynomialExponents"]];,
  {s, 2}];

(* ----- focus on the complex sector (ray +1 = sector 1) ----- *)
sCx = 1;
aP  = sd[[sCx]]["NewExponents"][[1]];
pf  = sd[[sCx]]["Prefactor"];
Print[""];
Print["Complex sector: a = ", aP, "   (Im a = ", Im[aP], ")"];

(* spiral: image points y = (y')^(1/a) as y' runs 1 -> 0 *)
Print[""];
Print["Contour image y = (y')^(1/a) for y' = 1, .5, .1, .01, .001:"];
Do[Print["   y' = ", yp, "  ->  y = ", N[yp^(1/aP)],
         "   (|y|=", N[Abs[yp^(1/aP)]], ", arg=", N[Arg[yp^(1/aP)]], ")"],
   {yp, {1, 0.5, 0.1, 0.01, 0.001}}];

(* pre-flatten (real-axis) integrand g(y) and flattened f(y') for sector + *)
g[y_]  := y^(aP - 1) * (1 + y)^Bc;                 (* on real y in (0,1] *)
f[yp_] := pf * (1 + yp^(1/aP))^Bc;                 (* flattened, y' in (0,1] *)

(* sanity: integrals must match each other and be finite *)
IgReal = NIntegrate[g[y], {y, 0, 1}, MaxRecursion -> 40, WorkingPrecision -> 20];
IfFlat = NIntegrate[f[yp], {yp, 0, 1}, MaxRecursion -> 40, WorkingPrecision -> 20];
Print[""];
Print["Sector+ value, real contour   Int_0^1 g  = ", N[IgReal]];
Print["Sector+ value, rotated contour Int_0^1 f = ", N[IfFlat]];

(* total over both sectors via the flattened integrands the pipeline samples *)
gOther[y_] := (sd[[2]]["Prefactor"]) *
   (1 + y^(1/sd[[2]]["NewExponents"][[1]]))^Bc;
ItotFlat = IfFlat + NIntegrate[gOther[y], {y,0,1}, MaxRecursion->40, WorkingPrecision->20];
Print["Total (both flattened sectors)           = ", N[ItotFlat]];
Print["Exact 1/(1 - I b)                        = ", N[exact]];
Print["Relative error                           = ", N[Abs[(ItotFlat-exact)/exact]]];

(* ----- variance seen by uniform MC: real contour vs rotated contour ----- *)
SeedRandom[42];
ys  = RandomReal[{0,1}, 200000];
gv  = g /@ ys;     (* would-be samples on the real (un-rotated) contour *)
fv  = f /@ ys;     (* actual flattened samples *)
Print[""];
Print["Per-sample statistics over 2e5 uniform points (sector +):"];
Print["   real contour  g : mean=", N[Mean[gv]],
      "  Var[Re]=", N[Variance[Re[gv]]], "  Var[Im]=", N[Variance[Im[gv]]]];
Print["   rotated (flat) f: mean=", N[Mean[fv]],
      "  Var[Re]=", N[Variance[Re[fv]]], "  Var[Im]=", N[Variance[Im[fv]]]];
Print["   |g| range = ", {N@Min@Abs@gv, N@Max@Abs@gv},
      "    |f| range = ", {N@Min@Abs@fv, N@Max@Abs@fv}];
Print["   Var reduction Re: ", N[Variance[Re[gv]]/Variance[Re[fv]]],
      "   Im: ", N[Variance[Im[gv]]/Variance[Im[fv]]]];
