(* Contour-rotation figure for the manual's complex-exponent example.
   Sector + of  I = Int_0^inf (1+x)^(-2+6I) dx,  a = 1 - 6 I.
   (A) flattening y=(y')^{1/a} maps the real segment [0,1] onto a log spiral.
   (B) what uniform MC sees vs its sample variable: un-rotated Re g oscillates
       with diverging frequency near 0; rotated (flattened) Re f is smooth. *)

b  = 6;
Bc = -2 + I b;
a  = 1 - I b;

g[y_]  := y^(a - 1) (1 + y)^Bc;       (* real-axis (un-rotated) integrand *)
f[yp_] := (1/a) (1 + yp^(1/a))^Bc;    (* rotated / flattened integrand    *)

red  = RGBColor[0.78, 0.15, 0.15];
lbl  = Directive[FontFamily -> "Times", FontSize -> 13, Black];

(* --- Panel A: parametrize by s = -log y' so the spiral winds visibly --- *)
yof[s_] := Exp[-(s/37) (1 + 6 I)];     (* = (y')^{1/a} with y'=e^{-s} *)
spiralPts = Table[{Re[#], Im[#]} &@yof[s], {s, 0, 60, 0.05}];

panelA = Show[
  Graphics[{
    (* unit circle + axes for reference *)
    Directive[Gray, Dotted], Circle[{0, 0}, 1],
    Directive[GrayLevel[0.5], Thickness[0.002]],
      Line[{{-1.15, 0}, {1.2, 0}}], Line[{{0, -1.2}, {0, 0.35}}],
    (* real integration segment [0,1] *)
    Directive[Gray, Dashed, Thickness[0.004]], Line[{{0, 0}, {1, 0}}],
    Text[Style["real segment [0,1]", lbl, Gray], {0.5, 0}, {0, -1.5}],
    (* rotated contour (spiral) -- solid *)
    Directive[red, Dashing[None], Thickness[0.006]], Line[spiralPts],
    (* markers *)
    Black, PointSize[0.022], Point[{1, 0}], Point[{0, 0}],
    Text[Style["y=1  (y'=1)", lbl], {1.02, 0}, {-1.05, 1.6}],
    Text[Style["y\[Rule]0  (y'\[Rule]0)", lbl], {0, 0}, {-1.15, 1.6}]
  }],
  PlotLabel -> Style["(A)  contour rotation onto a log spiral", lbl],
  Axes -> False, AspectRatio -> Automatic,
  PlotRange -> {{-1.18, 1.55}, {-1.22, 0.5}}, ImageSize -> 360];

(* --- Panel B: integrand vs uniform sample variable on [0,1] --- *)
uv = Subdivide[0.0004, 1, 2500];
panelB = ListLinePlot[
   {Transpose[{uv, Re[g /@ uv]}], Transpose[{uv, Re[f /@ uv]}]},
   PlotStyle -> {Directive[GrayLevel[0.55], Thickness[0.003]],
                 Directive[red, Thickness[0.006]]},
   PlotLegends -> Placed[
     {Style["before: Re g (real axis)", lbl],
      Style["after: Re f (rotated)", lbl]}, Below],
   Frame -> True, Axes -> False,
   FrameLabel -> {Style["uniform sample variable", lbl],
                  Style["Re integrand", lbl]},
   PlotLabel -> Style["(B)  what the Monte Carlo samples", lbl],
   LabelStyle -> lbl, ImageSize -> 350, AspectRatio -> 1,
   PlotRange -> {{0, 1}, {-1.05, 1.05}}];

fig = GraphicsRow[{panelA, panelB}, Spacings -> 20, ImageSize -> 730];
Export["/Users/aidanh/Desktop/TROPICAL_MONTE_CARLO2/TROPICAL_MONTE_CARLOv2/MANUAL/complex_flatten_spiral.pdf", fig];
Export["/tmp/complex_flatten_spiral.png", fig, ImageResolution -> 130];
Print["figure written; spiral |y| range = ", {Min[Norm /@ spiralPts], Max[Norm /@ spiralPts]}];
