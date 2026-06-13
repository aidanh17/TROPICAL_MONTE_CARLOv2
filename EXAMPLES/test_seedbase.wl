(* ============================================================================
   test_seedbase.wl

   Non-trivial test of the SeedBase RUNTIME-ARGV override (SUMMARY §4,
   GenerateCppMonteCarlo / EvaluateTropicalMC docs):

       tropical_mc <in> <out> [n_samples] [n_threads] [seed_base]

   The optional 5th argv overrides the compile-time SeedBase, and the
   per-kinematic-point RNG seed is seed_base + kp.  This lets a caller
   re-run the SAME compiled binary under different seeds WITHOUT recompiling.
   No EXAMPLES-level test exercised this contract before (only SANDBOX harness).

   Integrand:  Int_[0,Inf)^2 dx1 dx2 / (1 + x1^2 + x2^2)^3 = Pi/8  (no kinematics).

   Asserts:
     (A) determinism   : same seed_base => byte-identical result (run twice)
     (B) override wired : omitting argv[5] == passing argv[5] = compile-time 42
     (C) distinct stream: different seed_base => different MC estimate
     (D) correctness    : every stream is within 5 sigma of Pi/8

     wolframscript -file EXAMPLES/test_seedbase.wl
   ============================================================================ *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["=== Test SeedBase: runtime seed override of one compiled binary ==="];
Print["Int dx1 dx2 / (1 + x1^2 + x2^2)^3 = Pi/8 = ", N[Pi/8]];
Print[];

passAll = True;

Module[{poly, vars, spec, verts, fanData, dualVerts, simplices, sectors,
        cppFile, binary, kinFile, nSamples = 400000, exact = N[Pi/8],
        run, rNo5, r42, r42b, r99, r12345, within},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
           "PolynomialExponents" -> {-3}, "Variables" -> vars,
           "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;

  verts   = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dualVerts, simplices} = fanData;
  sectors = Select[
    Table[ProcessSector[spec, dualVerts, simplices[[s]], s], {s, Length[simplices]}],
    AssociationQ];
  Print["Fan: ", Length[dualVerts], " rays, ", Length[simplices], " sectors"];

  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile = FileNameJoin[{Directory[], "INTERFILES", "seedtest_mc_generated.cpp"}];
  binary  = FileNameJoin[{Directory[], "INTERFILES", "seedtest_mc"}];
  kinFile = FileNameJoin[{Directory[], "INTERFILES", "seedtest_kin.txt"}];

  (* Compile ONCE.  Default compile-time SeedBase = 42. *)
  GenerateCppMonteCarlo[sectors, {}, spec, cppFile, "NSamples" -> nSamples, "SeedBase" -> 42];
  If[CompileCpp[cppFile, binary, False] === $Failed,
    Print["Compilation FAILED"]; passAll = False; Abort[]];
  Print["Compiled once: ", FileNameTake[binary]];
  Export[kinFile, "1\n", "Text"];

  (* run[seedArgOrNone] -> {re, im, reErr, imErr}.  Reuses the SAME binary. *)
  run[seedArg_] := Module[{out, proc, tbl},
    out = FileNameJoin[{Directory[], "INTERFILES",
            "seedtest_res_" <> ToString[seedArg] <> ".txt"}];
    proc = RunProcess[Join[{binary, kinFile, out, ToString[nSamples], "4"},
             If[seedArg === None, {}, {ToString[seedArg]}]]];
    If[proc["ExitCode"] =!= 0, Print["  run FAILED: ", proc["StandardError"]]; Return[$Failed]];
    tbl = Import[out, "Table"];
    tbl[[1]]
  ];

  rNo5  = run[None];   (* no argv[5] -> compile-time 42 *)
  r42   = run[42];
  r42b  = run[42];
  r99   = run[99];
  r12345= run[12345];

  Print[];
  Print["  no-argv5 : ", rNo5];
  Print["  seed 42  : ", r42];
  Print["  seed 42' : ", r42b];
  Print["  seed 99  : ", r99];
  Print["  seed12345: ", r12345];
  Print[];

  (* (A) determinism: seed 42 twice identical *)
  If[r42 === r42b,
    Print["  (A) determinism (seed 42 == seed 42'): PASS"],
    Print["  (A) determinism: FAIL  ", r42, " vs ", r42b]; passAll = False];

  (* (B) override wired: omitting argv[5] equals passing 42 *)
  If[rNo5 === r42,
    Print["  (B) argv[5] override (no-arg == 42): PASS"],
    Print["  (B) argv[5] override: FAIL  ", rNo5, " vs ", r42]; passAll = False];

  (* (C) distinct streams: 42, 99, 12345 give different Re *)
  If[r42[[1]] =!= r99[[1]] && r42[[1]] =!= r12345[[1]] && r99[[1]] =!= r12345[[1]],
    Print["  (C) distinct streams (42/99/12345 differ): PASS"],
    Print["  (C) distinct streams: FAIL (seeds gave identical Re)"]; passAll = False];

  (* (D) correctness: each stream within 5 sigma of Pi/8 *)
  within[r_] := r[[3]] > 0 && Abs[r[[1]] - exact] <= 5 r[[3]];
  If[AllTrue[{rNo5, r42, r99, r12345}, within],
    Print["  (D) all streams within 5 sigma of Pi/8: PASS"],
    Print["  (D) correctness: FAIL  devs(sigma) = ",
      (Abs[#[[1]] - exact]/#[[3]]) & /@ {rNo5, r42, r99, r12345}]; passAll = False];
];

Print[];
Print["================================================================"];
Print["  SeedBase runtime-override test: ", If[passAll, "PASS", "FAIL"]];
Print["================================================================"];
