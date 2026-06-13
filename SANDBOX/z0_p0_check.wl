(* ============================================================================
   z0_p0_check.wl  --  P0 acceptance: reproduce Case A's 17.95x and gate metrics
   (benchmark_results.md) through the new sweepK harness.

   Set Z0_P0_FAST=1 in the environment for a quick shake-out (small N, 2 seeds).
   ============================================================================ *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "z0_sweep_common.wl"}]];

If[Environment["Z0_P0_FAST"] === "1",
  $Z0NSamples = 200000; $Z0Seeds = {42, 101};
  Print["[FAST MODE] N=", $Z0NSamples, " seeds=", $Z0Seeds]
];

caseA = <|
  "CaseId" -> "CaseA", "Tier" -> "P0",
  "Spec" -> <|
    "Polynomials"         -> {1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2},
    "MonomialExponents"   -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"           -> {x[1], x[2]},
    "KinematicSymbols"    -> {},
    "RegulatorSymbol"     -> None|>,
  "ProxyPoly"   -> 1 + x[1]^2 + x[2]^2 + x[1] x[2]^2,
  "PrimaryRule" -> <|"PolyIndex" -> 1, "ExponentVector" -> {2, 0}|>,
  "Magnitude"   -> 10^6, "Alpha" -> {2, 0},
  "FanProvider" -> Function[{k, ls}, tryAutoFan[ls]],
  "KList"       -> {2}|>;

Print["=== Case A unlifted baseline ==="];
baseA = unliftedBaseline[caseA];
Print["  unlifted sigmaSample(mean) = ", baseA["sigmaSampleMean"]];
Print["  unlifted sigmaTrue (strict)= ", baseA["sigmaTrue"]];
Print["  unlifted sigmaTrueConv     = ", baseA["sigmaTrueConv"]];
Print["  unlifted convMassFrac      = ", baseA["convMassFrac"]];
Print["  unlifted value(mean)       = ", baseA["valueMean"]];
Print["  unlifted nSectors          = ", baseA["nSectors"]];
Print["  truth                      = ", baseA["truth"]];

Print["\n=== Case A sweepK (k=2) ==="];
rowsA = sweepK[caseA, baseA];

(* k-level summary from the seed rows *)
Module[{k2, sampL, vrSample, vrTrue, vrTrueConv, gateOK, hcAny, nDrop, nSec,
        sigT, sigTC, nDiv, domSOM},
  k2 = Select[rowsA, #["k"] == 2 && NumericQ[#["sigmaSample"]] &];
  sampL = Mean[#["sigmaSample"] & /@ k2];
  sigT  = First[#["sigmaTrue"] & /@ k2];
  sigTC = First[#["sigmaTrueConv"] & /@ k2];
  nDiv  = First[#["nDivergent"] & /@ k2];
  domSOM= First[#["domSigOverMu"] & /@ k2];
  vrSample = (baseA["sigmaSampleMean"] / sampL)^2;
  vrTrue     = First[#["varRedTrue"] & /@ k2];
  vrTrueConv = First[#["varRedTrueConv"] & /@ k2];
  gateOK = First[#["gateOK"] & /@ k2];
  hcAny  = First[#["hasConstAny"] & /@ k2];
  nDrop  = First[#["nDropped"] & /@ k2];
  nSec   = First[#["nSectors"] & /@ k2];

  Print["  lifted k=2 sigmaSample(mean) = ", sampL];
  Print["  lifted k=2 sigmaTrue (strict)= ", sigT, "   nDivergent=", nDiv];
  Print["  lifted k=2 sigmaTrueConv     = ", sigTC];
  Print["  lifted k=2 dom sigma/mu      = ", domSOM];
  Print["  nSectors=", nSec, "  nDropped=", nDrop,
        "  survivors=", nSec - nDrop];
  Print["  gateOK=", gateOK, "  hasConstAny=", hcAny];
  Print[""];
  Print["  >>> sigma_sample  VarRed = ", N[vrSample, 4],
        "   (benchmark target ~ 17.95)"];
  Print["  >>> sigma_true    VarRed = ", vrTrue, "  (strict; 0 if heavy-tail)"];
  Print["  >>> sigma_true(conv) VarRed = ", vrTrueConv,
        "  (convergent sectors; phase1 reported ~4.2)"];
  Print[""];
  Print["  P0 sample-VR in [12,26]?  ",
        If[NumericQ[vrSample] && 12 <= vrSample <= 26, "PASS", "CHECK"]];
];

(* per-seed distinctness check (seed control is real, plan 6) *)
Module[{vals},
  vals = DeleteDuplicates[
    Cases[#["value"] & /@ Select[rowsA, #["k"] == 2 &], _?NumericQ]];
  Print["\n  distinct MC values across seeds: ", Length[vals],
        " of ", Length[$Z0Seeds],
        "  -> ", If[Length[vals] == Length[$Z0Seeds],
          "seed control VERIFIED", "WARNING: seeds not distinct"]];
  Print["  values: ", vals];
];

writeRows[rowsA, FileNameJoin[{DirectoryName[$InputFileName], "..", "AUXT",
  "z0_p0_check.csv"}]];
Print["\nz0_p0_check.wl COMPLETE"];
