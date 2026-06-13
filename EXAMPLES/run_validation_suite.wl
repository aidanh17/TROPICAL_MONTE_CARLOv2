(* ============================================================================
   run_validation_suite.wl

   Runs Examples 1-14 (by loading tropical_eval_examples.wl, which executes
   them inline and defines the RunTest1..7 / RunAllTests functions) and then
   invokes the Module-5 validation suite RunAllTests[].

   tropical_eval_examples.wl DEFINES RunAllTests[] but never CALLS it, so the
   6-test suite (Tests 1, 2, 3v2, 5, 6, 7) is never exercised by running that
   file directly.  This runner closes that gap.

     wolframscript -file EXAMPLES/run_validation_suite.wl
   ============================================================================ *)

SetDirectory[DirectoryName[$InputFileName]];

(* Get runs Examples 1-14 inline and brings RunAllTests[] / RunTest*[] into
   scope.  $InputFileName inside the Get resolves to the examples file, so its
   own relative package load works. *)
Get[FileNameJoin[{Directory[], "tropical_eval_examples.wl"}]];

Print[];
Print["################################################################"];
Print["#  Invoking RunAllTests[] (Module-5 validation suite)          #"];
Print["################################################################"];

results = RunAllTests[];

Print[];
Print["Per-test results: ", results];

(* Exit nonzero if any test failed, so the runner is CI-usable. *)
nFail = Count[results, {_, False}];
If[nFail > 0,
  Print[">>> VALIDATION SUITE: ", nFail, " test(s) FAILED"],
  Print[">>> VALIDATION SUITE: ALL TESTS PASSED"]
];
