#!/bin/bash
# ============================================================================
# TEST/run.sh
# Compile bench.cpp once per generated sector header (linking Boost + CUBA),
# run the budget sweep for each case, and collect all CSV rows into
# INTERFILES/results.csv.  Single-threaded throughout for fair timing.
# ============================================================================
set -u
cd "$(dirname "$0")"
PREFIX=/opt/homebrew
ID=INTERFILES
OUT=$ID/results.csv
LOG=$ID/run.log
: > "$LOG"

# header written once, from the first case
HEADER_DONE=0
: > "$OUT"

# choose a per-sector budget sweep based on dimension (bound wall time at high n)
budgets_for_dim() {
  local d=$1
  if   [ "$d" -le 4 ]; then echo "1000 3162 10000 31623 100000 316228 1000000"
  elif [ "$d" -le 6 ]; then echo "1000 3162 10000 31623 100000 316228"
  else                      echo "1000 3162 10000 31623 100000"
  fi
}

# iterate cases via the manifest (name,dim,nsectors,ref_re,ref_im,family)
# manifest fields are quoted by Mathematica's CSV export; strip the quotes
tail -n +2 "$ID/manifest.csv" | sed 's/"//g' | while IFS=, read -r name dim nsec rre rim fam; do
  [ -z "$name" ] && continue
  hdr="$ID/sectors_${name}.hpp"
  bin="$ID/bench_${name}"
  if [ ! -f "$hdr" ]; then echo "MISSING $hdr" | tee -a "$LOG"; continue; fi

  echo "=== compiling $name (dim=$dim, nsec=$nsec) ===" | tee -a "$LOG"
  g++ -O3 -std=c++17 -funroll-loops \
      -I"$PREFIX/include" \
      -DSECTOR_HEADER="\"$hdr\"" \
      bench.cpp \
      -L"$PREFIX/lib" -lcuba -lm -o "$bin" 2>>"$LOG"
  if [ $? -ne 0 ]; then echo "COMPILE FAILED: $name (see $LOG)"; continue; fi

  B=$(budgets_for_dim "$dim")
  echo "=== running   $name   budgets: $B ===" | tee -a "$LOG"
  if [ "$HEADER_DONE" -eq 0 ]; then
    "$bin" $B 2>>"$LOG" | tee -a "$OUT" >/dev/null
    HEADER_DONE=1
  else
    # drop the CSV header line from subsequent binaries
    "$bin" $B 2>>"$LOG" | tail -n +2 >> "$OUT"
  fi
  echo "    done $name" | tee -a "$LOG"
done

echo "ALL DONE -> $OUT"
wc -l "$OUT"
