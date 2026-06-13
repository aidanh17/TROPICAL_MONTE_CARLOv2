#!/usr/bin/env python3
"""
TEST/analyze.py  --  digest INTERFILES/results.csv

Produces, on stdout:
  * a per-case table of achieved relative error + wall time at the top budget
    for each method,
  * empirical convergence-rate exponents  p  in  rel_err ~ C * neval^(-p)
    (fit by least squares on log-log), for MC / QMC / each CUBA routine,
  * the QMC-canary rate check (Sobol must beat MC's 0.5 slope),
  * a "time to reach 1e-3 / 1e-4 relative error" summary,
  * dimension-scaling tables (rate & time vs dimension) for the simplex family.
Writes a compact markdown summary to INTERFILES/summary.md as well.
"""
import csv, math, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CSV  = os.path.join(HERE, "INTERFILES", "results.csv")

rows = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        for k in ("dim","nsec","budget","neval"):
            r[k] = int(float(r[k]))
        for k in ("wall_s","est_re","est_im","reported_err","abs_err","rel_err"):
            r[k] = float(r[k])
        rows.append(r)

# group: (case, method) -> list of rows sorted by neval
G = defaultdict(list)
for r in rows:
    G[(r["case"], r["method"])].append(r)
for k in G: G[k].sort(key=lambda r: r["neval"])

cases   = []
seen=set()
for r in rows:
    if r["case"] != "canary" and r["case"] not in seen:
        seen.add(r["case"]); cases.append((r["case"], r["dim"], r["nsec"], r["family"]))

METHODS = ["MC","QMC","VEGAS","SUAVE","DIVONNE","CUHRE"]

def fit_slope(pts):
    """least-squares slope of log10(rel_err) vs log10(neval); returns -p."""
    pts = [(n,e) for (n,e) in pts if e>0 and n>0]
    if len(pts) < 3: return float('nan')
    xs=[math.log10(n) for n,_ in pts]; ys=[math.log10(e) for _,e in pts]
    mx=sum(xs)/len(xs); my=sum(ys)/len(ys)
    num=sum((x-mx)*(y-my) for x,y in zip(xs,ys)); den=sum((x-mx)**2 for x in xs)
    return num/den if den else float('nan')

def time_to(case, method, target):
    """interpolate wall time to reach rel_err<=target on the (neval,err) curve."""
    pts=[(r["neval"], r["rel_err"], r["wall_s"]) for r in G[(case,method)]]
    for n,e,t in pts:
        if e<=target: return t
    return None

out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)

# ---------------------------------------------------------------- canary
P("\n================ QMC CANARY (smooth exp-sum, exact known) ================")
P("Confirms the Sobol generator delivers QMC convergence (slope ~ -1) vs MC (~ -0.5).")
P(f"{'dim':>4} {'MC slope':>10} {'QMC slope':>10}   (slope = d log10(relerr)/d log10 N)")
cdims=sorted({r['dim'] for r in rows if r['case']=='canary'})
for d in cdims:
    mc =[(r['neval'],r['rel_err']) for r in rows if r['case']=='canary' and r['method']=='MC'  and r['dim']==d]
    qmc=[(r['neval'],r['rel_err']) for r in rows if r['case']=='canary' and r['method']=='QMC' and r['dim']==d]
    P(f"{d:>4} {fit_slope(mc):>10.2f} {fit_slope(qmc):>10.2f}")

# ---------------------------------------------------------------- per case
P("\n================ PER-CASE: rel.err & wall(s) at TOP budget ================")
hdr = f"{'case':<20}{'dim':>4}{'nsec':>5} " + "".join(f"{m:>22}" for m in METHODS)
P(hdr)
P(f"{'':<29}" + "".join(f"{'relerr':>11}{'sec':>11}" for _ in METHODS))
for case,dim,nsec,fam in cases:
    line=f"{case:<20}{dim:>4}{nsec:>5} "
    for m in METHODS:
        g=G[(case,m)]
        if g:
            r=g[-1]; line+=f"{r['rel_err']:>11.2e}{r['wall_s']:>11.3f}"
        else: line+=f"{'-':>11}{'-':>11}"
    P(line)

# ---------------------------------------------------------------- slopes
P("\n================ CONVERGENCE RATE p  (rel_err ~ neval^-p) ================")
P(f"{'case':<20}{'dim':>4} " + "".join(f"{m:>9}" for m in METHODS))
for case,dim,nsec,fam in cases:
    line=f"{case:<20}{dim:>4} "
    for m in METHODS:
        pts=[(r["neval"],r["rel_err"]) for r in G[(case,m)]]
        s=fit_slope(pts)
        line+=f"{(-s):>9.2f}" if s==s else f"{'-':>9}"
    P(line)

# ---------------------------------------------------------------- time-to-acc
for target in (1e-3, 1e-4):
    P(f"\n================ WALL TIME (s) to reach rel.err <= {target:g} ================")
    P("(blank = not reached within the swept budget)")
    P(f"{'case':<20}{'dim':>4} " + "".join(f"{m:>9}" for m in METHODS))
    for case,dim,nsec,fam in cases:
        line=f"{case:<20}{dim:>4} "
        for m in METHODS:
            t=time_to(case,m,target)
            line+=f"{t:>9.3f}" if t is not None else f"{'--':>9}"
        P(line)

# ---------------------------------------------------------------- dim scaling
P("\n================ DIMENSION SCALING (simplex family A) ================")
P("rel.err at the largest COMMON budget, per dimension:")
A=[(c,d,ns,fam) for (c,d,ns,fam) in cases if fam=="simplex"]
A.sort(key=lambda t:t[1])
# largest budget common to all simplex cases
budsets=[set(r["budget"] for r in rows if r["case"]==c) for c,_,_,_ in A]
common=set.intersection(*budsets) if budsets else set()
cb=max(common) if common else None
P(f"(common budget per sector = {cb})")
P(f"{'dim':>4} " + "".join(f"{m:>11}" for m in METHODS))
for c,d,ns,fam in A:
    line=f"{d:>4} "
    for m in METHODS:
        rr=[r for r in G[(c,m)] if r["budget"]==cb]
        line+=f"{rr[0]['rel_err']:>11.2e}" if rr else f"{'-':>11}"
    P(line)
P("\nwall time (s) at that budget, per dimension:")
P(f"{'dim':>4} " + "".join(f"{m:>11}" for m in METHODS))
for c,d,ns,fam in A:
    line=f"{d:>4} "
    for m in METHODS:
        rr=[r for r in G[(c,m)] if r["budget"]==cb]
        line+=f"{rr[0]['wall_s']:>11.3f}" if rr else f"{'-':>11}"
    P(line)

with open(os.path.join(HERE,"INTERFILES","summary.txt"),"w") as f:
    f.write("\n".join(out))
print("\n[written INTERFILES/summary.txt]")
