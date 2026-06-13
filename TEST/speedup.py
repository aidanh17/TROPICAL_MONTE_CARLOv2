#!/usr/bin/env python3
"""
Speed-up to reach a TARGET relative accuracy, vs plain MC.

For each method we take the *measured* smallest swept evaluation count that hits
rel.err <= eps (no sub-budget extrapolation -> no artifacts).  Plain MC almost
never reaches the small targets within the sweep, so MC's required cost is
extrapolated from its best point at the *theoretical* 1/sqrt(N) rate:
    N_MC(eps) = N0_MC * (err0_MC / eps)^2 .
Speed-up = N_MC(eps) / N_method(eps).
  >X  : method already met the target at its smallest budget -> X is a lower bound
  --  : method never reached the target within the swept budget
"""
import csv, math, os
from collections import defaultdict

HERE=os.path.dirname(os.path.abspath(__file__))
rows=[]
with open(os.path.join(HERE,"INTERFILES","results.csv")) as f:
    for r in csv.DictReader(f):
        for k in ("dim","nsec","budget","neval"): r[k]=int(float(r[k]))
        for k in ("wall_s","rel_err"): r[k]=float(r[k])
        rows.append(r)

G=defaultdict(list)
for r in rows:
    if r["case"]!="canary": G[(r["case"],r["method"])].append(r)
for k in G: G[k].sort(key=lambda r:r["neval"])

METHODS=["QMC","VEGAS","SUAVE","DIVONNE","CUHRE"]

def measured_cross(g, eps, key="neval"):
    """smallest swept cost hitting err<=eps; (cost, at_first_point?) or None."""
    for i,r in enumerate(g):
        if r["rel_err"]<=eps:
            return r[key], (i==0)
    return None

def mc_extrap(g, eps, key="neval"):
    """MC cost to reach eps via theoretical p=0.5, anchored at the LARGEST
    budget (statistically most reliable; MC error is noisy point-to-point)."""
    if not g: return None
    anc=g[-1]                       # largest neval
    factor=(anc["rel_err"]/eps)**2.0
    return anc[key]*factor

cases=[]; seen=set()
for r in rows:
    if r["case"]!="canary" and r["case"] not in seen:
        seen.add(r["case"]); cases.append((r["case"],r["dim"]))

def report(eps, key, label):
    print(f"\n================  {label} speed-up vs plain MC, rel.err = {eps:g}  ================")
    print(f"{'case':<20}{'dim':>4} " + "".join(f"{m:>11}" for m in METHODS))
    geo=defaultdict(list)
    for case,dim in cases:
        cmc=mc_extrap(G[(case,"MC")],eps,key)
        line=f"{case:<20}{dim:>4} "
        for m in METHODS:
            mc=measured_cross(G[(case,m)],eps,key)
            if cmc and mc:
                cost,first=mc; su=cmc/cost
                geo[m].append(su)
                s=f"{su:,.0f}" if su>=1 else f"{su:.2f}"
                line+=f"{('>'+s) if first else s:>11}"
            else:
                line+=f"{'--':>11}"
        print(line)
    print(f"{'GEOMEAN (x)':<24} " + "".join(
        (f"{math.exp(sum(map(math.log,geo[m]))/len(geo[m])):>11,.0f}" if geo[m] else f"{'--':>11}")
        for m in METHODS))
    return geo

report(1e-3,"neval","EVALUATIONS")
report(1e-4,"neval","EVALUATIONS")
report(1e-3,"wall_s","WALL-TIME")
print("\n>X = lower bound (target met at the method's smallest budget); "
      "-- = not reached within the sweep.")
print("MC baseline extrapolated at the theoretical 1/sqrt(N); method costs are "
      "the measured crossing budgets.")
