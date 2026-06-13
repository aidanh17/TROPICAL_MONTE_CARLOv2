#!/usr/bin/env python3
"""Make convergence + dimension-scaling plots from INTERFILES/results.csv."""
import csv, math, os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV  = os.path.join(HERE, "INTERFILES", "results.csv")
rows = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        for k in ("dim","nsec","budget","neval"): r[k]=int(float(r[k]))
        for k in ("wall_s","reported_err","abs_err","rel_err"): r[k]=float(r[k])
        rows.append(r)

METHODS = ["MC","QMC","VEGAS","SUAVE","DIVONNE","CUHRE"]
STYLE = {"MC":("o","-"),"QMC":("s","-"),"VEGAS":("^","-"),
         "SUAVE":("v","--"),"DIVONNE":("D","--"),"CUHRE":("*","-")}

def curve(case, method, x="neval", y="rel_err"):
    pts=[(r[x],r[y]) for r in rows if r["case"]==case and r["method"]==method and r[y]>0]
    pts.sort()
    return [p[0] for p in pts],[p[1] for p in pts]

# ---- (1) convergence for a few representative cases ----
panel = ["A_simplex_n6","B_quad_n6","Af_fracsimplex_n4","D_product_n4"]
fig,axes=plt.subplots(2,2,figsize=(12,9))
for ax,case in zip(axes.flat,panel):
    for m in METHODS:
        xs,ys=curve(case,m)
        if xs:
            mk,ls=STYLE[m]
            ax.loglog(xs,ys,marker=mk,linestyle=ls,ms=5,label=m)
    # reference slopes
    if case=="A_simplex_n6":
        import numpy as np
        n0=2e4; ax.loglog([n0,n0*100],[2e-2,2e-2/10],'k:',lw=1)
        ax.text(n0*3,1.1e-2,"1/sqrt(N)",fontsize=8)
        ax.loglog([n0,n0*100],[2e-3,2e-3/100],'k-.',lw=1)
        ax.text(n0*3,3e-4,"1/N",fontsize=8)
    dim=[r["dim"] for r in rows if r["case"]==case][0]
    ax.set_title(f"{case}  (n={dim})")
    ax.set_xlabel("integrand evaluations"); ax.set_ylabel("relative error vs exact")
    ax.grid(True,which="both",alpha=0.3); ax.legend(fontsize=8)
fig.suptitle("Integrating the FLATTENED tropical sectors: error vs cost",fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(HERE,"INTERFILES","convergence.png"),dpi=130)
print("wrote convergence.png")

# ---- (2) dimension scaling for the simplex family at fixed per-sector budget ----
B=100000
fig2,(a1,a2)=plt.subplots(1,2,figsize=(13,5))
dims=sorted({r["dim"] for r in rows if r["family"]=="simplex"})
for m in METHODS:
    ds=[];es=[];ts=[]
    for d in dims:
        rr=[r for r in rows if r["family"]=="simplex" and r["dim"]==d
            and r["method"]==m and r["budget"]==B]
        if rr: ds.append(d); es.append(rr[0]["rel_err"]); ts.append(rr[0]["wall_s"])
    mk,ls=STYLE[m]
    if ds: a1.semilogy(ds,es,marker=mk,linestyle=ls,label=m)
    if ds: a2.plot(ds,ts,marker=mk,linestyle=ls,label=m)
a1.set_title(f"Accuracy vs dimension (simplex, {B:.0e} samples/sector)")
a1.set_xlabel("dimension n"); a1.set_ylabel("relative error vs exact")
a1.grid(True,which="both",alpha=0.3); a1.legend(fontsize=8)
a2.set_title(f"Wall time vs dimension (same budget)")
a2.set_xlabel("dimension n"); a2.set_ylabel("wall time (s)")
a2.grid(True,alpha=0.3); a2.legend(fontsize=8)
fig2.tight_layout()
fig2.savefig(os.path.join(HERE,"INTERFILES","dimension_scaling.png"),dpi=130)
print("wrote dimension_scaling.png")

# ---- (3) canary ----
fig3,ax=plt.subplots(figsize=(7,5))
for d in sorted({r["dim"] for r in rows if r["case"]=="canary"}):
    for m,c in (("MC","C0"),("QMC","C2")):
        xs,ys=[],[]
        for r in rows:
            if r["case"]=="canary" and r["dim"]==d and r["method"]==m and r["rel_err"]>0:
                xs.append(r["neval"]); ys.append(r["rel_err"])
        if xs:
            o=sorted(zip(xs,ys)); xs=[p[0] for p in o]; ys=[p[1] for p in o]
            ax.loglog(xs,ys,("o-" if m=="MC" else "s-"),color=c,alpha=0.4+0.07*d,ms=4)
ax.set_title("QMC canary: Sobol (green) vs pseudo-random MC (blue), dims 2-8\n"
             "smooth exp-sum integrand, exact known")
ax.set_xlabel("evaluations"); ax.set_ylabel("relative error")
ax.grid(True,which="both",alpha=0.3)
fig3.tight_layout(); fig3.savefig(os.path.join(HERE,"INTERFILES","canary.png"),dpi=130)
print("wrote canary.png")
