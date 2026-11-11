"""Per-replica opening frequency of residue 14, recomputed from the SASA arrays.

The time-averaged exposure the first draft reported is an average over an open-and-shut
residue: the median is zero in all twelve systems. The quantity that actually separates the
groups is how often the residue is open at all, so it is computed here per MD replica - the
unit of analysis - and checked against the group contrasts quoted in the manuscript.
"""
import glob
import os
import re
import numpy as np

SUF = ["wt", "g", "c", "s", "t", "v"]
GAM = {"wt": 0, "g": 0, "c": 1, "s": 1, "t": 2, "v": 2}
LAB = {"wt": "WT", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}
TAG = {"ILDR2": {s: "i" + s for s in SUF},
       "LSR": {s: "L" + s for s in SUF}}
EPS = 1e-9


def reps(pat):
    out = {}
    for f in sorted(glob.glob(pat)):
        out[int(re.search(r"_r(\d+)\.npz$", os.path.basename(f)).group(1))] = f
    return out


def per_system(rec, s):
    """free reference, and per replica: bound mean, exposure %, open fraction, open mean."""
    fr = reps("/home/yuki/sasa/free_L%s_r*.npz" % s)
    fv = float(np.mean([np.load(fr[r])["res14"].mean() for r in sorted(fr)]))
    bd = reps("/home/yuki/sasa/bound_%s_r*.npz" % TAG[rec][s])
    rows, pooled = [], []
    for r in sorted(bd):
        a = np.asarray(np.load(bd[r])["res14"], dtype=float)
        rel = 100.0 * a / fv
        op = rel > EPS
        rows.append({"rep": r, "n": len(a), "bound_mean": float(a.mean()), "free": fv,
                     "expo": float(rel.mean()), "openfrac": 100.0 * float(op.mean()),
                     "open_mean": float(rel[op].mean()) if op.any() else 0.0,
                     "median": float(np.median(rel))})
        pooled.append(rel)
    return fv, rows, np.concatenate(pooled)


def se_contrast(a, b):
    """separation of two group means in standard errors of the difference"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    s = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return (b.mean() - a.mean()) / s


DATA, POOL = {}, {}
for rec in ("ILDR2", "LSR"):
    for s in SUF:
        fv, rows, pooled = per_system(rec, s)
        DATA[(rec, s)] = rows
        POOL[(rec, s)] = pooled

print("                 frames  median%%  open%%   open-state%%  time-avg%%")
for rec in ("ILDR2", "LSR"):
    for s in SUF:
        r = DATA[(rec, s)]
        p = POOL[(rec, s)]
        print("%-6s %-5s g%d  %5d   %6.2f  %5.1f   %5.1f        %5.1f"
              % (rec, LAB[s], GAM[s], len(p), np.median(p),
                 np.mean([x["openfrac"] for x in r]),
                 np.mean([x["open_mean"] for x in r]),
                 np.mean([x["expo"] for x in r])))

print()
for rec in ("ILDR2", "LSR"):
    free = [np.mean([x["openfrac"] for x in DATA[(rec, s)]]) for s in ("wt", "g")]
    imp = [np.mean([x["openfrac"] for x in DATA[(rec, s)]]) for s in ("c", "s", "t", "v")]
    print("%-6s open fraction: no-gamma %.1f%%  gamma %.1f%%  "
          "+%.2f SE  overlap %s"
          % (rec, np.mean(free), np.mean(imp), se_contrast(free, imp),
             "yes" if max(free) > min(imp) else "no"))
    e_free = [np.mean([x["expo"] for x in DATA[(rec, s)]]) for s in ("wt", "g")]
    e_imp = [np.mean([x["expo"] for x in DATA[(rec, s)]]) for s in ("c", "s", "t", "v")]
    print("       time average : no-gamma %.1f%%  gamma %.1f%%  %+.2f SE"
          % (np.mean(e_free), np.mean(e_imp), se_contrast(e_free, e_imp)))
    opn = [np.mean([x["openfrac"] for x in DATA[(rec, s)]]) for s in SUF]
    print("       buried fraction range %.0f-%.0f%% of frames"
          % (100 - max(opn), 100 - min(opn)))
    om = [np.mean([x["open_mean"] for x in DATA[(rec, s)]]) for s in SUF]
    print("       open-state exposure range %.0f-%.0f%%" % (min(om), max(om)))

np.savez("/home/yuki/res14_open.npz",
         **{"%s_%s" % (rec, s): np.array(
             [[x["rep"], x["bound_mean"], x["free"], x["expo"], x["openfrac"],
               x["open_mean"], x["median"]] for x in DATA[(rec, s)]])
            for rec in ("ILDR2", "LSR") for s in SUF},
         **{"pool_%s_%s" % (rec, s): POOL[(rec, s)]
            for rec in ("ILDR2", "LSR") for s in SUF})
print("\nwrote /home/yuki/res14_open.npz")
