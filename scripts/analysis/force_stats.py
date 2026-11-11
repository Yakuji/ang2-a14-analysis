"""Force profile measures against the kinetics, replica as the unit.

Five measures per pull, averaged over a replica's six conformations, then the five variant
replica means tested against ka and kd by Spearman with an exact permutation p. Holm across
the five measures within each receptor and each rate constant.

The pull is far from equilibrium, so Fmax and work are not free energies. They are compared
only between systems run under one identical protocol, which is what the design provides.
"""
import numpy as np
import re
from itertools import permutations
from scipy.stats import spearmanr, ttest_ind

Z = np.load("/home/yuki/force_prof.npz")
LAB = ["Fmax", "work", "dx@max", "npeaks", "F_late", "dx_tot"]
KA = {"ig": 400000, "ic": 15000, "is": 8000, "it": 1700, "iv": 1400,
      "Lg": 2000000, "Lc": 1700000, "Ls": 2200000, "Lt": 2200000, "Lv": 1800000}
KD = {"ig": 0.0028, "ic": 0.017, "is": 0.011, "it": 0.02, "iv": 0.003,
      "Lg": 0.003, "Lc": 0.016, "Ls": 0.02, "Lt": 0.05, "Lv": 0.04}
G = {"ILDR2": ("iwt", ["ig", "ic", "is", "it", "iv"]),
     "LSR": ("Lwt", ["Lg", "Lc", "Ls", "Lt", "Lv"])}
NAME = {"g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}


def byrep(tag):
    d = {}
    for k in Z.files:
        t, conf = k.split("__")
        if t != tag:
            continue
        rep = int(re.match(r"r(\d+)_", conf).group(1))
        d.setdefault(rep, []).append(Z[k])
    return {r: np.array(v) for r, v in d.items()}


def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    a = np.empty(len(ps))
    run = 0.0
    for k, i in enumerate(o):
        run = max(run, (len(ps) - k) * ps[i])
        a[i] = min(1.0, run)
    return a


def exact_p(x, y):
    obs = abs(spearmanr(x, y).statistic)
    hit = tot = 0
    for p in permutations(range(len(y))):
        tot += 1
        if abs(spearmanr(x, [y[i] for i in p]).statistic) >= obs - 1e-12:
            hit += 1
    return hit / tot


for recep in ("ILDR2", "LSR"):
    wtag, tags = G[recep]
    print("=" * 92)
    print(recep)
    print("=" * 92)
    wt = byrep(wtag)
    wmean = np.array([wt[r].mean(0) for r in sorted(wt)])
    print("  系ごとのレプリカ平均（3レプリカの平均 ± レプリカ間SD）")
    hdr = "  %-6s" % "系"
    for L in LAB[:5]:
        hdr += "%16s" % L
    print(hdr)
    line = "  %-6s" % "wt"
    for j in range(5):
        line += "%9.1f±%-6.1f" % (wmean[:, j].mean(), wmean[:, j].std(ddof=1))
    print(line)
    M = {}
    for t in tags:
        b = byrep(t)
        bm = np.array([b[r].mean(0) for r in sorted(b)])
        M[t] = bm
        line = "  %-6s" % NAME[t[-1]]
        for j in range(5):
            line += "%9.1f±%-6.1f" % (bm[:, j].mean(), bm[:, j].std(ddof=1))
        print(line)
    print()
    for rate, RT in (("ka", KA), ("kd", KD)):
        ps, rs = [], []
        for j in range(5):
            m = np.array([M[t][:, j].mean() for t in tags])
            y = np.array([RT[t] for t in tags], float)
            rs.append(spearmanr(m, y).statistic)
            ps.append(exact_p(m, y))
        adj = holm(ps)
        print("  対 %s :" % rate, end="")
        for j in range(5):
            mk = " <-" if adj[j] < 0.05 else ""
            print("  %s rho=%+.2f p=%.3f Holm=%.3f%s" % (LAB[j], rs[j], ps[j], adj[j], mk), end="")
        print()
    print()
    print("  各変異体 対 野生型（SE倍、レプリカ間ばらつきから）")
    hdr = "  %-6s" % ""
    for L in LAB[:5]:
        hdr += "%12s" % L
    print(hdr)
    ps2 = []
    cl = []
    for t in tags:
        line = "  %-6s" % NAME[t[-1]]
        for j in range(5):
            a = M[t][:, j]
            w = wmean[:, j]
            se = np.sqrt(a.var(ddof=1) / len(a) + w.var(ddof=1) / len(w))
            sg = (a.mean() - w.mean()) / se if se > 0 else np.nan
            p = ttest_ind(a, w, equal_var=False).pvalue
            ps2.append(p)
            cl.append((NAME[t[-1]], LAB[j], a.mean() - w.mean(), sg))
            line += "%12s" % ("%+.1f" % sg)
        print(line)
    a2 = holm(ps2)
    print()
    print("  3 SE 以上:")
    any_h = False
    for (v, L, d, sg), x in zip(cl, a2):
        if abs(sg) >= 3:
            any_h = True
            print("    %-6s %-8s %+10.1f  %+5.1f SE  Holm p=%.3f" % (v, L, d, sg, x))
    if not any_h:
        print("    なし")
    print()
