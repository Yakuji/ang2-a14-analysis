"""Statistics for the SASA measurement: is position 14's side chain in the association path?

The steric reading of the gamma-atom result requires the side chain at position 14 to be
solvent exposed in the free state. A side chain packed against the ligand's own body cannot
obstruct an approaching receptor.

Three questions, in order of how much they bear on the hypothesis:

(1) Is the side chain exposed at all, and does its free-state exposure rise with gamma-atom
    count? A buried side chain falsifies the steric reading outright.
(2) Does free-state side-chain SASA order with ILDR2's ka, the way the gamma count does, and
    fail against LSR's flat ka?
(3) How much of position 14 is buried on binding, and does the buried area differ between
    receptors or between variants?

Replica as the unit, the replicate spread printed against the variant span, and LSR kept as
the internal control throughout.
"""
import os
import glob
import re
from itertools import permutations
import numpy as np
from scipy.stats import spearmanr, ttest_ind

D = "/home/yuki/sasa"
GAM = {"wt": 0, "g": 0, "c": 1, "s": 1, "t": 2, "v": 2}
RES = {"wt": "ALA", "g": "GLY", "c": "CYS", "s": "SER", "t": "THR", "v": "VAL"}
LAB = {"wt": "wt", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}
FREE = ["wt", "g", "c", "s", "t", "v"]
KA = {"Lwt": 1300000, "Lg": 2000000, "Lc": 1700000, "Ls": 2200000, "Lt": 2200000,
      "Lv": 1800000, "iwt": 450000, "ig": 400000, "ic": 15000, "is": 8000,
      "it": 1700, "iv": 1400}
KD = {"Lwt": 0.0025, "Lg": 0.003, "Lc": 0.016, "Ls": 0.02, "Lt": 0.05, "Lv": 0.04,
      "iwt": 0.013, "ig": 0.0028, "ic": 0.017, "is": 0.011, "it": 0.02, "iv": 0.003}
G = {"ILDR2": ["iwt", "ig", "ic", "is", "it", "iv"],
     "LSR": ["Lwt", "Lg", "Lc", "Ls", "Lt", "Lv"]}
SUF = {"iwt": "wt", "ig": "g", "ic": "c", "is": "s", "it": "t", "iv": "v",
       "Lwt": "wt", "Lg": "g", "Lc": "c", "Ls": "s", "Lt": "t", "Lv": "v"}


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


def load(which, tag):
    out = {}
    for f in sorted(glob.glob(os.path.join(D, "%s_%s_r*.npz" % (which, tag)))):
        out[int(re.search(r"_r(\d+)\.npz$", f).group(1))] = np.load(f)
    return out


def agg(which, tag, key):
    d = load(which, tag)
    return np.array([float(d[r][key].mean()) for r in sorted(d)])


print("=" * 96)
print("1. 遊離状態での位置14の溶媒露出（nm^2、3レプリカ）")
print("=" * 96)
print("  %-6s %-5s %-7s %-22s %-22s %s"
      % ("系", "残基", "gamma", "側鎖 SASA", "残基全体 SASA", "側鎖の露出率"))
sc_free, r14_free = {}, {}
for s in FREE:
    sc = agg("free", "L" + s, "sc14")
    r14 = agg("free", "L" + s, "res14")
    sc_free[s], r14_free[s] = sc, r14
    print("  %-6s %-5s %-7d %.4f +- %-13.4f %.4f +- %-13.4f %s"
          % (LAB[s], RES[s], GAM[s], sc.mean(), sc.std(ddof=1),
             r14.mean(), r14.std(ddof=1),
             "%.2f" % (sc.mean() / r14.mean()) if r14.mean() > 0 else "-"))

x = np.array([GAM[s] for s in FREE], float)
m = np.array([sc_free[s].mean() for s in FREE])
sd = np.array([sc_free[s].std(ddof=1) for s in FREE])
print()
print("  gamma原子数 対 側鎖SASA : rho=%+.3f  完全順列 p=%.4f"
      % (spearmanr(x, m).statistic, exact_p(x, m)))
print("  レプリカSD平均 %.4f  対  変異体間の幅 %.4f  比 %.2f  %s"
      % (sd.mean(), m.max() - m.min(),
         sd.mean() / (m.max() - m.min()) if m.max() > m.min() else np.inf,
         "分離できている" if sd.mean() < 0.5 * (m.max() - m.min()) else "ばらつきが覆う"))
nz = [s for s in FREE if sc_free[s].mean() > 0]
if nz:
    print("  側鎖が露出しているか: 最小 %s %.4f nm^2、最大 %s %.4f nm^2"
          % (LAB[min(nz, key=lambda s: sc_free[s].mean())],
             min(sc_free[s].mean() for s in nz),
             LAB[max(nz, key=lambda s: sc_free[s].mean())],
             max(sc_free[s].mean() for s in nz)))

print()
print("=" * 96)
print("2. 遊離状態の SASA 対 速度定数（6系、順列下限 2/720 = 0.0028）")
print("=" * 96)
for recep, tags in G.items():
    for key, name in (("sc14", "位置14 側鎖"), ("res14", "位置14 残基全体"),
                      ("iface", "界面18残基")):
        mm = np.array([agg("free", "L" + SUF[t], key).mean() for t in tags])
        ss = np.array([agg("free", "L" + SUF[t], key).std(ddof=1) for t in tags])
        span = mm.max() - mm.min()
        for rate, RT in (("ka", KA), ("kd", KD)):
            y = np.array([RT[t] for t in tags], float)
            print("  %-6s %-14s 対 %-3s: rho=%+.3f p=%.4f   SD %.4f 対 幅 %.4f 比 %.2f %s"
                  % (recep, name, rate, spearmanr(mm, y).statistic, exact_p(mm, y),
                     ss.mean(), span, ss.mean() / span if span > 0 else np.inf,
                     "分離" if ss.mean() < 0.5 * span else "覆う"))

print()
print("=" * 96)
print("3. 結合で埋まる面積（遊離 - 複合体、nm^2）")
print("=" * 96)
print("  %-6s %-6s %-20s %-20s %-20s"
      % ("受容体", "系", "位置14 側鎖 埋没", "位置14 残基 埋没", "界面18残基 埋没"))
bur = {}
for recep, tags in G.items():
    for t in tags:
        s = SUF[t]
        row = {}
        for key in ("sc14", "res14", "iface"):
            f = agg("free", "L" + s, key).mean()
            b = agg("bound", t, key)
            row[key] = f - b
        bur[t] = row
        print("  %-6s %-6s %+8.4f +- %-9.4f %+8.4f +- %-9.4f %+8.4f +- %-9.4f"
              % (recep, LAB[s],
                 row["sc14"].mean(), row["sc14"].std(ddof=1),
                 row["res14"].mean(), row["res14"].std(ddof=1),
                 row["iface"].mean(), row["iface"].std(ddof=1)))

print()
print("  埋没面積 対 速度定数、および野生型との差（Holm 補正、受容体ごと 5変異体 x 3量）")
for recep, tags in G.items():
    ps, cells = [], []
    for key, name in (("sc14", "側鎖"), ("res14", "残基"), ("iface", "界面")):
        for t in tags[1:]:
            a, w = bur[t][key], bur[tags[0]][key]
            dd = a.mean() - w.mean()
            se = np.sqrt(a.var(ddof=1) / len(a) + w.var(ddof=1) / len(w))
            ps.append(ttest_ind(a, w, equal_var=False).pvalue)
            cells.append((name, LAB[SUF[t]], dd, dd / se if se > 0 else np.nan))
    adj = holm(ps)
    hit = [(c, v) for c, v in zip(cells, adj) if abs(c[3]) >= 3]
    print("  %-6s 3 SE 超 %d / %d、Holm p<0.05 %d 件"
          % (recep, len(hit), len(cells), int((adj < 0.05).sum())))
    for (nmk, v, dd, sg), q in sorted(hit, key=lambda z: -abs(z[0][3]))[:8]:
        print("     %-6s %-6s %+8.4f  %+6.1f SE  Holm p=%.3f" % (nmk, v, dd, sg, q))
    for key, name in (("sc14", "側鎖"), ("res14", "残基"), ("iface", "界面")):
        mm = np.array([bur[t][key].mean() for t in tags])
        ss = np.array([bur[t][key].std(ddof=1) for t in tags])
        span = mm.max() - mm.min()
        for rate, RT in (("ka", KA), ("kd", KD)):
            y = np.array([RT[t] for t in tags], float)
            print("     %-6s 対 %-3s: rho=%+.3f p=%.4f  SD %.4f 対 幅 %.4f 比 %.2f %s"
                  % (name, rate, spearmanr(mm, y).statistic, exact_p(mm, y),
                     ss.mean(), span, ss.mean() / span if span > 0 else np.inf,
                     "分離" if ss.mean() < 0.5 * span else "覆う"))
