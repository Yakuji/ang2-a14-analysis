"""Contact occupancy per Ang-2 residue at three replicas.

Descriptive rather than a correlation hunt: which contacts exist, how permanent they are, and
whether the set differs between the two receptors or between the variants. Every distance
result in this project used a mean distance, which cannot distinguish a permanent contact from
an intermittent one.

  occ   fraction of frames with a heavy atom within 4.5 A of the receptor
  nsw   formed/broken transitions, a direct measure of intermittency
"""
import os
import glob
import re
from itertools import permutations
import numpy as np
from scipy.stats import spearmanr, ttest_ind

D = "/home/yuki/occ3"
KA = {"Lwt": 1300000, "Lg": 2000000, "Lc": 1700000, "Ls": 2200000, "Lt": 2200000,
      "Lv": 1800000, "iwt": 450000, "ig": 400000, "ic": 15000, "is": 8000,
      "it": 1700, "iv": 1400}
KD = {"Lwt": 0.0025, "Lg": 0.003, "Lc": 0.016, "Ls": 0.02, "Lt": 0.05, "Lv": 0.04,
      "iwt": 0.013, "ig": 0.0028, "ic": 0.017, "is": 0.011, "it": 0.02, "iv": 0.003}
G = {"ILDR2": ["iwt", "ig", "ic", "is", "it", "iv"],
     "LSR": ["Lwt", "Lg", "Lc", "Ls", "Lt", "Lv"]}
LAB = {"wt": "wt", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}


def nm(t):
    return LAB.get(t[1:], t)


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


def reps(tag):
    return {int(re.search(r"_r(\d+)\.npz$", f).group(1)): np.load(f)
            for f in sorted(glob.glob(os.path.join(D, tag + "_r*.npz")))}


data = {t: reps(t) for g in G.values() for t in g}
ids = data["iwt"][sorted(data["iwt"])[0]]["resids"]
# resnames was written as an object array by occ3.py; these are files written in this
# session by that script, so loading them with pickle enabled is reading our own output.
_f = sorted(glob.glob(os.path.join(D, "iwt_r*.npz")))[0]
rn = np.array([str(x) for x in np.load(_f, allow_pickle=True)["resnames"]])

print("=" * 96)
print("1. 接触占有率が 0.2 を超える Ang-2 残基（野生型、3レプリカ平均）")
print("=" * 96)
print("  %-6s %-6s %-18s %-18s %s" % ("残基", "名前", "ILDR2 占有率", "LSR 占有率", "両受容体"))
keep = []
for i, r in enumerate(ids):
    oi = np.array([data["iwt"][k]["occ"][i] for k in sorted(data["iwt"])])
    ol = np.array([data["Lwt"][k]["occ"][i] for k in sorted(data["Lwt"])])
    if max(oi.mean(), ol.mean()) < 0.2:
        continue
    keep.append(i)
    both = "両方" if min(oi.mean(), ol.mean()) >= 0.2 else (
        "ILDR2のみ" if oi.mean() >= 0.2 else "LSRのみ")
    print("  %-6d %-6s %.3f +- %-11.3f %.3f +- %-11.3f %s"
          % (r, rn[i], oi.mean(), oi.std(ddof=1), ol.mean(), ol.std(ddof=1), both))
print("  計 %d 残基" % len(keep))

print()
print("=" * 96)
print("2. 変異体と野生型の占有率の差（Holm 補正、受容体ごと 5変異体 x %d残基）" % len(keep))
print("=" * 96)
for recep, tags in G.items():
    ps, cells = [], []
    for t in tags[1:]:
        for i in keep:
            a = np.array([data[t][k]["occ"][i] for k in sorted(data[t])])
            w = np.array([data[tags[0]][k]["occ"][i] for k in sorted(data[tags[0]])])
            dd = a.mean() - w.mean()
            se = np.sqrt(a.var(ddof=1) / len(a) + w.var(ddof=1) / len(w))
            ps.append(ttest_ind(a, w, equal_var=False).pvalue)
            cells.append((nm(t), int(ids[i]), dd, dd / se if se > 0 else np.nan))
    adj = holm(ps)
    hit = [(c, x) for c, x in zip(cells, adj) if abs(c[3]) >= 3]
    print("  %-6s 3 SE 超 %d / %d、Holm p<0.05 %d 件"
          % (recep, len(hit), len(cells), int((adj < 0.05).sum())))
    for (v, r, dd, sg), x in sorted(hit, key=lambda z: -abs(z[0][3]))[:12]:
        print("     %-6s 残基%4d  %+7.3f  %+6.1f SE  Holm p=%.3f" % (v, r, dd, sg, x))

print()
print("=" * 96)
print("3. 占有率 対 速度定数（6系、順列下限 2/720 = 0.0028）")
print("=" * 96)
for recep, tags in G.items():
    for rate, RT in (("ka", KA), ("kd", KD)):
        y = np.array([RT[t] for t in tags], float)
        ps, rows, const = [], [], []
        for i in keep:
            m = np.array([np.mean([data[t][k]["occ"][i] for k in sorted(data[t])])
                          for t in tags])
            sd = np.array([np.std([data[t][k]["occ"][i] for k in sorted(data[t])], ddof=1)
                           for t in tags])
            # a residue whose occupancy is the same in all six systems has no rank order to
            # correlate; spearmanr returns nan there and must not be counted as a result
            if m.max() - m.min() < 1e-9:
                const.append(int(ids[i]))
                continue
            rows.append((int(ids[i]), spearmanr(m, y).statistic, sd.mean(),
                         m.max() - m.min()))
            ps.append(exact_p(m, y))
        if not rows:
            print("  %-6s 対 %-3s: 変動する残基なし" % (recep, rate))
            continue
        adj = holm(ps)
        b = int(np.argmin(ps))
        r, rho, sdm, span = rows[b]
        print("  %-6s 対 %-3s: 最良 残基%d rho=%+.2f p=%.4f Holm=%.3f  "
              "レプリカSD %.3f 対 幅 %.3f 比 %.2f %s"
              % (recep, rate, r, rho, ps[b], adj[b], sdm, span,
                 sdm / span if span > 0 else np.inf,
                 "分離" if sdm < 0.5 * span else "ばらつきが覆う"))
        print("          検定 %d 残基、Holm有意 %d 件、全系同一で除外 %d 残基 %s"
              % (len(rows), int((adj < 0.05).sum()), len(const), const))

print()
print("=" * 96)
print("4. 断続性: 形成/破断の遷移回数（野生型、1001フレーム中）")
print("=" * 96)
print("  %-6s %-6s %-14s %-14s" % ("残基", "名前", "ILDR2 遷移", "LSR 遷移"))
for i in keep:
    ni = np.mean([data["iwt"][k]["nsw"][i] for k in sorted(data["iwt"])])
    nl = np.mean([data["Lwt"][k]["nsw"][i] for k in sorted(data["Lwt"])])
    print("  %-6d %-6s %-14.1f %-14.1f" % (ids[i], rn[i], ni, nl))
