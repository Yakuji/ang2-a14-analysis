"""Everything kinetic, redone on the measured n = 3 data instead of figure readings.

Until now every rate constant in this project was a value read off a published figure, which
is why only their ORDER was ever tested and why the Phi values carried a hand-assumed +-20%.
These are three independent determinations per system per receptor, so the quantities can be
estimated with their own errors and the conclusions can be checked rather than ranked.

Rate constants are log-normally distributed, so means and errors are taken on ln k and
reported as geometric means with multiplicative confidence intervals. KA is formed per
replicate as ka/kd, which keeps the pairing within an experiment -- and is confirmed by the
reported KD column: ILDR2 wild type replicate 1 gives 0.0318/448000 = 71 nM against a
reported 71 nM.

ILDR2 A14T and A14V have KD entries censored at >10000 nM. They are not needed: ka and kd are
both reported for every replicate, so KA is computed directly.
"""
import numpy as np
from itertools import permutations
from scipy.stats import spearmanr

RT = 8.314e-3 * 298.15

KA = {
 "LSR": {"wt": [1270000, 1190000, 1950000], "g": [1410000, 1300000, 1570000],
         "c": [1700000, 1930000, 4810000], "s": [1470000, 2440000, 3820000],
         "t": [2180000,  925000, 4390000], "v": [1500000, 3150000, 3230000]},
 "ILDR2": {"wt": [448000, 305000, 173000], "g": [582000, 259000, 256000],
           "c": [10800, 3900, 32200], "s": [12100, 5430, 6280],
           "t": [2710, 2580, 2880], "v": [491, 2130, 2320]}}
KD = {
 "LSR": {"wt": [0.00234, 0.00258, 0.00316], "g": [0.00283, 0.00351, 0.00291],
         "c": [0.0106, 0.03, 0.0306], "s": [0.0214, 0.0329, 0.0547],
         "t": [0.0485, 0.0292, 0.111], "v": [0.0317, 0.0776, 0.0721]},
 "ILDR2": {"wt": [0.0318, 0.0156, 0.00948], "g": [0.00275, 0.00281, 0.00255],
           "c": [0.0276, 0.009, 0.0598], "s": [0.0287, 0.0105, 0.0114],
           "t": [0.0381, 0.00907, 0.0105], "v": [0.0198, 0.0075, 0.0081]}}

SUF = ["wt", "g", "c", "s", "t", "v"]
LAB = {"wt": "wt", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}
RES = {"wt": "ALA", "g": "GLY", "c": "CYS", "s": "SER", "t": "THR", "v": "VAL"}
GAM = {"wt": 0, "g": 0, "c": 1, "s": 1, "t": 2, "v": 2}
# the figure readings every earlier analysis used, for comparison
OLD_KA = {"LSR": dict(zip(SUF, [1300000., 2000000., 1700000., 2200000., 2200000., 1800000.])),
          "ILDR2": dict(zip(SUF, [450000., 400000., 15000., 8000., 1700., 1400.]))}
OLD_KD = {"LSR": dict(zip(SUF, [0.0025, 0.003, 0.016, 0.02, 0.05, 0.04])),
          "ILDR2": dict(zip(SUF, [0.013, 0.0028, 0.017, 0.011, 0.02, 0.003]))}
rng = np.random.default_rng(11)
N = 200000


def lstat(v):
    """mean and standard error of ln k, plus the geometric mean."""
    L = np.log(np.asarray(v, float))
    return L.mean(), L.std(ddof=1) / np.sqrt(len(L)), float(np.exp(L.mean()))


def exact_p(x, y):
    obs = abs(spearmanr(x, y).statistic)
    hit = tot = 0
    for p in permutations(range(len(y))):
        tot += 1
        if abs(spearmanr(x, [y[i] for i in p]).statistic) >= obs - 1e-12:
            hit += 1
    return hit / tot


print("=" * 100)
print("1. 実測 n=3 の幾何平均と 95% 信頼区間（対数正規として ln k 上で推定）")
print("=" * 100)
S = {}
for rec in ("ILDR2", "LSR"):
    print()
    print("  == %s" % rec)
    print("     %-6s %-5s %-3s %-28s %-26s %s"
          % ("系", "残基", "γ", "ka (1/Ms)", "kd (1/s)", "KD = kd/ka (nM)"))
    S[rec] = {}
    for s in SUF:
        ma, sa, ga = lstat(KA[rec][s])
        md, sd, gd = lstat(KD[rec][s])
        kdv = np.array(KD[rec][s]) / np.array(KA[rec][s]) * 1e9
        mk, sk, gk = lstat(kdv)
        S[rec][s] = dict(la=ma, sea=sa, ld=md, sed=sd, ga=ga, gd=gd, gk=gk)
        print("     %-6s %-5s %-3d %10.4g [%8.4g,%8.4g] %8.4g [%7.4g,%7.4g] %8.4g [%6.4g,%6.4g]"
              % (LAB[s], RES[s], GAM[s],
                 ga, np.exp(ma - 1.96 * sa), np.exp(ma + 1.96 * sa),
                 gd, np.exp(md - 1.96 * sd), np.exp(md + 1.96 * sd),
                 gk, np.exp(mk - 1.96 * sk), np.exp(mk + 1.96 * sk)))

print()
print("=" * 100)
print("2. 図読み値との比較（これまでの全解析が使っていた値）")
print("=" * 100)
for rec in ("ILDR2", "LSR"):
    print("  == %s" % rec)
    print("     %-6s %-22s %-22s %-20s %s"
          % ("系", "ka 図読み → 実測", "", "kd 図読み → 実測", ""))
    for s in SUF:
        print("     %-6s %10.4g → %-10.4g (x%.2f)  %8.4g → %-8.4g (x%.2f)"
              % (LAB[s], OLD_KA[rec][s], S[rec][s]["ga"],
                 S[rec][s]["ga"] / OLD_KA[rec][s],
                 OLD_KD[rec][s], S[rec][s]["gd"], S[rec][s]["gd"] / OLD_KD[rec][s]))
    print()

print("=" * 100)
print("3. gamma 原子数との関係（実測値で再検定）")
print("=" * 100)
x = np.array([GAM[s] for s in SUF], float)
for rec in ("ILDR2", "LSR"):
    print()
    print("  == %s" % rec)
    for nm, key in (("ka", "ga"), ("kd", "gd")):
        m = np.array([S[rec][s][key] for s in SUF])
        rho = spearmanr(x, m).statistic
        p = exact_p(x, m)
        g = {}
        for s in SUF:
            g.setdefault(GAM[s], []).append(S[rec][s][key])
        ks = sorted(g)
        mono = all(min(g[a]) > max(g[b]) for a, b in zip(ks, ks[1:]))
        anti = all(max(g[a]) < min(g[b]) for a, b in zip(ks, ks[1:]))
        sep = "単調・重なりなし" if (mono or anti) else "重なりあり"
        span = max(max(g[k]) for k in ks) / min(min(g[k]) for k in ks)
        print("     対 %-3s rho=%+.3f p=%.4f  %s  全幅 %.0f倍" % (nm, rho, p, sep, span))
        print("        群: " + "   ".join(
            "%d個 %s" % (k, " ".join("%.4g" % v for v in g[k])) for k in ks))
        if mono or anti:
            fac = [max(min(g[a]) / max(g[b]), max(g[b]) / min(g[a]))
                   for a, b in zip(ks, ks[1:])]
            print("        群間の比: " + " / ".join("%.1f倍" % f for f in fac))

print()
print("=" * 100)
print("4. Phi 値 — 実測の誤差を伝播（モンテカルロ %d 回）" % N)
print("=" * 100)
print("  Phi = ddG+ / ddG0 = [ln ka_m - ln ka_w] / [ln KA_m - ln KA_w],  KA = ka/kd")
for rec in ("ILDR2", "LSR"):
    print()
    print("  == %s" % rec)
    print("     %-6s %-3s %-18s %-18s %-24s %s"
          % ("系", "γ", "ddG+ (kJ/mol)", "ddG0 (kJ/mol)", "Phi [68%区間]", "判定"))
    w = S[rec]["wt"]
    for s in SUF[1:]:
        m = S[rec][s]
        # independent draws for each measured log-mean
        aw = rng.normal(w["la"], w["sea"], N)
        dw = rng.normal(w["ld"], w["sed"], N)
        am = rng.normal(m["la"], m["sea"], N)
        dm = rng.normal(m["ld"], m["sed"], N)
        A = am - aw
        B = (am - dm) - (aw - dw)
        dgp, dg0 = -RT * A, -RT * B
        ok = np.abs(B) > 0.5
        phi_pt = (m["la"] - w["la"]) / ((m["la"] - m["ld"]) - (w["la"] - w["ld"]))
        if ok.mean() < 0.6:
            v = "KA がほぼ動かず Phi 不定"
            q = (np.nan, np.nan)
        else:
            p = A[ok] / B[ok]
            q = np.percentile(p, [16, 84])
            v = ("遷移状態で形成済み" if phi_pt > 0.7 else
                 "遷移状態では未形成" if phi_pt < 0.3 else "部分的")
        print("     %-6s %-3d %+7.1f +- %-7.1f %+7.1f +- %-7.1f %7s [%6s,%6s]  %s"
              % (LAB[s], GAM[s], -RT * (m["la"] - w["la"]), RT * np.std(A),
                 -RT * ((m["la"] - m["ld"]) - (w["la"] - w["ld"])), RT * np.std(B),
                 "%.2f" % phi_pt,
                 "nan" if np.isnan(q[0]) else "%.2f" % q[0],
                 "nan" if np.isnan(q[1]) else "%.2f" % q[1], v))

print()
print("=" * 100)
print("5. gamma 原子 1 個あたりの障壁上昇（ILDR2）と複合体の不安定化（LSR）")
print("=" * 100)
for rec, key, nm in (("ILDR2", "la", "会合障壁 ddG+"), ("LSR", None, "複合体 ddG0")):
    w = S[rec]["wt"]
    print("  == %s  %s" % (rec, nm))
    for gg in (1, 2):
        ss = [s for s in SUF if GAM[s] == gg]
        if rec == "ILDR2":
            d = [-RT * (S[rec][s]["la"] - w["la"]) for s in ss]
        else:
            d = [-RT * ((S[rec][s]["la"] - S[rec][s]["ld"]) - (w["la"] - w["ld"]))
                 for s in ss]
        print("     γ%d個 (%s): %s kJ/mol   平均 %+.1f"
              % (gg, ", ".join(LAB[s] for s in ss),
                 " / ".join("%+.1f" % v for v in d), float(np.mean(d))))
