"""Tyr99 across the six variants in both receptors.

The question is whether A14G raises Tyr99's contact, and whether the other variants move it
too. Occupancy cannot answer it: Tyr99 is in contact in 1.000 of frames in all twelve systems
with zero on-off switches, which is why the grouping sweep excluded it as degenerate. So
contact is counted - heavy-atom pairs within 4.5 A, and the number of distinct receptor
residues touched - which has range where a saturated indicator has none.

Replica is the unit, as everywhere else. Three things are reported for each quantity, by the
standard this project adopted: whether the replicate scatter is small enough to read the
ordering at all, the correlation with the gamma count, and the margin for A14G alone against
a null obtained by singling out each of the other five systems in turn.
"""
import os
import glob
import re
import numpy as np
from itertools import permutations

SRC = "/home/yuki/y99"
LSR = ["Lwt", "Lg", "Lc", "Ls", "Lt", "Lv"]
ILD = ["iwt", "ig", "ic", "is", "it", "iv"]
LAB = ["wt", "A14G", "A14C", "A14S", "A14T", "A14V"]
GAM = [0, 0, 1, 1, 2, 2]
KD_LSR = [0.00267, 0.00307, 0.0214, 0.0338, 0.0540, 0.0562]
KA_ILD = [2.87e5, 3.38e5, 1.11e4, 7445, 2721, 1344]
GUARD = 0.5


def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[s[j + 1]] == v[s[i]]:
            j += 1
        a = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[s[k]] = a
        i = j + 1
    return r


def sp(a, b):
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = sum((x - ma) ** 2 for x in ra) ** 0.5
    db = sum((y - mb) ** 2 for y in rb) ** 0.5
    return num / (da * db) if da and db else float("nan")


def pp(a, b):
    r0 = abs(sp(a, b))
    if r0 != r0:
        return float("nan")
    ge = sum(1 for p in permutations(range(6))
             if abs(sp([a[i] for i in p], b)) >= r0 - 1e-12)
    return ge / 720.0


def per_rep(tag, key):
    out = {}
    for f in sorted(glob.glob(os.path.join(SRC, "%s_r*.npz" % tag))):
        rep = int(re.search(r"_r(\d+)\.npz$", f).group(1))
        out[rep] = float(np.load(f)[key].mean())
    return [out[r] for r in sorted(out)]


def block(key, title, unit):
    print()
    print("=" * 94)
    print("%s  [%s]" % (title, unit))
    print("=" * 94)
    res = {}
    for recep, tags, rate, rname in (("LSR", LSR, KD_LSR, "kd"),
                                     ("ILDR2", ILD, KA_ILD, "ka")):
        m, sd, n = [], [], []
        for t in tags:
            v = per_rep(t, key)
            if not v:
                m.append(float("nan")); sd.append(float("nan")); n.append(0); continue
            m.append(float(np.mean(v)))
            sd.append(float(np.std(v, ddof=1)) if len(v) > 1 else float("nan"))
            n.append(len(v))
        span = max(m) - min(m)
        ratio = float(np.nanmean(sd)) / span if span else float("nan")
        res[recep] = (m, sd, ratio)
        print("  %-6s %s" % (recep, "  ".join("%-9s" % l for l in LAB)))
        print("         %s" % "  ".join("%-9.1f" % x for x in m))
        print("     ±SD %s" % "  ".join(
            ("%-9.1f" % x) if x == x else "%-9s" % "-" for x in sd))
        print("     n   %s" % "  ".join("%-9d" % x for x in n))
        print("     レプリカSD / 変異体間の幅 = %.2f  %s"
              % (ratio, "" if ratio <= GUARD else "<<< 0.5 超、順序は読めない"))
        print("     γ原子数との相関 rho = %+.3f (p %.4f)   実測 %s との相関 rho = %+.3f (p %.4f)"
              % (sp(m, GAM), pp(m, GAM), rname, sp(m, rate), pp(m, rate)))
        # A14G alone, against the null of singling out each other system
        gaps = []
        for i in range(6):
            others = [m[j] for j in range(6) if j != i]
            pooled = float(np.nanmean([s / np.sqrt(max(k, 1))
                                       for s, k in zip(sd, n)]))
            gaps.append((m[i] - float(np.mean(others))) / pooled if pooled else float("nan"))
        rk = sorted(range(6), key=lambda i: -abs(gaps[i])).index(1) + 1
        print("     A14G 単独の余裕 %+.2f SE   他5系で同じ計算をしたときの順位 %d / 6"
              % (gaps[1], rk))
        print("        " + "  ".join("%s %+.2f" % (LAB[i], gaps[i]) for i in range(6)))
    return res


print("Tyr99 — 複合体MD、重原子 4.5 A 以内、レプリカ単位")
print("（参考: 接触占有率は12系×3レプリカすべてで 1.000、スイッチ 0 回。飽和している）")
block("npair", "Tyr99 と受容体の接触ペア数", "重原子対")
block("nres", "Tyr99 が接触する受容体残基の数", "残基")
block("dmin", "Tyr99 の受容体への最接近距離", "A")
