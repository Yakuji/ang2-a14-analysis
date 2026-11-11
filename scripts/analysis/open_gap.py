"""Score the opening frequency with the project's own separation statistic.

The manuscript's "standard errors" are grouping.py's gap2: the gap between the nearest
members of the two groups, divided by the mean standard error of a system mean over its MD
replicas. A quick two-sample t-style contrast gives a different number (+5.22 rather than
+3.77 here) because it is driven by the spread between systems, which is the effect, not the
noise. The figure, the workbook and the text all have to use the one definition, so it is
taken from grouping.py here.
"""
import numpy as np

D = np.load("/home/yuki/res14_open.npz")
SUF = ["wt", "g", "c", "s", "t", "v"]
LAB = {"wt": "WT", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}
COL = {"expo": 3, "openfrac": 4, "open_mean": 5}


def gap2(per_sys):
    m = np.array([np.mean(v) for v in per_sys])
    sd = np.array([np.std(v, ddof=1) for v in per_sys])
    sem = sd / np.sqrt([len(v) for v in per_sys])
    semm = sem.mean()
    lo, hi = m[[0, 1]], m[[2, 3, 4, 5]]
    g = min(hi) - max(lo)
    if g < 0:
        g = min(lo) - max(hi)
    return g / semm, m, sd, sem, semm, (min(hi) > max(lo) or min(lo) > max(hi))


for key, name in (("openfrac", "open fraction (% of frames)"),
                  ("expo", "time-averaged exposure (% of free)"),
                  ("open_mean", "exposure when open (% of free)")):
    print("==", name)
    for rec in ("ILDR2", "LSR"):
        per = [D["%s_%s" % (rec, s)][:, COL[key]] for s in SUF]
        g, m, sd, sem, semm, nover = gap2(per)
        res = sd.mean() / (m.max() - m.min())
        print("   %-6s  no-gamma %5.1f  gamma %5.1f  gap %+5.2f SE  "
              "overlap %-3s  res_sd/span %.3f"
              % (rec, m[[0, 1]].mean(), m[[2, 3, 4, 5]].mean(), g,
                 "no" if nover else "yes", res))
        print("          means " + "  ".join("%s %5.1f+-%.1f" % (LAB[s], mi, si)
                                             for s, mi, si in zip(SUF, m, sem)))
