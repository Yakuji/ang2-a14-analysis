"""Source data for the figures, laid out for Prism.

One sheet per panel. Every sheet carries the individual determinations where they exist, so
Prism can draw the points rather than only the summary.
"""
import os
import glob
import re
import numpy as np
from scipy.stats import spearmanr, linregress
import openpyxl
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

OUT = "/mnt/c/Users/y-niw/OneDrive/Desktop/Ang2_MD/FIGURE_DATA_2026-10-24.xlsx"
SUF = ["wt", "g", "c", "s", "t", "v"]
LAB = {"wt": "WT", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}
RES = {"wt": "Ala", "g": "Gly", "c": "Cys", "s": "Ser", "t": "Thr", "v": "Val"}
GAM = {"wt": 0, "g": 0, "c": 1, "s": 1, "t": 2, "v": 2}
VOL = {"wt": 27.5, "g": 0.0, "c": 44.6, "s": 29.3, "t": 51.3, "v": 71.5}
ILD = ["iwt", "ig", "ic", "is", "it", "iv"]
LSRT = ["Lwt", "Lg", "Lc", "Ls", "Lt", "Lv"]
RT = 8.314462618e-3 * 298.15
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

wb = openpyxl.Workbook()
wb.remove(wb.active)


def sheet(name, header, rows, note=None):
    ws = wb.create_sheet(name[:31])
    if note:
        ws.append([note])
        ws["A1"].font = Font(italic=True, size=9)
        ws.append([])
    ws.append(header)
    for c in ws[ws.max_row]:
        c.font = Font(bold=True)
    for r in rows:
        ws.append(list(r))
    for i, h in enumerate(header, 1):
        ws.column_dimensions[get_column_letter(i)].width = max(11, min(30, len(str(h)) + 3))
    return ws


# ---------------------------------------------------------------- 1 raw kinetics
rows = []
for rec in ("ILDR2", "LSR"):
    for s in SUF:
        for i in range(3):
            ka, kd = KA[rec][s][i], KD[rec][s][i]
            rows.append([rec, LAB[s], RES[s], GAM[s], i + 1, ka, kd, kd / ka * 1e9])
sheet("1_raw_kinetics",
      ["Receptor", "Variant", "Residue", "gamma atoms", "Replicate",
       "ka (1/Ms)", "kd (1/s)", "KD (nM)"], rows,
      "Three independent SPR determinations per system. KD computed per replicate as kd/ka.")

# ---------------------------------------------------------------- 2 summary
rows = []
for rec in ("ILDR2", "LSR"):
    for s in SUF:
        a, d = np.log(KA[rec][s]), np.log(KD[rec][s])
        k = a - d
        f = lambda v: (np.exp(v.mean()), np.exp(v.mean() - 1.96 * v.std(ddof=1) / 3 ** .5),
                       np.exp(v.mean() + 1.96 * v.std(ddof=1) / 3 ** .5))
        ga, la, ha = f(a)
        gd, ld, hd = f(d)
        gk, lk, hk = f(-k)
        rows.append([rec, LAB[s], GAM[s], ga, la, ha, gd, ld, hd,
                     gk * 1e9, lk * 1e9, hk * 1e9])
sheet("2_summary_geomean",
      ["Receptor", "Variant", "gamma atoms",
       "ka geomean", "ka 95% lo", "ka 95% hi",
       "kd geomean", "kd 95% lo", "kd 95% hi",
       "KD geomean (nM)", "KD 95% lo", "KD 95% hi"], rows,
      "Geometric means with 95% intervals, estimated on ln k (rate constants are log-normal).")

# ---------------------------------------------------------------- 3 grouped by gamma
rows = []
for rec, const, D in (("ILDR2", "ka", KA["ILDR2"]), ("LSR", "kd", KD["LSR"]),
                      ("ILDR2", "kd", KD["ILDR2"]), ("LSR", "ka", KA["LSR"])):
    for s in SUF:
        for i in range(3):
            rows.append([rec, const, GAM[s], LAB[s], i + 1, D[s][i],
                         "determining" if (rec, const) in (("ILDR2", "ka"), ("LSR", "kd"))
                         else "control"])
sheet("3_by_gamma_count",
      ["Receptor", "Constant", "gamma atoms", "Variant", "Replicate", "Value", "Role"], rows,
      "For the grouped scatter/bar: x = gamma-atom count, y = rate constant, log axis.")

# ---------------------------------------------------------------- 4 Phi
rows = []
for rec in ("ILDR2", "LSR"):
    M = {}
    for s in SUF:
        a = np.log(KA[rec][s])
        k = a - np.log(KD[rec][s])
        M[s] = (a.mean(), a.std(ddof=1) / 3 ** .5, k.mean(), k.std(ddof=1) / 3 ** .5)
    w = M["wt"]
    rng = np.random.default_rng(7)
    for s in SUF[1:]:
        m = M[s]
        A, B = m[0] - w[0], m[2] - w[2]
        aa = rng.normal(m[0], m[1], 400000) - rng.normal(w[0], w[1], 400000)
        bb = rng.normal(m[2], m[3], 400000) - rng.normal(w[2], w[3], 400000)
        ok = np.abs(bb) > 0.5
        lo, hi = (np.percentile(aa[ok] / bb[ok], [16, 84]) if ok.mean() > 0.6
                  else (np.nan, np.nan))
        rows.append([rec, LAB[s], GAM[s], -RT * A, RT * aa.std(), -RT * B, RT * bb.std(),
                     A / B, lo, hi])
sheet("4_phi_values",
      ["Receptor", "Variant", "gamma atoms", "ddG_TS (kJ/mol)", "SE", "ddG_eq (kJ/mol)",
       "SE", "Phi", "Phi 68% lo", "Phi 68% hi"], rows,
      "Phi = ddG_TS / ddG_eq. Errors by Monte Carlo on the measured log-means; "
      "KA is paired within each replicate.")

# ---------------------------------------------------------------- 5 Broensted plot
rows = []
for rec in ("ILDR2", "LSR"):
    M = {}
    for s in SUF:
        a = np.log(KA[rec][s])
        M[s] = (a.mean(), (a - np.log(KD[rec][s])).mean())
    w = M["wt"]
    for s in SUF[1:]:
        rows.append([rec, LAB[s], GAM[s],
                     -RT * (M[s][1] - w[1]), -RT * (M[s][0] - w[0])])
sheet("5_broensted_plot",
      ["Receptor", "Variant", "gamma atoms", "ddG_eq (x axis)", "ddG_TS (y axis)"], rows,
      "Broensted plot: slope = Phi. A slope of 1 means the interaction is fully formed at "
      "the transition state; 0 means it is not formed there.")

# ---------------------------------------------------------------- 6 predictor regression
rows = []
for rec, const, D in (("ILDR2", "ka", KA["ILDR2"]), ("LSR", "kd", KD["LSR"])):
    y = np.array([np.log(D[s]).mean() for s in SUF])
    for nm, P in (("gamma atoms", GAM), ("side-chain volume", VOL)):
        x = np.array([P[s] for s in SUF], float)
        r = linregress(x, y)
        for s in SUF:
            rows.append([rec, const, nm, LAB[s], P[s], np.log(D[s]).mean(),
                         np.exp(np.log(D[s]).mean()), r.rvalue ** 2])
sheet("6_predictor_regression",
      ["Receptor", "Constant", "Predictor", "Variant", "Predictor value",
       "ln k", "k (geomean)", "R2 of this predictor"], rows,
      "Regression of ln k on each side-chain predictor. Gamma count R2 0.940 and 0.908; "
      "volume 0.752 and 0.640.")


def reps(p):
    return {int(re.search(r"_r(\d+)\.npz$", f).group(1)): f for f in sorted(glob.glob(p))}


# ---------------------------------------------------------------- 7 MD residue 14 exposure
rows = []
for rec, tags in (("ILDR2", ILD), ("LSR", LSRT)):
    for t, s in zip(tags, SUF):
        fr = reps("/home/yuki/sasa/free_L%s_r*.npz" % s)
        fv = np.mean([float(np.load(fr[r])["res14"].mean()) for r in sorted(fr)])
        bd = reps("/home/yuki/sasa/bound_%s_r*.npz" % t)
        for j, r in enumerate(sorted(bd)):
            b = float(np.load(bd[r])["res14"].mean())
            rows.append([rec, LAB[s], GAM[s], j + 1, b, fv, 100 * b / fv])
sheet("7_MD_res14_exposure",
      ["Receptor", "Variant", "gamma atoms", "MD replica",
       "bound SASA (nm2)", "free SASA (nm2)", "remaining exposure (%)"], rows,
      "Residue 14 in the complex. Remaining exposure = bound SASA / that ligand's own free "
      "SASA, which removes the side-chain size difference.")

# ---------------------------------------------------------------- 8 MD Fmax
FP = np.load("/home/yuki/force_prof.npz")
rows = []
for rec, tags in (("LSR", LSRT), ("ILDR2", ILD)):
    for t, s in zip(tags, SUF):
        acc = {}
        for k in FP.files:
            if k.startswith(t + "__"):
                acc.setdefault(int(re.search(r"^r(\d+)_", k[len(t) + 2:]).group(1)),
                               []).append(float(FP[k][0]))
        for j, r in enumerate(sorted(acc)):
            rows.append([rec, LAB[s], GAM[s], j + 1, float(np.mean(acc[r]))])
sheet("8_MD_Fmax",
      ["Receptor", "Variant", "gamma atoms", "MD replica", "Fmax (kJ/mol/nm)"], rows,
      "Peak pulling force, averaged over the six starting conformations of each MD replica. "
      "LSR Fmax against measured kd: rho = -0.943, p = 0.0167.")

# ---------------------------------------------------------------- 9 negative control
rows = []
for rec, tags in (("LSR", LSRT), ("ILDR2", ILD)):
    f0 = reps("/home/yuki/recv/%s_r*.npz" % tags[0])
    d0 = np.load(f0[sorted(f0)[0]])
    ca, ids = d0["ca_resids"], d0["resids"]
    occ = np.array([np.load(x)["occ"] for x in f0.values()]).mean(0)
    dmn = np.array([np.load(x)["dmin"] for x in f0.values()]).mean(0)
    for const, D in (("ka", KA[rec]), ("kd", KD[rec])):
        y = np.array([np.exp(np.log(D[s]).mean()) for s in SUF])
        n = 0
        for i in range(len(ca)):
            v = np.array([np.mean([float(np.load(reps("/home/yuki/recv/%s_r*.npz" % t)[r])
                                         ["rmsf"][i])
                                   for r in sorted(reps("/home/yuki/recv/%s_r*.npz" % t))])
                          for t in tags])
            rho = spearmanr(v, y).statistic
            if abs(rho) >= 0.9999:
                j = np.where(ids == int(ca[i]))[0]
                rows.append([rec, const, int(ca[i]), rho,
                             float(dmn[j[0]]) / 10 if len(j) else np.nan,
                             float(occ[j[0]]) if len(j) else np.nan])
                n += 1
        rows.append([rec, const, "TOTAL", n, len(ca), round(len(ca) * 2 / 720, 2)])
sheet("9_negative_control",
      ["Receptor", "Constant", "Receptor residue", "rho (or count)",
       "distance to Ang-2 (nm) / n residues", "contact occupancy / expected by chance"], rows,
      "Every receptor alpha carbon scanned against both rate constants. The TOTAL rows give "
      "the number of perfect correlations found, the number of residues scanned, and the "
      "number expected by chance at the 2/720 permutation floor.")

wb.save(OUT)
print("written:", OUT)
for ws in wb.worksheets:
    print("  %-24s %d rows" % (ws.title, ws.max_row))
