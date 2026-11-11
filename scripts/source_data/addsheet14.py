# -*- coding: utf-8 -*-
"""Sheet 14: the receptor-side false-positive control, on the measured n = 3 kinetics.

Sheet 9 reports this control as four totals and two hits. For a figure the distribution
itself is wanted - every receptor alpha carbon's rank correlation against each rate
constant, with how far that residue sits from Ang-2 and whether it touches it - because the
point of the control is not that two hits exist but that the hits are remote residues with
zero contact, and that two is what 336 scans at the 2/720 floor produce anyway.

It also re-scores the control on the three-determination kinetics; sheet 9 was built on the
single figure-read constants those replaced. The outcome is unchanged: 2 hits against 0.93
expected, both on LSR kd, both with occupancy 0.000.
"""
import os
import re
import glob
import numpy as np
import openpyxl
from openpyxl.styles import Font
from scipy.stats import spearmanr

XL = "/mnt/c/Users/y-niw/OneDrive/Desktop/Ang2_MD/FIGURE_DATA_2026-10-24.xlsx"
D = "/home/yuki/recv"
SUF = ["wt", "g", "c", "s", "t", "v"]
KA = {"LSR": {"wt": [1270000, 1190000, 1950000], "g": [1410000, 1300000, 1570000],
              "c": [1700000, 1930000, 4810000], "s": [1470000, 2440000, 3820000],
              "t": [2180000, 925000, 4390000], "v": [1500000, 3150000, 3230000]},
      "ILDR2": {"wt": [448000, 305000, 173000], "g": [582000, 259000, 256000],
                "c": [10800, 3900, 32200], "s": [12100, 5430, 6280],
                "t": [2710, 2580, 2880], "v": [491, 2130, 2320]}}
KD = {"LSR": {"wt": [0.00234, 0.00258, 0.00316], "g": [0.00283, 0.00351, 0.00291],
              "c": [0.0106, 0.03, 0.0306], "s": [0.0214, 0.0329, 0.0547],
              "t": [0.0485, 0.0292, 0.111], "v": [0.0317, 0.0776, 0.0721]},
      "ILDR2": {"wt": [0.0318, 0.0156, 0.00948], "g": [0.00275, 0.00281, 0.00255],
                "c": [0.0276, 0.009, 0.0598], "s": [0.0287, 0.0105, 0.0114],
                "t": [0.0381, 0.00907, 0.0105], "v": [0.0198, 0.0075, 0.0081]}}
FLOOR = 2.0 / 720.0


def load(tag):
    return {int(re.search(r"_r(\d+)\.npz$", f).group(1)): np.load(f)
            for f in sorted(glob.glob(os.path.join(D, tag + "_r*.npz")))}


rows = []
for rec, pre in (("ILDR2", "i"), ("LSR", "L")):
    sysd = {s: load(pre + s) for s in SUF}
    wt = sysd["wt"]
    resids = np.asarray(wt[sorted(wt)[0]]["ca_resids"])
    n = len(resids)
    rmsf = np.array([[np.mean([sysd[s][r]["rmsf"][k] for r in sorted(sysd[s])])
                      for s in SUF] for k in range(n)])
    dmin = np.array([np.mean([wt[r]["dmin"][k] for r in sorted(wt)]) for k in range(n)])
    occ = np.array([np.mean([wt[r]["occ"][k] for r in sorted(wt)]) for k in range(n)])
    for const, tab in (("ka", KA[rec]), ("kd", KD[rec])):
        y = [float(np.exp(np.log(tab[s]).mean())) for s in SUF]
        for k in range(n):
            if not np.all(np.isfinite(rmsf[k])) or np.ptp(rmsf[k]) == 0:
                continue
            rows.append([rec, const, int(resids[k]),
                         round(float(spearmanr(rmsf[k], y).statistic), 4),
                         round(float(dmin[k]) / 10.0, 3),
                         round(float(occ[k]), 4)])

w = openpyxl.load_workbook(XL)
if "14_receptor_scan" in w.sheetnames:
    del w["14_receptor_scan"]
s = w.create_sheet("14_receptor_scan")
s.cell(row=1, column=1).value = (
    "False-positive control. Every receptor alpha carbon's RMSF is correlated across the six "
    "systems with each measured rate constant. With six systems the exact permutation floor "
    "is 2/720, so 336 scans per constant produce about one perfect correlation by chance "
    "whatever the biology.")
s.cell(row=2, column=1).value = (
    "Two turn up, both against LSR kd: receptor residues 266 and 281, rho = -1.000, sitting "
    "1.68 and 3.25 nm from Ang-2 with contact occupancy 0.000. They never touch the ligand. "
    "Plot rho against distance, or against occupancy, to show that the perfect correlations "
    "live where nothing can be transmitted.")
s.cell(row=3, column=1).value = (
    "Re-scored here on the three-determination kinetics; sheet 9 used the single figure-read "
    "constants those replaced, and gives the same answer.")
hdr = ["Receptor", "Constant", "Receptor residue", "rho (RMSF vs constant)",
       "distance to Ang-2 (nm)", "contact occupancy", "perfect correlation"]
for j, h in enumerate(hdr, 1):
    c = s.cell(row=5, column=j)
    c.value = h
    c.font = Font(bold=True)
for i, r in enumerate(rows):
    for j, v in enumerate(r + ["yes" if abs(r[3]) >= 0.9999 else ""], 1):
        s.cell(row=6 + i, column=j).value = v
r0 = 6 + len(rows) + 1
for rec in ("ILDR2", "LSR"):
    for const in ("ka", "kd"):
        sub = [x for x in rows if x[0] == rec and x[1] == const]
        hits = sum(1 for x in sub if abs(x[3]) >= 0.9999)
        for j, v in enumerate([rec, const, "TOTAL", hits, len(sub),
                               round(FLOOR * len(sub), 3),
                               "found / scanned / expected by chance"], 1):
            s.cell(row=r0, column=j).value = v
        r0 += 1
for col, wd in zip("ABCDEFG", (9, 9, 17, 22, 21, 17, 34)):
    s.column_dimensions[col].width = wd
w.save(XL)
print("sheet 14 written, %d rows; sheets now %d" % (len(rows), len(w.sheetnames)))
