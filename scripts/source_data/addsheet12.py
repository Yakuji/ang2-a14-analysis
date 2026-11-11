# -*- coding: utf-8 -*-
"""Add the opening-frequency sheets to the Prism workbook, and flag sheet 7 as an average.

Sheet 7 holds the time-averaged exposure, which is what the first draft plotted. That
average is taken over a residue whose median exposure is zero in all twelve systems, so it
describes no structure. The two sheets added here carry what Figure 4 now shows: the
fraction of frames in which residue 14 has any exposure at all (panel A), and the per-frame
distribution behind it (panel C).
"""
import numpy as np
import openpyxl
from openpyxl.styles import Font

XL = "/mnt/c/Users/y-niw/OneDrive/Desktop/Ang2_MD/FIGURE_DATA_2026-10-24.xlsx"
SUF = ["wt", "g", "c", "s", "t", "v"]
GAM = {"wt": 0, "g": 0, "c": 1, "s": 1, "t": 2, "v": 2}
LAB = {"wt": "WT", "g": "A14G", "c": "A14C", "s": "A14S", "t": "A14T", "v": "A14V"}
D = np.load("/home/yuki/res14_open.npz")
B = Font(bold=True)

w = openpyxl.load_workbook(XL)
for n in ("12_MD_res14_open_fraction", "13_MD_res14_distribution"):
    if n in w.sheetnames:
        del w[n]

# sheet 7 gets a line saying what the average is an average of
s7 = w["7_MD_res14_exposure"]
s7.cell(row=2, column=1).value = (
    "The median of this quantity is zero in all twelve systems: residue 14 is shut for "
    "51-88% of frames and open to 38-66% the rest of the time. These are therefore averages "
    "over an open-and-shut residue. Sheet 12 has the frequency, sheet 13 the distribution.")

# ---------------------------------------------------------------- sheet 12
s = w.create_sheet("12_MD_res14_open_fraction")
s.cell(row=1, column=1).value = (
    "How often residue 14 has any exposure at all in the complex, per MD replica. This is "
    "Figure 4A, and it is the quantity that separates the groups: ILDR2 15.2% for the two "
    "gamma-free systems against 41.2% for the four impaired ones, +2.20 SE with no overlap "
    "(p = 0.0667, the permutation floor for a binary split of six systems).")
s.cell(row=2, column=1).value = (
    "LSR gives no reading: the groups overlap (28.0% against 35.2%) and the replicate "
    "scatter is 0.635 of the spread across variants, past the 0.50 resolution guard. "
    "SE here is the project's gap statistic - the gap between the nearest group members "
    "over the mean standard error of a system mean.")
hdr = ["Receptor", "Variant", "gamma atoms", "MD replica", "frames",
       "frames with exposure", "open fraction (%)", "exposure when open (% of free)",
       "time-averaged exposure (% of free)", "median exposure (% of free)"]
for j, h in enumerate(hdr, 1):
    c = s.cell(row=4, column=j)
    c.value = h
    c.font = B
r = 5
for rec in ("ILDR2", "LSR"):
    for suf in SUF:
        a = D["%s_%s" % (rec, suf)]
        n = len(D["pool_%s_%s" % (rec, suf)]) // len(a)
        for row in a:
            vals = [rec, LAB[suf], GAM[suf], int(row[0]), n,
                    int(round(n * row[4] / 100.0)), round(float(row[4]), 2),
                    round(float(row[5]), 2), round(float(row[3]), 2),
                    round(float(row[6]), 3)]
            for j, v in enumerate(vals, 1):
                s.cell(row=r, column=j).value = v
            r += 1
for col, wd in zip("ABCDEFGHIJ", (9, 8, 12, 11, 8, 20, 16, 27, 30, 24)):
    s.column_dimensions[col].width = wd

# ---------------------------------------------------------------- sheet 13
s = w.create_sheet("13_MD_res14_distribution")
s.cell(row=1, column=1).value = (
    "Per-frame distribution of residue 14's exposure, pooled over the three replicas of "
    "each system (3003 frames each). This is Figure 4C: a spike at zero and a separate "
    "open mode, with nothing in between - which is why the time average in sheet 7 is not "
    "a description of a structure.")
s.cell(row=2, column=1).value = (
    "Each column is the percentage of that system's frames falling in the bin. The first "
    "row is exposure exactly zero; the remaining rows are 2.5%-wide bins; the last row "
    "collects everything above 100%. A bound frame can exceed 100% because the reference "
    "is the free ligand's time average, and 4-12% of frames in the gamma-bearing systems "
    "do (the per-frame maximum is 328%).")
bins = np.linspace(0, 100, 41)
cols = [("%s %s" % (rec, LAB[suf]), D["pool_%s_%s" % (rec, suf)])
        for rec in ("ILDR2", "LSR") for suf in SUF]
s.cell(row=4, column=1).value = "bin from (%)"
s.cell(row=4, column=2).value = "bin to (%)"
for j, (name, _) in enumerate(cols, 3):
    s.cell(row=4, column=j).value = name
for j in range(1, len(cols) + 3):
    s.cell(row=4, column=j).font = B
s.cell(row=5, column=1).value = 0.0
s.cell(row=5, column=2).value = 0.0
for j, (_, p) in enumerate(cols, 3):
    s.cell(row=5, column=j).value = round(100.0 * float((p <= 1e-9).mean()), 3)
for i in range(len(bins) - 1):
    lo, hi = bins[i], bins[i + 1]
    s.cell(row=6 + i, column=1).value = round(float(lo), 2)
    s.cell(row=6 + i, column=2).value = round(float(hi), 2)
    for j, (_, p) in enumerate(cols, 3):
        m = (p > 1e-9) & (p > lo) & (p <= hi) if i else (p > 1e-9) & (p <= hi)
        s.cell(row=6 + i, column=j).value = round(100.0 * float(m.mean()), 3)
row = 6 + len(bins) - 1
s.cell(row=row, column=1).value = 100.0
s.cell(row=row, column=2).value = "and above"
for j, (_, p_) in enumerate(cols, 3):
    s.cell(row=row, column=j).value = round(100.0 * float((p_ > 100.0).mean()), 3)
s.column_dimensions["A"].width = 12
s.column_dimensions["B"].width = 12
for j in range(3, len(cols) + 3):
    s.column_dimensions[openpyxl.utils.get_column_letter(j)].width = 13

w.save(XL)
print("sheets now:", ", ".join(w.sheetnames))
# the columns must each sum to 100
for name, p in cols:
    tot = 100.0 * float((p <= 1e-9).mean()) + sum(
        100.0 * float((((p > 1e-9) & (p > bins[i]) & (p <= bins[i + 1])) if i
                       else ((p > 1e-9) & (p <= bins[1]))).mean())
        for i in range(len(bins) - 1)) + 100.0 * float((p > 100.0).mean())
    assert abs(tot - 100.0) < 0.05, (name, tot)
print("all 12 distribution columns sum to 100%")
