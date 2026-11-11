"""Force profiles of all 216 pulls, an observable the breaking-order analysis discarded.

A steered pull at 0.5 nm/ns is far from equilibrium, so none of these numbers is a free
energy and absolute values carry no thermodynamic meaning. What they can do is compare
systems run under one identical protocol.

Per pull:
  Fmax      the largest restoring force, the height of the resistance to separation
  work      the integral of F dx, the non-equilibrium work done pulling the complex apart
  dx_at_max the separation at which resistance peaks
  npeaks    local maxima above half of Fmax, at least 0.1 nm apart, which asks whether the
            complex comes apart in one step or several
  F_late    mean force over the last quarter, after the main rupture

By microscopic reversibility the last contacts to break are the first to form, so the late
part of the profile is the part that bears on association and the peak on dissociation.

GROMACS reports the force on the pull coordinate, negative while the spring is stretched
against the complex; magnitudes are used.
"""
import os
import numpy as np

ROOT = "/home/yuki/vpull/runs"
OUT = "/home/yuki/force_prof.npz"


def read_xvg(p):
    t, v = [], []
    with open(p) as f:
        for line in f:
            if line[:1] in ("#", "@"):
                continue
            a = line.split()
            if len(a) >= 2:
                t.append(float(a[0]))
                v.append(float(a[1]))
    return np.array(t), np.array(v)


def peaks(x, F, frac=0.5, sep=0.1):
    thr = frac * F.max()
    out = []
    for i in range(1, len(F) - 1):
        if F[i] >= thr and F[i] >= F[i - 1] and F[i] >= F[i + 1]:
            if not out or abs(x[i] - x[out[-1]]) >= sep:
                out.append(i)
            elif F[i] > F[out[-1]]:
                out[-1] = i
    return out


rows = {}
for tag in sorted(os.listdir(ROOT)):
    for conf in sorted(os.listdir(os.path.join(ROOT, tag))):
        d = os.path.join(ROOT, tag, conf)
        fx = os.path.join(d, "pullx.xvg")
        ff = os.path.join(d, "pullf.xvg")
        if not (os.path.exists(fx) and os.path.exists(ff)):
            continue
        _, x = read_xvg(fx)
        _, F = read_xvg(ff)
        n = min(len(x), len(F))
        x = x[:n]
        F = np.abs(F[:n])
        k = int(np.argmax(F))
        try:
            w = float(np.trapezoid(F, x))
        except AttributeError:
            w = float(np.trapz(F, x))
        q = int(0.75 * n)
        rows[tag + "__" + conf] = np.array(
            [F.max(), w, x[k] - x[0], len(peaks(x, F)), F[q:].mean(), x[-1] - x[0]])

np.savez(OUT, **rows)
print("profiles:", len(rows))
k0 = sorted(rows)[0]
r = rows[k0]
print("  example", k0)
print("    Fmax %.0f kJ/mol/nm  work %.0f kJ/mol  dx@max %.2f nm  peaks %d  F_late %.0f  dx %.2f nm"
      % (r[0], r[1], r[2], int(r[3]), r[4], r[5]))
