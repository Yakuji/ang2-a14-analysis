"""Tyr99 across all six variants in both receptors.

Occupancy cannot answer the question: Tyr99 is in contact in 1.000 of frames in all twelve
systems with zero on-off switches, which is why it was excluded from the grouping sweep as a
degenerate observable. A saturated indicator says nothing about how much contact there is.

So contact is counted rather than detected: the number of heavy-atom pairs within 4.5 A
between Tyr99 and the receptor, the number of distinct receptor residues it touches, and its
mean minimum distance. These have range even when occupancy does not.

Same trajectories, stride and contact definition as occ3.py, and the replica is the unit, so
these numbers sit beside the rest of the paper's analysis rather than beside nothing.
"""
import os
import re
import glob
import numpy as np
from multiprocessing import Pool
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

OUT = "/home/yuki/y99"
os.makedirs(OUT, exist_ok=True)
CUT = 4.5
STRIDE = 1
RESID = 99
R1 = {"Lwt": "v2wt", "Lg": "v2g", "Lc": "v2c", "Ls": "v2s", "Lt": "v2t",
      "Lv": "v2v", "iwt": "iwt", "ig": "ig", "ic": "ic", "is": "is",
      "it": "it", "iv": "iv"}


def grp(path, name):
    idx, grab = [], False
    for line in open(path):
        s = line.strip()
        if s.startswith("["):
            grab = (s.strip("[] ").strip() == name)
            continue
        if grab and s:
            idx += [int(x) - 1 for x in s.split()]
    return np.array(sorted(set(idx)), dtype=int)


def one(job):
    name, tpr, xtc, ndx = job
    dst = os.path.join(OUT, name + ".npz")
    if os.path.exists(dst):
        return (name, "cached")
    try:
        u = mda.Universe(tpr, xtc)
        A = u.atoms[grp(ndx, "ChainA")]
        B = u.select_atoms("protein").subtract(A)
        y = [r for r in A.residues if int(r.resid) == RESID]
        if not y:
            return (name, "no residue %d" % RESID)
        yh = y[0].atoms.select_atoms("not name H*")
        assert y[0].resname in ("TYR",), y[0].resname
        Bres = [r for r in B.residues]
        Bh = [r.atoms.select_atoms("not name H*") for r in Bres]
        npair, nres, dmin = [], [], []
        for ts in u.trajectory[::STRIDE]:
            per = [distance_array(yh.positions, b.positions, box=ts.dimensions)
                   for b in Bh]
            cnt = sum(int((d < CUT).sum()) for d in per)
            touched = sum(1 for d in per if d.min() < CUT)
            npair.append(cnt)
            nres.append(touched)
            dmin.append(min(d.min() for d in per))
        np.savez_compressed(dst, npair=np.array(npair), nres=np.array(nres),
                            dmin=np.array(dmin), resname=y[0].resname,
                            nframes=len(npair))
        return (name, "ok %d frames, pairs %.1f" % (len(npair), np.mean(npair)))
    except Exception as e:
        return (name, "FAIL %s: %s" % (type(e).__name__, e))


def jobs():
    J = []
    for tag in R1:
        ndx = "/home/yuki/pmf/smd_an/%s/chainA_prot.ndx" % tag
        p = "/home/yuki/md_project_%s" % R1[tag]
        if os.path.exists(p + "/protein_cluster.xtc"):
            J.append((tag + "_r1", "/home/yuki/pmf/smd_an/%s/prot.tpr" % tag,
                      p + "/protein_cluster.xtc", ndx))
        dirs = sorted(glob.glob("/home/yuki/variant_replicas/runs/%s_r*" % tag))
        dirs += sorted(glob.glob("/home/yuki/md_replicas/runs/%s_r*" % tag))
        for d in dirs:
            rep = re.search(r"_r(\d+)$", d).group(1)
            t = d + "/prot_only.tpr"
            if not os.path.exists(t):
                t = d + "/prot.tpr"
            J.append((tag + "_r" + rep, t, d + "/whole.xtc", ndx))
    return J


if __name__ == "__main__":
    J = jobs()
    print("trajectories:", len(J), flush=True)
    ok = 0
    with Pool(6) as p:
        for name, st in p.imap_unordered(one, J):
            if st.startswith("ok") or st == "cached":
                ok += 1
            else:
                print(" ", name, st, flush=True)
    print("done %d / %d" % (ok, len(J)))
