"""The receptor side of the interface, which no analysis in this project has looked at.

Every observable so far has been an Ang-2 residue. If the association barrier lives in a
rearrangement of the receptor rather than the ligand, it would show up here and nowhere else
that has been measured.

Per complex trajectory, for every receptor residue:
  dmin   mean closest heavy-atom distance to chain A
  occ    fraction of frames within 4.5 A of chain A
  rmsf   alpha-carbon RMSF after superposition on chain B alpha carbons

There is no free-receptor simulation, so whether binding rigidifies the receptor cannot be
asked. What can be asked is whether the receptor behaves differently depending on which
Ang-2 variant is bound, with the other receptor as control.

Receptor residue numbering restarts at 1 in these builds, so chain B is taken as everything
protein that is not in the chain A index, never by residue number.
"""
import os
import glob
import re
import numpy as np
from multiprocessing import Pool
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis import align

STRIDE = 20
CUT = 4.5
OUT = "/home/yuki/recv"
os.makedirs(OUT, exist_ok=True)
R1 = {"Lwt": "v2wt", "Lg": "v2g", "Lc": "v2c", "Ls": "v2s", "Lt": "v2t", "Lv": "v2v",
      "iwt": "iwt", "ig": "ig", "ic": "ic", "is": "is", "it": "it", "iv": "iv"}


def grp(path, name):
    idx = []
    grab = False
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
        Ah = A.select_atoms("not name H*")
        sel = [r.atoms.select_atoms("not name H*") for r in B.residues]
        ca = B.select_atoms("name CA")
        rows = []
        pos = []
        for ts in u.trajectory[::STRIDE]:
            rows.append([distance_array(s.positions, Ah.positions,
                                        box=ts.dimensions).min() for s in sel])
            pos.append(ca.positions.copy())
        D = np.asarray(rows)
        X = np.asarray(pos)
        ref = X.mean(0)
        fit = np.empty_like(X)
        for i in range(len(X)):
            R, _ = align.rotation_matrix(X[i] - X[i].mean(0), ref - ref.mean(0))
            fit[i] = (X[i] - X[i].mean(0)) @ R.T
        rmsf = np.sqrt(((fit - fit.mean(0)) ** 2).sum(2).mean(0))
        np.savez_compressed(dst, resids=B.residues.resids, dmin=D.mean(0),
                            occ=(D < CUT).mean(0), rmsf=rmsf,
                            ca_resids=ca.resids, nframes=len(D))
        return (name, "ok")
    except Exception as e:
        return (name, "FAIL " + type(e).__name__ + ": " + str(e))


def jobs():
    J = []
    for tag in R1:
        ndx = "/home/yuki/pmf/smd_an/" + tag + "/chainA_prot.ndx"
        p = "/home/yuki/md_project_" + R1[tag]
        if os.path.exists(p + "/protein_cluster.xtc"):
            J.append((tag + "_r1", "/home/yuki/pmf/smd_an/" + tag + "/prot.tpr",
                      p + "/protein_cluster.xtc", ndx))
        dirs = sorted(glob.glob("/home/yuki/variant_replicas/runs/" + tag + "_r*"))
        dirs += sorted(glob.glob("/home/yuki/md_replicas/runs/" + tag + "_r*"))
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
            if st in ("ok", "cached"):
                ok += 1
                if ok % 6 == 0:
                    print("  ...", ok, "/", len(J), flush=True)
            else:
                print(" ", name, st, flush=True)
    print("done", ok, "/", len(J))
