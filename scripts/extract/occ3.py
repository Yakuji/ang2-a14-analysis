"""Per-residue contact occupancy on the Ang-2 side, all 120 residues, three replicas.

Descriptive, not a correlation hunt: which contacts actually exist, and for what fraction of
the trajectory. A mean distance hides whether a contact is permanent or intermittent, and
every distance-based result in this project so far has used the mean.

Per complex trajectory, for every chain A residue:
  occ    fraction of frames with any heavy atom within 4.5 A of the receptor
  dmin   mean closest heavy-atom distance to the receptor
  nsw    number of formed/broken transitions, as a measure of intermittency

Minimum image is applied through distance_array(box=), so this is independent of how the
trajectory was made whole.
"""
import os
import glob
import re
import numpy as np
from multiprocessing import Pool
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

STRIDE = 10
CUT = 4.5
OUT = "/home/yuki/occ3"
os.makedirs(OUT, exist_ok=True)
R1 = {"Lwt": "v2wt", "Lg": "v2g", "Lc": "v2c", "Ls": "v2s", "Lt": "v2t", "Lv": "v2v",
      "iwt": "iwt", "ig": "ig", "ic": "ic", "is": "is", "it": "it", "iv": "iv"}


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
        Bh = B.select_atoms("not name H*")
        sel = [r.atoms.select_atoms("not name H*") for r in A.residues]
        rows = []
        for ts in u.trajectory[::STRIDE]:
            rows.append([distance_array(s.positions, Bh.positions,
                                        box=ts.dimensions).min() for s in sel])
        D = np.asarray(rows)
        inc = D < CUT
        nsw = np.abs(np.diff(inc.astype(np.int8), axis=0)).sum(0)
        np.savez_compressed(dst, resids=A.residues.resids,
                            resnames=np.array(A.residues.resnames),
                            occ=inc.mean(0), dmin=D.mean(0),
                            dsd=D.std(0), nsw=nsw, nframes=len(D))
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
