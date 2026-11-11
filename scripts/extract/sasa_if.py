"""Solvent-accessible surface area: a direct test of the steric hypothesis for ka.

The gamma-atom count at position 14 orders ILDR2's ka monotonically. If that works through
steric obstruction of the association path, the gamma atoms have to be in the path -- solvent
exposed in the free state, not packed against the ligand's own body. If instead position 14's
side chain is buried in the free protein, it cannot obstruct an incoming receptor and the
correlation needs a different mechanism.

Computed with gmx sasa, which is the Eisenhaber algorithm, on:
  free state    surface = the single chain; output = residue 14 whole residue, residue 14
                side chain, and the 18 interface residues
  bound state   surface = the whole complex, so the receptor occludes the epitope; same
                output groups, all on chain A

The differences give the buried area on binding. Glycine has no side-chain heavy atom beyond
CA, so its side-chain group is empty and recorded as zero rather than run through gmx.

Prediction, stated before the numbers: side-chain SASA at position 14 in the free state should
RISE with gamma-atom count (0 for Gly, then Ala, Cys/Ser, Thr/Val) and should correlate with
ILDR2's ka in the same direction the gamma count does, negatively. A buried side chain, or a
SASA that does not order with ka, removes the steric reading.
"""
import os
import glob
import re
import subprocess
import numpy as np
from multiprocessing import Pool
import MDAnalysis as mda

OUT = "/home/yuki/sasa"
os.makedirs(OUT, exist_ok=True)
GMX = "/home/yuki/gromacs-2024.3/build/bin/gmx"
R1 = {"Lwt": "v2wt", "Lg": "v2g", "Lc": "v2c", "Ls": "v2s", "Lt": "v2t", "Lv": "v2v",
      "iwt": "iwt", "ig": "ig", "ic": "ic", "is": "is", "it": "it", "iv": "iv"}
FREE = {"Lwt": "wt", "Lg": "g", "Lc": "c", "Ls": "s", "Lt": "t", "Lv": "v"}
IFACE = [12, 13, 14, 16, 18, 20, 21, 22, 23, 25, 27, 72, 73, 98, 99, 118, 119, 120]
DT = 100  # ps between analysed frames


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


def write_ndx(path, groups):
    with open(path, "w") as f:
        for name, idx in groups:
            f.write("[ %s ]\n" % name)
            a = np.asarray(idx) + 1
            for i in range(0, len(a), 15):
                f.write(" ".join("%d" % v for v in a[i:i + 15]) + "\n")


def parse_xvg(path):
    rows = []
    for line in open(path):
        if line[0] in "#@":
            continue
        p = line.split()
        if len(p) >= 2:
            rows.append([float(x) for x in p])
    return np.asarray(rows)


def one(job):
    name, tpr, xtc, which, chainA_ndx = job
    dst = os.path.join(OUT, name + ".npz")
    if os.path.exists(dst):
        return (name, "cached")
    wd = os.path.join(OUT, "w_" + name)
    os.makedirs(wd, exist_ok=True)
    try:
        u = mda.Universe(tpr, xtc)
        if which == "bound":
            A = u.atoms[grp(chainA_ndx, "ChainA")]
            surf = u.select_atoms("protein")
        else:
            A = u.atoms
            surf = u.select_atoms("protein")
        r14 = A.residues[13]
        assert int(r14.resid) == 14, "residue 14 is not where expected"
        res14 = r14.atoms.select_atoms("not name H*")
        sc14 = r14.atoms.select_atoms("not name H* and not name N CA C O")
        iface = A.residues[[i - 1 for i in IFACE]].atoms.select_atoms("not name H*")
        groups = [("Surf", surf.indices), ("Res14", res14.indices),
                  ("Iface", iface.indices)]
        outs = ["Res14", "Iface"]
        if len(sc14) > 0:
            groups.append(("SC14", sc14.indices))
            outs.append("SC14")
        ndx = os.path.join(wd, "g.ndx")
        write_ndx(ndx, groups)
        xvg = os.path.join(wd, "a.xvg")
        cmd = [GMX, "sasa", "-s", tpr, "-f", xtc, "-n", ndx, "-surface", "Surf",
               "-output"] + outs + ["-o", xvg, "-dt", str(DT), "-nopbc"]
        r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True)
        if not os.path.exists(xvg):
            return (name, "FAIL gmx: " + (r.stderr or "")[-400:])
        d = parse_xvg(xvg)
        # columns: time, total surface area, then one per -output group
        res = dict(t=d[:, 0], total=d[:, 1], res14=d[:, 2], iface=d[:, 3])
        res["sc14"] = d[:, 4] if len(outs) == 3 else np.zeros(len(d))
        res["nsc"] = len(sc14)
        res["resname14"] = str(r14.resname)
        np.savez_compressed(dst, **res)
        return (name, "ok")
    except Exception as e:
        return (name, "FAIL " + type(e).__name__ + ": " + str(e))


def jobs():
    J = []
    for tag in R1:
        ndx = "/home/yuki/pmf/smd_an/" + tag + "/chainA_prot.ndx"
        p = "/home/yuki/md_project_" + R1[tag]
        if os.path.exists(p + "/protein_cluster.xtc"):
            J.append(("bound_" + tag + "_r1",
                      "/home/yuki/pmf/smd_an/" + tag + "/prot.tpr",
                      p + "/protein_cluster.xtc", "bound", ndx))
        dirs = sorted(glob.glob("/home/yuki/variant_replicas/runs/" + tag + "_r*"))
        dirs += sorted(glob.glob("/home/yuki/md_replicas/runs/" + tag + "_r*"))
        for d in dirs:
            rep = re.search(r"_r(\d+)$", d).group(1)
            t = d + "/prot_only.tpr"
            if not os.path.exists(t):
                t = d + "/prot.tpr"
            J.append(("bound_" + tag + "_r" + rep, t, d + "/whole.xtc", "bound", ndx))
    for tag in FREE:
        f = FREE[tag]
        J.append(("free_" + tag + "_r1",
                  "/home/yuki/md_project_mono_" + f + "/prot.tpr",
                  "/home/yuki/md_project_mono_" + f + "/prot_whole.xtc", "free", None))
        for d in sorted(glob.glob("/home/yuki/free_replicas/runs/mono_" + f + "_r*")):
            rep = re.search(r"_r(\d+)$", d).group(1)
            J.append(("free_" + tag + "_r" + rep, d + "/prot_only.tpr",
                      d + "/whole.xtc", "free", None))
    return J


if __name__ == "__main__":
    J = jobs()
    print("trajectories:", len(J), flush=True)
    ok = 0
    with Pool(6) as p:
        for name, st in p.imap_unordered(one, J):
            if st in ("ok", "cached"):
                ok += 1
                if ok % 8 == 0:
                    print("  ...", ok, "/", len(J), flush=True)
            else:
                print(" ", name, st, flush=True)
    print("done", ok, "/", len(J))
