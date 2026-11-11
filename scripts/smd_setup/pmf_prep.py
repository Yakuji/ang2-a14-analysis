"""Prepare one system for the association PMF.

  usage: pmf_prep.py SRC_GRO NDX OUTDIR PAD PULL

Takes the PBC-whole final frame of the equilibrium MD, rotates the complex so the
chain A -> chain B centre-of-mass vector lies along +z, writes the protein-only
coordinates in topology order, and emits the box dimensions, an index file and the
pbcatom choices needed by the pull code.

Selections are by ATOM INDEX from index_chains.ndx: chain B restarts its residue
numbering at 1, so residue-based selections match atoms in both chains.
"""
import sys, os
import numpy as np
import MDAnalysis as mda

SRC, NDX, OUT, PAD, PULL = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
os.makedirs(OUT, exist_ok=True)


def read_ndx(path, name):
    idx, grab = [], False
    for line in open(path):
        s = line.strip()
        if s.startswith("["):
            grab = (s.strip("[] ").strip() == name)
            continue
        if grab and s:
            idx += [int(x) - 1 for x in s.split()]
    return np.array(idx, dtype=int)


u = mda.Universe(SRC)
ia, ib = read_ndx(NDX, "ChainA"), read_ndx(NDX, "ChainB")
assert ia.max() < ib.min(), "chain A must precede chain B in the coordinate file"
prot = u.atoms[np.concatenate([ia, ib])]

A, B = u.atoms[ia], u.atoms[ib]
v = B.center_of_mass() - A.center_of_mass()
z = v / np.linalg.norm(v)
t = np.array([1.0, 0.0, 0.0])
if abs(np.dot(t, z)) > 0.9:
    t = np.array([0.0, 1.0, 0.0])
x = np.cross(t, z); x /= np.linalg.norm(x)
y = np.cross(z, x)
R = np.vstack([x, y, z])                      # rows are the new basis vectors

pos = prot.positions @ R.T                    # rotate into the pull frame
pos -= pos.min(axis=0)
prot.positions = pos

span = (pos.max(axis=0) - pos.min(axis=0)) / 10.0
box = span + 2 * PAD
box[2] += PULL

prot.write(os.path.join(OUT, "protein.gro"))

nA, nB = len(ia), len(ib)
with open(os.path.join(OUT, "chains.ndx"), "w") as f:
    for name, lo, hi in (("ChainA", 1, nA), ("ChainB", nA + 1, nA + nB)):
        f.write(f"[ {name} ]\n")
        nums = list(range(lo, hi + 1))
        for i in range(0, len(nums), 15):
            f.write(" ".join(f"{k:6d}" for k in nums[i:i + 15]) + "\n")
    f.write("[ Protein ]\n")
    nums = list(range(1, nA + nB + 1))
    for i in range(0, len(nums), 15):
        f.write(" ".join(f"{k:6d}" for k in nums[i:i + 15]) + "\n")

# pbcatom for each pull group: the atom closest to that group's centre of mass
def pbcatom(sel, offset):
    p = sel.positions @ R.T
    p -= p.mean(axis=0)
    return int(np.argmin(np.linalg.norm(p, axis=1))) + 1 + offset

pa, pb = pbcatom(A, 0), pbcatom(B, nA)

# alpha carbons of chain B, for the reference-chain position restraints
ca = [i + 1 for i, at in enumerate(u.atoms[ib]) if at.name == "CA"]
with open(os.path.join(OUT, "posre_B_ca.itp"), "w") as f:
    f.write("; alpha-carbon restraints on the reference chain, applied during pulling\n")
    f.write("; and umbrella sampling so the pull axis stays defined. Side chains are free.\n")
    f.write("[ position_restraints ]\n")
    for i in ca:
        f.write(f"{i:6d}     1  1000  1000  1000\n")

# The rotation puts the mass-weighted COM separation exactly along +z, so the z-component
# the pull code will see equals |v|. Taken from v, which was captured before
# prot.positions was overwritten.
d0 = float(np.linalg.norm(v)) / 10.0
with open(os.path.join(OUT, "geom.txt"), "w") as f:
    f.write(f"BOXX {box[0]:.3f}\nBOXY {box[1]:.3f}\nBOXZ {box[2]:.3f}\n")
    f.write(f"NA {nA}\nNB {nB}\nPBCA {pa}\nPBCB {pb}\nD0 {d0:.4f}\nNCA {len(ca)}\n")
print(f"box {box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} nm   d0(z) {d0:.3f} nm   "
      f"nA {nA} nB {nB}  pbcatom {pa}/{pb}  CA {len(ca)}")
