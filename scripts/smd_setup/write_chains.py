"""ChainA/ChainB in protein-only numbering: chain A holds the first nA atoms and chain B the
rest, because convert-tpr preserves topology order and chain A precedes chain B."""
import sys
na, ntot, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(out, "w") as f:
    for name, lo, hi in (("ChainA", 1, na), ("ChainB", na + 1, ntot)):
        f.write(f"[ {name} ]\n")
        nums = list(range(lo, hi + 1))
        for i in range(0, len(nums), 15):
            f.write(" ".join(f"{k:6d}" for k in nums[i:i + 15]) + "\n")
