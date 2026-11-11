#!/bin/bash
# Package for the variant conformational pulls.
#
#   usage: mkpull_var.sh TAG [TAG ...]
#
# Rebuilt from scratch each time so adding a system cannot leave a stale copy behind. The
# rotation, box and chain indices were computed locally by frames_var.sh; what ships is the
# protein coordinates plus the topology, and the pod only needs GROMACS.
#
# GROMACS on the pod must be 2024.3. The wild-type 51-conformation set was computed entirely
# under 2024.3 and these pulls are compared against it; an earlier set run under 2026.3
# confounded a receptor comparison with the version difference and had to be set aside.
set -eu
OUT=/home/yuki/runpod_vpull
CONF=/home/yuki/confs_var
rm -rf $OUT
mkdir -p $OUT/scripts $OUT/mdp

cp /home/yuki/pmf/mdp/em.mdp /home/yuki/pmf/mdp/nvt.mdp /home/yuki/pmf/mdp/npt.mdp $OUT/mdp/
cp /home/yuki/runpod_pull/scripts/build_pull.sh $OUT/scripts/
cp /home/yuki/mdpkg/rp_vpull_queue.sh $OUT/scripts/
cp /home/yuki/mdpkg/vpull_stat.sh $OUT/scripts/stat.sh
cp /home/yuki/runpod_vmd_E/scripts/setup.sh $OUT/scripts/ 2>/dev/null \
  || cp /home/yuki/mdpkg/pod_setup.sh $OUT/scripts/setup.sh
chmod +x $OUT/scripts/*.sh

r1dir() {
  case "$1" in
    Lwt) echo v2wt ;; Lg) echo v2g ;; Lc) echo v2c ;;
    Ls)  echo v2s  ;; Lt) echo v2t ;; Lv) echo v2v ;;
    *)   echo "$1" ;;
  esac
}

for TAG in "$@"; do
  [ -d $CONF/$TAG ] || { echo "  $TAG: no conformations"; continue; }
  SRC=/home/yuki/md_project_$(r1dir $TAG)
  D=$OUT/systems/$TAG
  mkdir -p $D
  # topology header only; the molecules section is written fresh by build_pull.sh, then
  # solvate and genion append the solvent and ion counts
  sed '/^\[ molecules \]/,$d' $SRC/topol.top > $D/topol_head.top
  cp $SRC/topol_Protein_chain_A.itp $SRC/topol_Protein_chain_B.itp $D/
  cp $SRC/posre_Protein_chain_A.itp $SRC/posre_Protein_chain_B.itp $D/
  # the reference-chain alpha-carbon restraints are per conformation, so chain B's topology
  # has to be able to include them
  cat >> $D/topol_Protein_chain_B.itp <<'ITP'

#ifdef POSRES_B_CA
#include "posre_B_ca.itp"
#endif
ITP
  n=0
  for C in $CONF/$TAG/r*ns; do
    name=$(basename $C)
    O=$OUT/confs/$TAG/$name
    mkdir -p $O
    cp $C/protein.gro $C/chains.ndx $C/geom.txt $C/posre_B_ca.itp $O/
    n=$((n+1))
  done
  echo "  $TAG: $n conformations, $(du -sh $OUT/confs/$TAG | cut -f1)"
done
echo
du -sh $OUT
