#!/bin/bash
# Starting conformations for the variant pulls: six time points from each of three MD
# replicas, per system.
#
#   usage: frames_var.sh TAG [TAG ...]
#
# Six points at 20, 35, 50, 65, 80, 95 ns are a subgrid of the 5 ns grid the wild-type set
# used, so the matching eighteen wild-type conformations already exist and wt needs no new
# pulls. Eighteen rather than sixteen because sixteen does not divide by three and would
# leave the replicas unbalanced; eighteen gives SEM 12/sqrt(18) = 2.83 against a target of 3.
#
# Every frame is dumped with -pbc cluster, not -pbc whole. `whole` rebuilds molecules that
# the periodic boundary split but leaves the two chains in whatever images they ended up in;
# six of the first seventy-two conformations came out with the chains an image apart, giving
# a chain A to chain B separation of 8-11 nm instead of 3.5. The pull geometry is derived
# from that vector, so those conformations would have been set up around empty space.
#
# The d0 guard below is the check that catches it without anyone reading the numbers. A bound
# complex sits at 3.2-3.8 nm; anything outside 3.0-5.0 means the frame is not a bound complex
# and the conformation is discarded rather than shipped.
source /home/yuki/gromacs-2024.3/build/scripts/GMXRC.bash
PY=/home/yuki/miniforge3/envs/aptamd/bin/python
OUT=/home/yuki/confs_var
# Six points at 15 ns spacing by default. Both can be overridden per system from the
# environment, which is how A14C/ILDR2 is built: its r1 is being set aside, so it draws
# eight points from r2 and r4 instead of six from three replicas.
POINTS="${POINTS:-20 35 50 65 80 95}"
ONLY_REPS="${ONLY_REPS:-}"
D0MIN=3.0
D0MAX=5.0
mkdir -p $OUT

r1dir() {
  case "$1" in
    Lwt) echo v2wt ;; Lg) echo v2g ;; Lc) echo v2c ;;
    Ls)  echo v2s  ;; Lt) echo v2t ;; Lv) echo v2v ;;
    *)   echo "$1" ;;
  esac
}

BAD=0
for TAG in "$@"; do
  A=/home/yuki/pmf/smd_an/$TAG
  NDX=$OUT/${TAG}_chains.ndx
  P1=/home/yuki/md_project_$(r1dir $TAG)
  [ -f $P1/protein_cluster.xtc ] || { echo "=== $TAG: no r1 trajectory at $P1"; continue; }

  NA=$($PY -I /home/yuki/mdpkg/count_ndx.py $A/chainA_prot.ndx)
  NTOT=$(gmx check -f $P1/protein_cluster.xtc 2>&1 | awk '/# Atoms/{print $3}')
  $PY -I /home/yuki/mdpkg/write_chains.py "$NA" "$NTOT" "$NDX"
  echo "=== $TAG  protein $NTOT atoms, chain A $NA, chain B $((NTOT-NA))"

  # Replica numbers are discovered, not assumed. A14C against ILDR2 has r1, r2 and r4:
  # its r3 was lost with the pod that was running it, and a hard-coded 1 2 3 would have
  # silently produced twelve conformations for that system instead of eighteen.
  REPS="1"
  # The two wild types' replicas 2 and 3 were harvested into md_replicas, the variants' into
  # variant_replicas. Both are searched rather than assuming one tree.
  for d in /home/yuki/variant_replicas/runs/${TAG}_r* /home/yuki/md_replicas/runs/${TAG}_r*; do
    [ -d "$d" ] || continue
    REPS="$REPS $(basename "$d" | sed "s/.*_r//")"
  done
  REPS=$(echo $REPS | tr ' ' '
' | sort -n -u | tr '
' ' ')
  # Three replicas per system, lowest numbers first, and never more. The design is six
  # time points from each of three replicas for every system; a system that happens to
  # own a fourth would otherwise get 24 conformations and a precision the others do not
  # have. For A14C those three are r1, r2 and r4 - r3 was lost with its pod - and r11 is
  # the post-hoc sensitivity replica, which stays out of the primary set on both the MD
  # and the pull side so the two rest on the same trajectories.
  if [ -n "$ONLY_REPS" ]; then
    REPS="$ONLY_REPS"
  else
    REPS=$(echo $REPS | tr ' ' '
' | sort -n | head -3 | tr '
' ' ')
  fi
  echo "  replicas: $REPS(capped at 3)"
  for REP in $REPS; do
    if [ "$REP" = "1" ]; then
      X=$P1/protein_cluster.xtc; T=$A/prot.tpr
    else
      R=/home/yuki/variant_replicas/runs/${TAG}_r${REP}
      [ -d "$R" ] || R=/home/yuki/md_replicas/runs/${TAG}_r${REP}
      X=$R/whole.xtc
      T=$R/prot_only.tpr
      [ -f "$T" ] || T=$R/prot.tpr
    fi
    if [ ! -f "$X" ]; then echo "  r$REP: trajectory missing, skipped"; continue; fi
    XA=$(gmx check -f "$X" 2>&1 | awk '/# Atoms/{print $3}')
    if [ "$XA" != "$NTOT" ]; then
      echo "  r$REP: $XA atoms, expected $NTOT - SKIPPED (chain index would not apply)"
      continue
    fi
    for NS in $POINTS; do
      C=$OUT/$TAG/r${REP}_${NS}ns
      rm -rf $C && mkdir -p $C
      printf 'Protein\nProtein\n' | gmx trjconv -s "$T" -f "$X" -dump $((NS*1000)) \
          -pbc cluster -o $C/frame.gro > $C/trjconv.log 2>&1
      if [ ! -f $C/frame.gro ]; then
        echo "  r$REP ${NS}ns trjconv FAILED"; tail -3 $C/trjconv.log; rm -rf $C; BAD=$((BAD+1)); continue
      fi
      $PY /home/yuki/pmf/pmf_prep.py $C/frame.gro "$NDX" $C 1.6 1.6 > $C/prep.log 2>&1 \
        || { echo "  r$REP ${NS}ns prep FAILED"; tail -3 $C/prep.log; rm -rf $C; BAD=$((BAD+1)); continue; }
      D0=$(awk '/^D0/{print $2}' $C/geom.txt)
      BZ=$(awk '/^BOXZ/{print $2}' $C/geom.txt)
      OK=$(awk -v d="$D0" -v lo="$D0MIN" -v hi="$D0MAX" 'BEGIN{print (d>=lo && d<=hi) ? "y" : "n"}')
      if [ "$OK" != "y" ]; then
        printf "  r%s %3s ns  d0 %6s nm  REJECTED (outside %s-%s nm, chains not a bound complex)\n" \
          "$REP" "$NS" "$D0" "$D0MIN" "$D0MAX"
        rm -rf $C; BAD=$((BAD+1)); continue
      fi
      printf "  r%s %3s ns  d0 %6s nm  boxz %6s nm\n" "$REP" "$NS" "$D0" "$BZ"
    done
  done
  echo "  -> $(ls -1d $OUT/$TAG/r*ns 2>/dev/null | wc -l)/18 conformations, $(du -sh $OUT/$TAG | cut -f1)"
  echo
done
echo "rejected or failed: $BAD"
du -sh $OUT
