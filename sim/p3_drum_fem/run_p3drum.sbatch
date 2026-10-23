#!/bin/bash
# P3 drum-inductor FEM: run Palace (magnetostatic) for every run in a list file and keep the
# inductance matrices.
#
# Submit from this directory on the cluster (account, partition and install path at submit time):
#   sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=<N> -t <HH:MM:SS> \
#          --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_p3drum.sbatch lists/<round>.txt
# configs/<run>.json and mesh/<mesh>.msh must exist (made locally by drum_geometry.py).
# Output: results/<run>_M.csv (terminal-M.csv), results/<round>-timings.txt, logs/<run>.log.
#
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -J p3-drum
#SBATCH -o slurm-%j.out
set -uo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"
cd "$SLURM_SUBMIT_DIR"
LIST=${1:?list file}
ROUND=$(basename "$LIST" .txt)
mkdir -p logs results
names=$(tr -d '\r' < "$LIST" | grep -v '^\s*#' | awk 'NF {print $1}')
echo "ntasks=$SLURM_NTASKS list=$LIST" | tee -a "results/$ROUND-timings.txt"
for n in $names; do
  t0=$(date +%s.%N)
  srun -n "$SLURM_NTASKS" palace-x86_64.bin "configs/$n.json" > "logs/$n.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  awk -v c="$n" -v a="$t0" -v b="$t1" -v rc="$rc" -v np="$SLURM_NTASKS" \
    'BEGIN { printf "%s rc=%d ntasks=%d wall_s=%.1f\n", c, rc, np, b - a }' | tee -a "results/$ROUND-timings.txt"
  cp "postpro/$n/terminal-M.csv" "results/${n}_M.csv" 2>/dev/null || echo "$n: no terminal-M.csv" | tee -a "results/$ROUND-timings.txt"
  grep -E -i "unknowns|iterations|Elapsed Time|error|warning" "logs/$n.log" | head -n 30 > "results/$n-log-excerpt.txt"
done
