#!/bin/bash
# S2 MNG slab: run Palace for every run listed in a list file, then summarize (inside the job).
#
# Submit from this directory on LANTA (account, partition and install path are given at submit time):
#   sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=<N> -t <HH:MM:SS> \
#          --export=ALL,PALACE_ENV=<PROJECT_DIR>/software/palace/env.sh run_s2.sbatch lists/<round>.txt
# The list file holds one run name per line; configs/<run>.json, runs/<run>.json, cases/<mesh>.json
# and mesh/<mesh>.msh must exist (made locally by s2_geometry.py).
# Output: postpro/<run>/ (Palace CSVs), results/<run>_*.csv (summaries, S-parameters).
#
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -J s2-slab
#SBATCH -o slurm-%j.out
set -uo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"
cd "$SLURM_SUBMIT_DIR"
LIST=${1:?list file}
ROUND=$(basename "$LIST" .txt)
mkdir -p logs results
names=$(grep -v '^\s*#' "$LIST" | awk 'NF {print $1}')
echo "ntasks=$SLURM_NTASKS palace=${PALACE_VERSION:-?} list=$LIST" | tee -a "results/$ROUND-timings.txt"
for n in $names; do
  t0=$(date +%s.%N)
  srun -n "$SLURM_NTASKS" palace-x86_64.bin "configs/$n.json" > "logs/$n.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  awk -v c="$n" -v a="$t0" -v b="$t1" -v rc="$rc" -v np="$SLURM_NTASKS" \
    'BEGIN { printf "%s rc=%d ntasks=%d wall_s=%.1f\n", c, rc, np, b - a }' | tee -a "results/$ROUND-timings.txt"
  # keep the solver summary (sample count of the adaptive sweep, solve times) next to the results
  grep -E -i -A3 "adaptive sampling|unknowns|Elapsed Time Report|error|warning" "logs/$n.log" | grep -v "^It " | head -n 80 > "results/$n-log-excerpt.txt"
  python3 postprocess_s2.py "$n"
done
