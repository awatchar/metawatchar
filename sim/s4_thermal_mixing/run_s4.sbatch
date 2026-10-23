#!/bin/bash
# S4 thermal mixing: run Palace (eigenmode) for every run in a list file, then summarize inside
# the job.
#
# Submit from this directory on the cluster (account, partition, resources and install path are
# given at submit time, or by a site-specific wrapper that sets CASE_DIR and calls this script):
#   sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=<N> -t <HH:MM:SS> \
#          --export=ALL,PALACE_ENV=<PALACE_INSTALL>/env.sh run_s4.sbatch lists/<round>.txt
# The list file holds one run name per line; configs/<run>.json, runs/<run>.json and
# mesh/<mesh>.msh must exist (made locally by s4_geometry.py / build_s4.py).
# Output: postpro/<run>/ (Palace CSVs), results/<run>_eig.csv, results/<round>-timings.txt.
#
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -J s4-mix
#SBATCH -o slurm-%j.out
set -uo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"
cd "${CASE_DIR:-$SLURM_SUBMIT_DIR}"
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
  grep -E -i "unknowns|converged eigenvalues|Elapsed Time Report|error|warning" "logs/$n.log" | head -n 40 > "results/$n-log-excerpt.txt"
  python3 postprocess_s4.py "$n"
done
