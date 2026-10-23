#!/bin/bash
# S2 fixture: run Palace for every run listed in a list file, then summarize (inside the job).
#
# Submit from this directory on LANTA (account, partition, resources and install path are given at
# submit time, or by a site-specific wrapper that sets CASE_DIR and calls this script):
#   sbatch -A <ACCOUNT> -p <PARTITION> --ntasks=<N> -t <HH:MM:SS> \
#          --export=ALL,PALACE_ENV=<PROJECT_DIR>/software/palace/env.sh run_fixture.sbatch lists/<round>.txt
# The list file holds one run name per line, optionally followed by the number of MPI ranks for that
# run (default: all of the job's tasks); RUN_TIMEOUT (seconds, optional) caps each run; configs/<run>.json, runs/<run>.json and the meshes
# (made locally by fixture_geometry.py) must exist.
# Output: postpro/<run>/ (Palace CSVs), results/<run>_*.csv (S-parameters, probes, summaries).
#
#SBATCH --cpus-per-task=1
#SBATCH -J s2-fixture
#SBATCH -o slurm-%j.out
set -uo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"
cd "${CASE_DIR:-$SLURM_SUBMIT_DIR}"
LIST=${1:?list file}
ROUND=$(basename "$LIST" .txt)
mkdir -p logs results
lines=$(grep -v '^\s*#' "$LIST" | awk 'NF {print $1 ":" ($2 ? $2 : 0)}')
echo "ntasks=$SLURM_NTASKS palace=${PALACE_VERSION:-?} list=$LIST" | tee -a "results/$ROUND-timings.txt"
for l in $lines; do
  n=${l%%:*}; np=${l##*:}; [[ $np == 0 ]] && np=$SLURM_NTASKS
  t0=$(date +%s.%N)
  ${RUN_TIMEOUT:+timeout $RUN_TIMEOUT} srun -n "$np" palace-x86_64.bin "configs/$n.json" > "logs/$n.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  awk -v c="$n" -v a="$t0" -v b="$t1" -v rc="$rc" -v np="$np" \
    'BEGIN { printf "%s rc=%d ntasks=%d wall_s=%.1f\n", c, rc, np, b - a }' | tee -a "results/$ROUND-timings.txt"
  grep -E -i -A3 "adaptive sampling|unknowns|Elapsed Time Report|error|warning" "logs/$n.log" | grep -v "^It " | head -n 80 > "results/$n-log-excerpt.txt"
  python3 postprocess_fixture.py "$n"
done
