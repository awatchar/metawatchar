#!/bin/bash
# LGR-A driven/eigenmode study: run every Palace case listed in a list file, NPAR cases at a time
# (each with RANKS MPI ranks) inside one allocation, then summarize (postprocess.py).
#
# Submit from this directory on LANTA (account and install path are given at submit time):
#   sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PROJECT_DIR>/software/palace/env.sh \
#          run_cases.sbatch lists/<round>.txt
# Optional environment: RANKS (default 32), NPAR (default SLURM_NTASKS / RANKS).
# mesh/<mesh>.msh and configs/<name>.json must exist (made locally by build_cases.py).
#
#SBATCH -p compute-limited
#SBATCH -N 1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH -J lgr-a-driven
#SBATCH -o slurm-%j.out
set -uo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"
cd "$SLURM_SUBMIT_DIR"
LIST=${1:?list file}
ROUND=$(basename "$LIST" .txt)
RANKS=${RANKS:-32}
NPAR=${NPAR:-$(( SLURM_NTASKS / RANKS ))}
mkdir -p logs results
names=$(grep -v '^\s*#' "$LIST" | awk 'NF {print $1}')
echo "ntasks=$SLURM_NTASKS ranks=$RANKS npar=$NPAR palace=$PALACE_VERSION list=$LIST" | tee "results/$ROUND-timings.txt"

run_one() {
  local n=$1 t0 t1 rc
  t0=$(date +%s.%N)
  srun --exact -n "$RANKS" -c 1 palace-x86_64.bin "configs/$n.json" > "logs/$n.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  awk -v c="$n" -v a="$t0" -v b="$t1" -v rc="$rc" 'BEGIN { printf "%s rc=%d wall_s=%.1f\n", c, rc, b - a }' \
    >> "results/$ROUND-timings.txt"
}

for n in $names; do
  while (( $(jobs -rp | wc -l) >= NPAR )); do wait -n; done
  run_one "$n" &
  sleep 2
done
wait
# Post-processing (light, inside the job): eigenmode summaries and S-parameter CSVs
python3 postprocess.py "$ROUND" $names >> "results/$ROUND-timings.txt" 2>&1
