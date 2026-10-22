#!/bin/bash
# LGR-A size sweep: run Palace eigenmode for every case listed in a list file, then summarize.
#
# Submit from this directory on LANTA (account and install path are given at submit time):
#   sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PROJECT_DIR>/software/palace/env.sh \
#          run_sweep.sbatch lists/<round>.txt
# The list file holds one case name per line; mesh/<name>.msh, configs/<name>.json and
# cases/<name>.json must exist (made locally by lgr_geometry.py). Output: results/<round>.csv.
#
#SBATCH -p compute-limited
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH -J lgr-sweep
#SBATCH -o slurm-%j.out
set -uo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"
cd "$SLURM_SUBMIT_DIR"
LIST=${1:?list file}
ROUND=$(basename "$LIST" .txt)
mkdir -p logs results
names=$(grep -v '^\s*#' "$LIST" | awk 'NF {print $1}')
echo "ntasks=$SLURM_NTASKS palace=$PALACE_VERSION list=$LIST" | tee "results/$ROUND-timings.txt"
for n in $names; do
  t0=$(date +%s.%N)
  srun -n "$SLURM_NTASKS" palace-x86_64.bin "configs/$n.json" > "logs/$n.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  awk -v c="$n" -v a="$t0" -v b="$t1" -v rc="$rc" 'BEGIN { printf "%s rc=%d wall_s=%.1f\n", c, rc, b - a }' \
    | tee -a "results/$ROUND-timings.txt"
  # resolved config, timing and solver summary are small: keep them next to the results
done
# Post-processing (light, inside the job)
python3 postprocess.py "results/$ROUND.csv" $names
