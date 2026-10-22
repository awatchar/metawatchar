#!/bin/bash
# Palace validation: eigenmodes of a Teflon-filled cylindrical cavity (PEC and copper walls),
# compared against analytic resonant frequencies and Q factors (compare_cavity.py).
#
# Submit from this directory on LANTA (account and install path are given at submit time):
#   sbatch -A <ACCOUNT> --export=ALL,PALACE_ENV=<PROJECT_DIR>/software/palace/env.sh run_cavity.sbatch
#
#SBATCH -p compute-limited
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH -t 00:45:00
#SBATCH -J palace-cavity
#SBATCH -o slurm-%j.out
set -euo pipefail
: "${PALACE_ENV:?set PALACE_ENV to the env.sh of the Palace install}"
source "$PALACE_ENV"            # loads modules, puts palace on PATH, sets PALACE_EXAMPLES
cd "$SLURM_SUBMIT_DIR"

# Meshes come from the Palace source tree (examples/cylinder/mesh, Apache-2.0).
mkdir -p mesh logs
cp "$PALACE_EXAMPLES/cylinder/mesh/cylinder_hex.msh" "$PALACE_EXAMPLES/cylinder/mesh/cylinder_prism.msh" mesh/

{
  echo "ntasks=$SLURM_NTASKS palace=$PALACE_VERSION"
  for cfg in cavity_pec cavity_impedance; do
    t0=$(date +%s.%N)
    srun -n "$SLURM_NTASKS" palace-x86_64.bin "$cfg.json" > "logs/$cfg.log" 2>&1
    t1=$(date +%s.%N)
    awk -v c="$cfg" -v a="$t0" -v b="$t1" 'BEGIN { printf "%s wall_s=%.1f\n", c, b - a }'
  done
} | tee timings.txt

# Post-processing (light, inside the job): analytic comparison table.
python3 compare_cavity.py postpro/cavity_pec postpro/cavity_impedance | tee results.md
