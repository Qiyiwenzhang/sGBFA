#!/bin/bash
#SBATCH --job-name=rexample
#SBATCH --nodes=1
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --qos=Short
#SBATCH --mem=64G
#SBATCH --cpus-per-task=6
#SBATCH --time=24:00:00
#SBATCH --array=47
#SBATCH --output=logs/my_r_code_%j.log   # %j = job id

set -euo pipefail
module purge
module load r/4.5.0


mkdir -p logs results
srun Rscript --vanilla sim_outcome_low_outcome_10.R "$SLURM_ARRAY_TASK_ID"
