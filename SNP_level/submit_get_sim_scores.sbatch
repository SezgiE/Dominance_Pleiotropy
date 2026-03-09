#!/bin/bash
#SBATCH --job-name=sim_scores
#SBATCH --time=01:00:00
#SBATCH --array=1-22%22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=32G               
#SBATCH --output=logs/sim_scores/job_%A_%a.out
#SBATCH --error=logs/sim_scores/job_%A_%a.err

PYTHON_EXEC="$HOME/.conda/envs/main-py3/bin/python"

# Set up paths
HOME_DIR="/home/sercan/dominance/overlapped_SNPs"
OUT_DIR="$HOME_DIR/significant_SNPs/similarity_scores"

# Ensure the final output directory exists safely before the jobs run wild
mkdir -p "$OUT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID starting..."

cd "$HOME_DIR"

"$PYTHON_EXEC" ./scripts/get_sim_scores.py $SLURM_ARRAY_TASK_ID "$OUT_DIR"
