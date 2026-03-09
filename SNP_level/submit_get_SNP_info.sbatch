#!/bin/bash
#SBATCH --job-name=sumStats_process
#SBATCH --time=01:00:00
#SBATCH --array=1-1060%60
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=28G               
#SBATCH --output=logs/sumstats_process/job_%A_%a.out
#SBATCH --error=logs/sumstats_process/job_%A_%a.err

PYTHON_EXEC="$HOME/.conda/envs/main-py3/bin/python"

# Set up paths
HOME_DIR="/home/sercan/dominance/overlapped_SNPs"
OUT_DIR="$HOME_DIR/significant_SNPs/sig_sumstats_merged"

# Ensure the final output directory exists safely before the jobs run wild
mkdir -p "$OUT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID starting..."

cd "$HOME_DIR"

p_threshold=$(awk 'BEGIN {print 5e-8 / 1060}') 

# Pass the $TMPDIR and the final output directory to your Python script
"$PYTHON_EXEC" ./scripts/get_SNP_info.py $SLURM_ARRAY_TASK_ID "$TMPDIR" "$OUT_DIR" "$p_threshold"
