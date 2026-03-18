#!/bin/bash
#SBATCH --job-name=get_loci
#SBATCH --time=01:00:00
#SBATCH --array=1-85%85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=28G               
#SBATCH --output=logs/get_loci/job_%A_%a.out
#SBATCH --error=logs/get_loci/job_%A_%a.err

PYTHON_EXEC="$HOME/.conda/envs/main-py3/bin/python"

# Set up paths
HOME_DIR="/home/sercan/dominance/dominance_pleiotropy"
LD_DIR="$HOME_DIR/loci_level/ld_files"
OUT_DIR="$HOME_DIR/loci_level/sig_loci"

# Ensure the final output directory exists
mkdir -p "$OUT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID starting..."

cd "$HOME_DIR"

# Pass the $TMPDIR and the final output directory to your Python script
"$PYTHON_EXEC" ./scripts/loci_level/get_sumStats_QC.py "$SUMSTATS_PATH" "$LD_DIR" "$OUT_DIR"
