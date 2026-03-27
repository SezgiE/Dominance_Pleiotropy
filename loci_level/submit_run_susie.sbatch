#!/bin/bash
#SBATCH --job-name=run_susie
#SBATCH --time=00:30:00
#SBATCH --array=1-75%25
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=28G               
#SBATCH --output=logs/run_susie/job_%A_%a.out
#SBATCH --error=logs/run_susie/job_%A_%a.err

# R env
module purge
module load 2025
module load R/4.5.1-gfbf-2025a

# Python env
PYTHON_EXEC="$HOME/.conda/envs/main-py3/bin/python"

# Set up paths
HOME_DIR="/home/sercan/dominance/dominance_pleiotropy"
SUMSTATS_DIR="$HOME_DIR/loci_level/sumstats_QCed"
PHEN_DICT_PATH="$HOME_DIR/UKB_sumstats_Neale/phen_dict.xlsx"
LOCI_DIR="$HOME_DIR/loci_level/sig_loci"
LD_DIR="$HOME_DIR/loci_level/ld_files"
R_SCRIPT_PATH="$HOME_DIR/scripts/loci_level/run_SuSie_core.R"
OUT_DIR="$HOME_DIR/loci_level/susie_results"

# Ensure the final output directory exists
mkdir -p "$OUT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID starting..."

cd "$HOME_DIR"

# Pass the $TMPDIR and the final output directory to your Python script
"$PYTHON_EXEC" ./scripts/loci_level/run_susie.py "$SLURM_ARRAY_TASK_ID" "$SUMSTATS_DIR" "$PHEN_DICT_PATH" "$LOCI_DIR" "$LD_DIR" "$R_SCRIPT_PATH" "$OUT_DIR"
