#!/bin/bash
#SBATCH --job-name=check_mult
#SBATCH --time=00:30:00 
#SBATCH --array=1-75%25         
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=28G
#SBATCH --output=logs/check_mult/job_%A_%a.out
#SBATCH --error=logs/check_mult/job_%A_%a.err

PYTHON_EXEC="$HOME/.conda/envs/main-py3/bin/python"

HOME_DIR="/home/sercan/dominance/dominance_pleiotropy"

SUMSTATS_DIR="$HOME_DIR/loci_level/sumstats_QCed"
PHEN_DICT_PATH="$HOME_DIR/UKB_sumstats_Neale/phen_dict.xlsx"
LOCI_DIR="$HOME_DIR/loci_level/sig_loci"
LD_DIR="$HOME_DIR/loci_level/ld_files"
OUT_DIR="$HOME_DIR/loci_level/mult_results"

echo "Task $SLURM_ARRAY_TASK_ID starting..."

cd "$HOME_DIR"

# Pass the $TMPDIR and the final output directory to your Python script
"$PYTHON_EXEC" ./scripts/loci_level/check_mult.py "$SLURM_ARRAY_TASK_ID" "$SUMSTATS_DIR" "$PHEN_DICT_PATH" "$LOCI_DIR" "$LD_DIR" "$OUT_DIR"