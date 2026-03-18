#!/bin/bash
#SBATCH --job-name=ld_coverage
#SBATCH --time=00:10:00 
#SBATCH --array=1-1%1         
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=28G
#SBATCH --output=logs/ld_coverage/job_%A_%a.out
#SBATCH --error=logs/ld_coverage/job_%A_%a.err

PYTHON_EXEC="$HOME/.conda/envs/main-py3/bin/python"

HOME_DIR="/home/sercan/dominance/dominance_pleiotropy"
SUMSTATS_PATH="$HOME_DIR/loci_level/sumstats_QCed"
LD_DIR="$HOME_DIR/loci_level/ld_files"
ALL_SIG_FILE="$HOME_DIR/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
    
OUT_DIR="$HOME_DIR/loci_level"

echo "Task $SLURM_ARRAY_TASK_ID starting..."

cd "$HOME_DIR"

# Pass the $TMPDIR and the final output directory to your Python script
"$PYTHON_EXEC" ./scripts/loci_level/coverage.py "$SUMSTATS_PATH" "$LD_DIR" "$OUT_DIR" "$ALL_SIG_FILE"