#!/bin/bash
#SBATCH --job-name=ld_download
#SBATCH --time=00:20:00 
#SBATCH --array=1-64%64          
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=28G
#SBATCH --output=logs/ld_down/job_%A_%a.out
#SBATCH --error=logs/ld_down/job_%A_%a.err


HOME_DIR="/home/sercan/dominance/dominance_pleiotropy"
LIST_PATH="$HOME_DIR/loci_level/required_LD_chunks.txt"
OUT_DIR="$HOME_DIR/loci_level/ld_files"

mkdir -p "$OUT_DIR"

# Extract the specific line from your specified path
BASENAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_PATH")

echo "Job ${SLURM_ARRAY_TASK_ID}: Downloading ${BASENAME} to ${OUT_DIR}"

# Download to the specified output directory
wget -q -nc "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/${BASENAME}.npz" -P "$OUT_DIR"
wget -q -nc "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/${BASENAME}.gz" -P "$OUT_DIR"

echo "Job ${SLURM_ARRAY_TASK_ID} finished."