#!/usr/bin/env bash
#SBATCH --job-name=l1ou_array
#SBATCH --output=logs/l1ou_%A_%a.out
#SBATCH --error=logs/l1ou_%A_%a.err
#SBATCH --array=1-$(ls /path/to/rds_folder/*.rds | wc -l)
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# -------------------------------
# CONFIGURATION
# -------------------------------
DIR="/path/to/rds_folder"       # folder with RDS models
TREE="/path/to/tree_file.tre"   # tree file
NITER=100                       # number of bootstrap iterations
SCRIPT="/path/to/l1ou_bootstrap_Rscript.R" # the R script you posted above

# Create logs dir for stdout
mkdir -p logs

# -------------------------------
# SELECT RDS FILE FOR THIS TASK
# -------------------------------
RDS_FILE=$(ls ${DIR}/*.rds | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "Task ${SLURM_ARRAY_TASK_ID}: Processing ${RDS_FILE}"

# -------------------------------
# RUN
# -------------------------------
Rscript "${SCRIPT}" "${TREE}" "${RDS_FILE}" "${NITER}"
