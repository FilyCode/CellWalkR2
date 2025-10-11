#!/bin/bash -l

# --------------------------------------------------------------------------------
# SGE Job Directives
# --------------------------------------------------------------------------------
# Project name
#$ -P agedisease

# Hard time limit (hh:mm:ss). Adjust as needed.
#$ -l h_rt=48:00:00

# Job name
#$ -N RAgeSignatureCalculation

# Merge stdout and stderr into a single file (.o#jobID).
#$ -j y

# Email notification: job ends (e)
#$ -m e

# Your email address for notifications
#$ -M phitro@bu.edu

# Request cores.
#$ -pe omp 8

# Request total memory for the job.
#$ -l mem_total=8G

# --------------------------------------------------------------------------------
# Setup R Environment
# --------------------------------------------------------------------------------
echo "Loading R module..."
module load R/4.4.0 

# --------------------------------------------------------------------------------
# Prepare Working Directory
# --------------------------------------------------------------------------------
echo "Changing to submission directory: $SGE_O_WORKDIR"
cd $SGE_O_WORKDIR

# --------------------------------------------------------------------------------
# Run R Script
# --------------------------------------------------------------------------------
echo "Job started on host: $(hostname)"
echo "Requested CPU cores (NSLOTS): $NSLOTS" 

# Execute your R script using Rscript
Rscript 3_get_signatures_for_each_tissue_sapiens.R

# --------------------------------------------------------------------------------
# Post-job Cleanup
# --------------------------------------------------------------------------------
echo "R script finished."
