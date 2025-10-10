#!/bin/bash -l

# --------------------------------------------------------------------------------
# SGE Job Directives
# --------------------------------------------------------------------------------
# Project name
#$ -P agedisease

# Hard time limit (hh:mm:ss). Adjust as needed.
#$ -l h_rt=12:00:00

# Job name
#$ -N RCensusDataFetch

# Merge stdout and stderr into a single file (.o#jobID).
#$ -j y

# Email notification: job ends (e)
#$ -m e

# Your email address for notifications
#$ -M phitro@bu.edu

# Request cores.
#$ -pe omp 2

# Request total memory for the job.
#$ -l mem_total=60G

# --------------------------------------------------------------------------------
# Setup R Environment
# --------------------------------------------------------------------------------
echo "Loading R module..."
module load R/4.4.0 

# If you installed your R packages in a non-standard location and R can't find them,
# you might need to uncomment and set this environment variable:
# export R_LIBS_USER="/path/to/your/custom/R/library"

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
Rscript 2_download_data_tabula_sapiens.R

# --------------------------------------------------------------------------------
# Post-job Cleanup
# --------------------------------------------------------------------------------
echo "R script finished."
