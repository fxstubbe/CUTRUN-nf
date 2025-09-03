#!/bin/bash
#SBATCH --job-name=CUTRUN_nextflow_pipeline
#SBATCH --output=logs/nextflow_%j.log
#SBATCH --error=logs/nextflow_%j.err
#SBATCH --partition=shared-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=francois-xavier.stubbe@unige.ch

# -----------------------------------
# Load required modules
# -----------------------------------
module purge
module load Java/17
module load Nextflow

# -----------------------------------
# Environment for Singularity / Apptainer binding
# -----------------------------------
# Bind scratch directory so the container can access workDir

#export APPTAINER_BIND="/home/users/s/stubbe/scratch"

# -----------------------------------
# Prepare directories for logs and reports
# -----------------------------------
#mkdir -p ../logs
mkdir -p ../reports

# -----------------------------------
# Run the pipeline
# -----------------------------------
nextflow run ./main.nf \
    -profile unige_baobab \
    -params-file ./params.yml \
    -with-trace ../reports/trace.txt \
    -with-report ../reports/report.html \
    -with-timeline ../reports/timeline.html \
    -with-dag ../reports/flowchart.png \

