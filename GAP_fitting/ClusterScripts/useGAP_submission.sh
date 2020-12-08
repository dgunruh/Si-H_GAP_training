#!/bin/bash
#
#! -cwd
#! -j y
#! -S /bin/bash

# Name of the job
#SBATCH --job-name=QUIP
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH --partition=high2                 # Use the high partition
#SBATCH -t 1-00:00                      # Runtime in D-HH:MM format
#SBATCH -o 'farm_outputs/gap_fitting-%j.output' #File to which STDOUT will be written
#SBATCH --mail-user="dgunruh@ucdavis.edu"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# run one thread for each one the user asks the queue for
# hostname is just for debugging
# hostname
export OMP_NUM_THREADS=64
# export t=$SLURM_ARRAY_TASK_ID
export j=$SLURM_JOB_ID

module load openmpi

# Folders for data: checkGAPvsDFT_postMDData/Iterative_Training_round2.xyz; trainingGAP_iterationData/round2Iteration.xyz
# Folders in quip_results: checkGAPvsDFT_afterMD (NA); checkGAPvsDFT_afterTraining_initialSet (recheck); checkGAPvsDFT_afterTraining_fullSet (crosscheck)

srun ~/QUIP/build/linux_x86_64_gfortran_openmp/quip E=T F=T V=T atoms_filename=checkGAPvsDFT_postMDData/Iterative_Training_round7.xyz param_filename=potentials/GAP_glue_soap_round6/GAP_glue_soap_round6.xml | grep AT | sed 's/AT//' | sed -e 's/^[ \t]*//' > quip_results/checkGAPvsDFT_afterMD/quip_iteration7.xyz

# sed -e 's/^[ \t]*//'
