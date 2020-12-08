#!/bin/bash
#
#! -cwd
#! -j y
#! -S /bin/bash

# Name of the job
#SBATCH --job-name=GAPe0
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH --partition=high2                 # Use the high partition
#SBATCH -t 2-00:00                      # Runtime in D-HH:MM format
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

set lc="{"
set rc="}"

srun ~/QUIP/build/linux_x86_64_gfortran_openmp/gap_fit \
gap={soap Z=14 \
	n_species=2 \
	species_Z="${lc}1 14${rc}" \
	cutoff=5.0 \
	cutoff_transition_width=1.0 \
	atom_sigma=0.5 \
	l_max=6 \
	n_max=12 \
	zeta=4 \
	delta=3 \
	covariance_type=dot_product \
	sparse_method=cur_points \
	n_sparse=6000:\
     soap Z=1 \
	n_species=2 \
	species_Z="${lc}1 14${rc}" \
	cutoff=3.0 \
	cutoff_transition_width=0.5 \
	atom_sigma=0.3 \
	l_max=6 \
	n_max=12 \
	zeta=4 \
	delta=1 \
	covariance_type=dot_product \
	sparse_method=cur_points \
	n_sparse=3000} \
e0_method=isolated \
energy_parameter_name=dft_energy force_parameter_name=dft_force virial_parameter_name=dft_virial config_type_parameter_name=config_type \
do_copy_at_file=F \
sparse_separate_file=T \
at_file=trainingGAP_iterationData/round6Iteration.xyz \
default_sigma={0.001 0.1 0.05 0} \
config_type_sigma={isolated_atom:0.0001:0.01:0.0001:0.0:liquid:0.003:0.15:0.2:0.0:amorph:0.01:0.2:0.4:0.0} \
sparse_jitter=1.0e-8 \
core_param_file=potentials/pairPotentials/SiH_pairpot.xml \
core_ip_args={IP Glue} \
gp_file=GAP_glue_soap_round6.xml 2>&1 | grep -v FoX
