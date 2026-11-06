#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --exclusive
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH -N 1
#SBATCH -J "aln_kp18"
#SBATCH -p largemem
#SBATCH --account=hpc-prf-rqsued
#SBATCH -D /scratch/hpc-prf-metdyn/metdyn07_Ritaj/benchmark_2d/aln/pr_kp18
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=true
ulimit -s unlimited
module load chem/CP2K/2026.2-foss-2025b-cpu
export LD_LIBRARY_PATH=/pc2/groups/hpc-prf-metdyn/metdyn07_Ritaj/implementation/github/cp2k-noctua/libwignernj-0.8.0/lib64:$LD_LIBRARY_PATH
export CP2K_DATA_DIR=/scratch/hpc-prf-metdyn/metdyn07_Ritaj/periodic_3d_test/cp2k-pr6447/data
srun /scratch/hpc-prf-metdyn/metdyn07_Ritaj/periodic_3d_test/cp2k-pr6447/build/bin/cp2k.psmp aln_kp18.inp &> aln_kp18.out
