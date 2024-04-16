#!/bin/bash -l

#SBATCH --cluster=wice
#SBATCH --account=lp_dreesenlab
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=72
#SBATCH --time=00:59:00

module load NONMEM/7.5.0-GCC-10.3.0-MPICH-3.4.2
module load PsN/5.3.0-foss-2021a

cd $SLURM_SUBMIT_DIR
execute -nodes=2 -threads=72 -nm_version=nm750 -parafile=mpilinux8.pnm 1.mod