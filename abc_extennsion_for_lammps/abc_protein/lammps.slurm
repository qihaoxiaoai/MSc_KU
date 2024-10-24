#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1000:00:00
#SBATCH --job-name=lam_abc_1
#SBATCH --partition=X96C

module load gcc/5.4.0 fftw/3.3.8 openmpi/3.1.4-gnu-5.4.0

mpirun -np 4 lmp_mpi_abc -in in.gro2lam
