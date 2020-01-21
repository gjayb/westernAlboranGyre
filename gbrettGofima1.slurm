#!/bin/sh
#SBATCH -p sched_mit_darwin
#SBATCH -N 4
#SBATCH -n 56
#SBATCH --mem-per-cpu 4000
#SBATCH -t 12:00:00
#SBATCH -J gbAF

module add harvard/centos6/netcdf-4.3.0_gcc-4.4.7
module add mvapich2/gcc/64/2.0b
module add harvard/centos6/hdf5-1.8.11_gcc-4.4.7
mpirun -np 56 ./mitgcmuvFull56 > output1.txt
