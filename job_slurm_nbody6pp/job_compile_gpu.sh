#!/bin/bash
#SBATCH --job-name=nbody6++
#SBATCH --partition=<partition>
#SBATCH --nodelist=<node>
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --output=./log/test_nbody6pp.out
#SBATCH --mail-type=TIME_LIMIT,END,FAIL
#SBATCH --mail-user=<your email>

# Use available modules
module purge
module load gcc/9.4.0
module load cuda/11.7
module load hdf5/1.10.7-intel
module load inteloneapi/2021.2

PREFIX=/path/to/Nbody6ppGPU/
INSTALLPATH=/path/to/install/folder/

RUNCONFIG=./configure

echo "Starting..."

cd ${PREFIX}

echo "Configure NBODY6++"

${RUNCONFIG} --prefix=${INSTALLPATH} --enable-simd=avx --enable-gpu --enable-hdf5 --disable-mpi --with-nmax=1048576 --with-kmax=65536 --with-lmax=600 --with-mmax=1024 

echo "Create Makefile and install..."

make clean
make
make install

echo "Ending..."