#!/bin/bash
#SBATCH --job-name=mw14
#SBATCH --partition=<partition>
#SBATCH --gres=gpu
#SBATCH --nodelist=<node>
#SBATCH --exclusive
#SBATCH --time=4-04:00:00
#SBATCH --output=./log/run_mw2014.out
#SBATCH --mail-type=TIME_LIMIT,END,FAIL
#SBATCH --mail-user=<your email>

# Use available modules
module purge
module load gcc/9.4.0
module load cuda/11.7
module load hdf5/1.10.7-intel
module load inteloneapi/2021.2


NBODY=/path/to/install/bin/nbody6++.avx.gpu.hdf5
INPUT=/path/to/input/N_10000_MW2014.input
PREFIX=/path/to/output/folder/
IC=/path/to/ic/dat.10

mkdir -p ${PREFIX}
cd ${PREFIX}

cp ${IC} ${PREFIX}dat.10

ulimit -s unlimited
export OMP_NUM_THREADS=32

time ${NBODY} < ${INPUT} > OUT