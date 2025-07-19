#!/bin/bash
#SBATCH --job-name=mw14_save
#SBATCH --partition=<partition>
#SBATCH --nodelist=<node>
#SBATCH --exclusive
#SBATCH --time=2-00:00:00
#SBATCH --output=./log/save_mw2014.out
#SBATCH --mail-type=TIME_LIMIT,END,FAIL
#SBATCH --mail-user=<your email>

SNAP_INIT=0
SNAP_LAST=100

module purge
PYTHON=/path/to/python

PREFIX=/path/to/output/folder/
SAVE=/path/to/scripts/ReadNbodyOutput.py

cd ${PREFIX}
mkdir -p output

${PYTHON} ${SAVE} -snap_init ${SNAP_INIT} -snap_last ${SNAP_LAST}