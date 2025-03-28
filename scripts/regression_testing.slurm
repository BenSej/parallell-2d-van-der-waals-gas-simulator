#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name mp3-regression
#SBATCH -p secondary
#PBS -S /projects/eng/shared/cs484/sing_shell.sh

#TODO: create a nodefile and populate PBS_NUM_NODES/PBS_NUM_PPN if not running in torque.

## If not started with PBS, figure out where we are relative to the build directory
#####Snippet from:   http://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#####end snippet
#IF SLURM_SUBMIT_DIR is not set, we are not running in SLURM, choose directory relative to script.
SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}/..}

#moves to the directory the user was in when they ran qsub
cd ${SLURM_SUBMIT_DIR} #assumed to be the source tree

export PBS_O_WORKDIR=${SLURM_SUBMIT_DIR}
srun -n 1 /projects/eng/shared/cs484/sing_exec.sh ${SLURM_SUBMIT_DIR}/scripts/regression_testing.bash
