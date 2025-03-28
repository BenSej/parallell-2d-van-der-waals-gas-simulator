#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:ppn=4
#PBS -l naccesspolicy=shared
#PBS -N mp3-regression
#PBS -j oe
#PBS -q secondary
#PBS -S /projects/eng/shared/cs484/sing_shell.sh

#TODO: create a nodefile and populate PBS_NUM_NODES/PBS_NUM_PPN if not running in torque.

## If not started with PBS, figure out where we are relative to the build directory
#####Snippet from:   http://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#####end snippet
#IF PBS_O_WORKDIR is not set, we are not running in PBS, choose directory relative to script.
PBS_O_WORKDIR=${PBS_O_WORKDIR:-${SLURM_SUBMIT_DIR}}
PBS_O_WORKDIR=${PBS_O_WORKDIR:-${SCRIPT_DIR}/..}

#moves to the directory the user was in when they ran qsub
cd ${PBS_O_WORKDIR} #assumed to be the source tree

#check that the script was submit from the right place.
if [ -d "./cmake" ] && [ -d "./tests" ] && [ -d "./writeup" ]
then
	echo "We seem to be in the right place."
else
	echo "Not submit from the right place! Submit from the root of your repository."
	exit 1
fi

#creates an out-of-tree build directory for CMake and moves to it
mkdir -p ${PBS_O_WORKDIR}/build
pushd ${PBS_O_WORKDIR}/build

#build the programs (into the build directory, IE, the current directory)
#then benchmark them. Quit early on failure.
echo "Compiling (if necessary)"
cmake ${PBS_O_WORKDIR} && make part3 xyz_sorter || exit 1
popd

bash ${PBS_O_WORKDIR}/tests/p3regress.bash ${PBS_O_WORKDIR}/build > ${PBS_O_WORKDIR}/writeup/regression.txt
exit $?
