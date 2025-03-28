#!/bin/bash
#SBATCH --time=01:20:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --job-name cs484-mp3-experiments
#SBATCH -p eng-instruction

export MAXIMUM_DECOMPOSITION=256
export USE_DECOMPOSITION=100
export TOTAL_CPUS=$(( ${SLURM_JOB_NUM_NODES} * ${SLURM_CPUS_ON_NODE} ))
export PARTICLE_RUNSPEC="-s 48 --high-density 0.3 --low-density 0.1 --exchange 5 --steps 10000 -i 500"
export CHARM_BALANCER_OPTIONS="+balancer GreedyRefine +isomalloc_sync --process-migrate 1000"

## If not started with PBS, figure out where we are relative to the build directory
#####Snippet from:   http://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#####end snippet
#IF SLURM_SUBMIT_DIR is not set, we are not running in PBS, choose directory relative to script.
SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}/..}

#moves to the directory the user was in when they ran qsub
cd ${SLURM_SUBMIT_DIR} #assumed to be the source tree

#check that the script was submit from the right place.
if [ -d "./cmake" ] && [ -d "./tests" ] && [ -d "./writeup" ]
then
	echo "We seem to be in the right place."
else
	echo "Not submit from the right place! Submit from the root of your repository."
	exit 1
fi

set -e

#creates an out-of-tree build directory for CMake and moves to it
mkdir -p ${SLURM_SUBMIT_DIR}/build
pushd ${SLURM_SUBMIT_DIR}/build

#build the programs (into the build directory, IE, the current directory)
#then benchmark them. Quit early on failure.
echo "Compiling"
srun --mpi=pmi2 --ntasks 1 /projects/eng/shared/cs484/sing_exec.sh cmake ${SLURM_SUBMIT_DIR} && srun --mpi=pmi2 --ntasks 1 /projects/eng/shared/cs484/sing_exec.sh make

#Google test tests
echo "Testing"
srun --mpi=pmi2 --ntasks 4 /projects/eng/shared/cs484/sing_exec.sh ./bin/run_tests || ( echo "No use benchmarking an incorrect program." ; exit 1 )



OUTFILE_MPI_BENCHMARK="${SLURM_SUBMIT_DIR}/writeup/benchmark_mpi-${SLURM_JOB_ID}.txt"
OUTFILE_AMPI_BENCHMARK="${SLURM_SUBMIT_DIR}/writeup/benchmark_ampi-${SLURM_JOB_ID}.txt"
OUTFILE_CHARM_BENCHMARK="${SLURM_SUBMIT_DIR}/writeup/benchmark_charm-${SLURM_JOB_ID}.txt"


echo "BEGIN _VARIES_"             >> ${OUTFILE_MPI_BENCHMARK}
echo "BEGIN ${USE_DECOMPOSITION}" >> ${OUTFILE_AMPI_BENCHMARK}
echo "BEGIN ${USE_DECOMPOSITION}" >> ${OUTFILE_CHARM_BENCHMARK}

for N_CPUS in 1 4 9 16 25 36
do

	#Individual MPI run
	echo "Doing MPI CPUS = ${N_CPUS} "
	echo "Doing MPI CPUS = ${N_CPUS} " >> ${OUTFILE_MPI_BENCHMARK}
	srun --mpi=pmi2 --ntasks-per-node ${SLURM_NTASKS_PER_NODE} --cpu-bind=cores --ntasks ${N_CPUS} /projects/eng/shared/cs484/sing_exec.sh ./bin/part1 ${PARTICLE_RUNSPEC} | grep "Simulation_Duration" >> ${OUTFILE_MPI_BENCHMARK}

	echo "Doing aMPI CPUS = ${N_CPUS} "
	echo "Doing aMPI CPUS = ${N_CPUS} " >> ${OUTFILE_AMPI_BENCHMARK}
	srun --mpi=pmi2 --ntasks-per-node ${SLURM_NTASKS_PER_NODE} --cpu-bind=cores --ntasks ${N_CPUS} /projects/eng/shared/cs484/sing_exec.sh ./bin/part2 ${PARTICLE_RUNSPEC} +vp ${USE_DECOMPOSITION} ${CHARM_BALANCER_OPTIONS} | grep "Simulation_Duration" >> ${OUTFILE_AMPI_BENCHMARK}

#	echo "Doing CHARM CPUS = ${N_CPUS} "
#	echo "Doing CHARM CPUS = ${N_CPUS} " >> ${OUTFILE_CHARM_BENCHMARK}
#	srun --mpi=pmi2 --ntasks-per-node ${SLURM_NTASKS_PER_NODE} --cpu-bind=cores --ntasks ${N_CPUS} /projects/eng/shared/cs484/sing_exec.sh ./bin/part3 ${PARTICLE_RUNSPEC} --ranks ${USE_DECOMPOSITION} ${CHARM_BALANCER_OPTIONS} | grep "Simulation_Duration" >> ${OUTFILE_CHARM_BENCHMARK}

done
