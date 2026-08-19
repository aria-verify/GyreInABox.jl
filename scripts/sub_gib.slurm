#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=oceananigans
#SBATCH --time=02:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --qos=serial

#SBATCH --account=n02-VERIFY
#SBATCH --partition=serial

# Setup the job environment (this module needs to be loaded before any other modules)
module load PrgEnv-cray
module load cray-mpich

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Ensure the cpus-per-task option is propagated to srun commands
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

TEMPERATURE_RESTORING=20.0
SURFACE_EVAPORATION=-9.25e-8
#SIMULATION_YEARS=.2
OUTPUT_INTERVAL_DAYS=30
PLOT_START_TIME_DAYS=1800  # Needs to be a multiple of OUTPUT_INTERVAL_DAYS
PLOT_END_TIME_DAYS=$((SIMULATION_YEARS * 365))
#PLOT_END_TIME_DAYS=2
OUTPUT_FORMAT=NetCDF
#GRID="40 60 40"
#GRID="10 10 10"

export WORK="/work/n02/n02/ryapat30"
PROJECT_DIR=$WORK/Oceananigans_clean/projects/GyreInABox.jl
EXPERIMENT_DIR=$PROJECT_DIR/$SLURM_JOB_NAME/job-$SLURM_JOB_ID
mkdir -p $EXPERIMENT_DIR

export JULIA_DEPOT_PATH="$WORK/.julia"
export PATH="$PATH:$WORK/julia-1.12.6/bin"
export PATH="$PATH:$JULIA_DEPOT_PATH/bin"
# Define some paths
#export WORK=/work/n02/n02/ryapat30
#
#export JULIA="$WORK/julia-1.12.6/bin/julia"  # The julia executable
#export PATH="$PATH:$WORK/julia-1.12.6/bin"  # The folder of the julia executable
#export JULIA_DEPOT_PATH="$WORK/.julia"
export MPIEXECJL="$JULIA_DEPOT_PATH/bin/mpiexecjl"  # The path to the mpiexexjl executable

julia --project=$WORK/Oceananigans_clean/projects/GyreInABox.jl/scripts \
        	$WORK/Oceananigans_clean/projects/GyreInABox.jl/scripts/spall_dg.jl\
		-O $EXPERIMENT_DIR \
	       	-R 20.\
	       	-E -9.25e-8 \
	       	-G 400 600 40 \
	       	-Y 1 \
                -F NetCDF \
	       	--cpu
