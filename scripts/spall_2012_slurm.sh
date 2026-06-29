#!/bin/bash

#SBATCH --job-name=spall-2012
#SBATCH --output=%x-%j.out
#SBATCH --gpus=1
#SBATCH --ntasks-per-gpu=1
#SBATCH --time=3:00:00

# Run Spall (2012) model configuration and post process outputs

# Set up environment variables - project directory variable may need to be
# adjusted to wherever GyreInABox.jl repository is cloned to

TEMPERATURE_RESTORING=20.0
SURFACE_EVAPORATION=-9.25e-8
SIMULATION_YEARS=10
OUTPUT_INTERVAL_DAYS=30
PLOT_START_TIME_DAYS=1800  # Needs to be a multiple of OUTPUT_INTERVAL_DAYS
PLOT_END_TIME_DAYS=$((SIMULATION_YEARS * 365))

PROJECT_DIR=$HOME/projects/GyreInABox.jl
EXPERIMENT_DIR=$SCRATCH/$SLURM_JOB_NAME/job-$SLURM_JOB_ID
SCRIPTS_DIR=$PROJECT_DIR/scripts
RUN_SCRIPT_PATH=$SCRIPTS_DIR/spall_2012.jl
POST_PROCESS_SCRIPT_PATH=$SCRIPTS_DIR/spall_2012_post_process.jl

mkdir -p $EXPERIMENT_DIR

echo "Writing outputs to $EXPERIMENT_DIR..."

# Copy run and post process scripts to experiment directory

cp $RUN_SCRIPT_PATH $EXPERIMENT_DIR
cp $POST_PROCESS_SCRIPT_PATH $EXPERIMENT_DIR

# Run simulation outputting to experiment directory

echo "Running simulation for $SIMULATION_YEARS years..."

julia --project=$SCRIPTS_DIR $RUN_SCRIPT_PATH \
    -R $TEMPERATURE_RESTORING -E $SURFACE_EVAPORATION \
    -Y $SIMULATION_YEARS -I $OUTPUT_INTERVAL_DAYS -O $EXPERIMENT_DIR

# Post process simulation outputs in experiment directory

echo "Post processing simulation outputs..."

julia --project=$SCRIPTS_DIR $POST_PROCESS_SCRIPT_PATH $EXPERIMENT_DIR \
    -t $PLOT_START_TIME_DAYS $OUTPUT_INTERVAL_DAYS $PLOT_END_TIME_DAYS

