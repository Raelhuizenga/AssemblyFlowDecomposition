#!/bin/sh
#SBATCH --partition=general   # Request partition
#SBATCH --qos=short           # Request Quality of Service
#SBATCH --time=0:05:00        # Request run time (wall-clock)
#SBATCH --ntasks=1            # Request number of parallel tasks per job
#SBATCH --cpus-per-task=2     # Request number of CPUs (threads) per task
#SBATCH --mem=1GB             # Request memory (MB) per node
#SBATCH --mail-type=END       # Notify when the job ends
#SBATCH --output=slurm_%j.out # Set name of output log
#SBATCH --error=slurm_%j.err  # Set name of error log

# Define container folder and name
export APPTAINER_ROOT="."
export APPTAINER_NAME="vg-flow-image.sif"

# Get current working directory
export WORKDIR="/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition"
export SCRIPT_DIR="$WORKDIR/scripts"

# Set Gurobi license file path
export GRB_LICENSE_FILE="/tudelft.net/staff-umbrella/FlowDecomposition"

# Run the first script inside the container
srun apptainer exec \
  -B "$WORKDIR:$WORKDIR" \
  -B "$HOME:$HOME" \
  "$APPTAINER_ROOT/$APPTAINER_NAME" \
  python "$SCRIPT_DIR/exact_flow_decomposition.py"