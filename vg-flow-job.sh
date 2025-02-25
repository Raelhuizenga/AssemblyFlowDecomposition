#!/bin/sh
#SBATCH --partition=general   # Request partition
#SBATCH --qos=short           # Request Quality of Service
#SBATCH --time=2:00:00        # Request run time (wall-clock)
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
export WORKDIR="/tudelft.net/staff-umbrella/FlowDecomposition"
export SCRIPT_DIR="$WORKDIR/flow_assembly/AssemblyFlowDecomposition/scripts"

# Set Gurobi license file path
export GRB_LICENSE_FILE="$WORKDIR/gurobi.lic"

# Run the first script inside the container
srun apptainer exec \
  -B "$WORKDIR:$WORKDIR" \
  -B "$HOME:$HOME" \
  "$APPTAINER_ROOT/$APPTAINER_NAME" \
  python "$SCRIPT_DIR/vg-flow.py" -m 1 -c 1 -d 0 "output/vg-flow/data/abundance_10001_4.txt" "output/vg-flow/data/graph_10001_4.gfa" --max_strains 10 --trim 0

if [ $? -eq 0 ]; then
  echo "Script finished successfully."
fi