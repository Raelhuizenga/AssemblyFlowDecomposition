#!/bin/sh
#SBATCH --partition=general   # Request partition
#SBATCH --qos=short           # Request Quality of Service
#SBATCH --time=2:00:00        # Request run time (wall-clock)
#SBATCH --ntasks=1            # Request number of parallel tasks per job
#SBATCH --cpus-per-task=8     # Request number of CPUs (threads) per task
#SBATCH --mem=2GB             # Request memory (MB) per node
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

export genome_size=1000
export num_haps=2

# Run the first script inside the container
srun apptainer exec \
  -B "$WORKDIR:$WORKDIR" \
  -B "$HOME:$HOME" \
  "$APPTAINER_ROOT/$APPTAINER_NAME" \
  python "$SCRIPT_DIR/simulate_data.py" -g $genome_size -m 0.005 -k $num_haps -f simulated_data/gridsearch -c 500

# Check if the first script ran successfully
if [ $? -eq 0 ]; then
  echo "First script finished successfully. Running the optimization script..."

  # Run the second script inside the container
  srun apptainer exec \
    -B "$WORKDIR:$WORKDIR" \
    -B "$HOME:$HOME" \
    "$APPTAINER_ROOT/$APPTAINER_NAME" \
    python "$SCRIPT_DIR/exact_flow_decomposition.py"  -t 0 -g graph_${genome_size}_${num_haps}.gfa -a abundances_${genome_size}_${num_haps}.txt -o output/simulated_data/gridsearch -i simulated_data/gridsearch
else
  echo "First script failed. Aborting second script."
  exit 1
fi