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
export WORKDIR="/tudelft.net/staff-umbrella/FlowDecomposition/vg-flow"
export SCRIPT_DIR="$WORKDIR/files_vg_flow/scripts"

# Set Gurobi license file path
export GRB_LICENSE_FILE="/tudelft.net/staff-umbrella/FlowDecomposition"

# Run the first script inside the container
srun apptainer exec \
  -B "$WORKDIR:$WORKDIR" \
  -B "$HOME:$HOME" \
  "$APPTAINER_ROOT/$APPTAINER_NAME" \
  python "$SCRIPT_DIR/build_graph_msga.py" -f $WORKDIR/files_vg_flow/example/forward.fastq -r $WORKDIR/files_vg_flow/example/reverse.fastq -c $WORKDIR/files_vg_flow/example/input.fasta -vg $WORKDIR//vg -t 8

# Check if the first script ran successfully
if [ $? -eq 0 ]; then
  echo "First script finished successfully. Running vg-flow.py..."

  # Run the second script inside the container
  srun apptainer exec \
    -B "$WORKDIR:$WORKDIR" \
    -B "$HOME:$HOME" \
    "$APPTAINER_ROOT/$APPTAINER_NAME" \
    python "$SCRIPT_DIR/vg-flow.py" -m 1 -c 2 node_abundance.txt contig_graph.final.gfa
else
  echo "First script failed. Aborting second script."
  exit 1
fi