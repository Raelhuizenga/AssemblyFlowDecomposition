# VG-flow: haplotype reconstruction and abundance estimation from mixed samples using flow variation graphs.

License: MIT (see LICENSE)

Current version: 0.0.4

### Motivation ###
The goal of haplotype-aware genome assembly is to reconstruct all individual
haplotypes from a mixed sample and to provide corresponding abundance estimates.
VG-flow provides a reference-genome-independent solution based on the
construction of a flow variation graph, capturing all quasispecies diversity
present in the sample. VG-Flow computes a solution to the contig abundance
estimation problem and employs a greedy algorithm to efficiently build full-length
haplotypes. Finally, VG-Flow calculates accurate frequency estimates for the
reconstructed haplotypes through linear programming techniques.

### Citation ###
The algorithms used by VG-flow are described in detail in our [preprint](https://www.biorxiv.org/content/10.1101/645721v3).
Please cite this preprint if you are using VG-flow.


### Installation and dependencies ###
[Download](https://bitbucket.org/jbaaijens/vg-flow/downloads/?tab=tags) the latest
release of this VG-flow repository. Let `/path/to/vg-flow/` denote the path to the
vg-flow root directory.

The **easiest and recommended** way to install all dependencies required to run
VG-Flow is using [(mini)conda](https://conda.io/miniconda.html):
```
conda create --name vg-flow-env
conda activate vg-flow-env
conda install -c bioconda -c conda-forge -c gurobi python=3 graph-tool minimap2 gurobi biopython numpy
```
You now have an active conda environment with all dependencies installed.
In addition, you need to make sure that you
have a [license for Gurobi](https://user.gurobi.com/download/licenses/free-academic) (free for academic use).
Note that the Gurobi installation is already managed by conda.

If you prefer not to use conda, please make sure to install the following
dependencies:
- [vg](https://github.com/vgteam/vg)
- [Gurobi](https://www.gurobi.com) + [free academic license](https://user.gurobi.com/download/licenses/free-academic)
- [minimap2]()
as well as the following python modules:
- graph-tool
- biopython
- numpy


### Input ###
VG-Flow requires as input a fasta file with pre-assembled contigs and two fastq
files (forward and reverse) with the original sequencing data. For the assembly
step, we recommend using [savage](https://bitbucket.org/jbaaijens/savage)
because it produces highly accurate, strain-specific contigs.

### Output ###
VG-Flow reconstructs the haplotypes present in the sample, along with their
relative abundance rates. It outputs a fasta file with the resulting haplotypes,
the corresponding frequencies added to the identifier field. In addition, the
final genome variation graph is written to a GFA file.

### Workflow ###
If not yet active, make sure to activate your conda virtual environment:
```
conda activate vg-flow-env
```

**Step 1: variation graph construction**

Using `build_graph_msga.py` you build your contig variation graph: first it performs
multiple sequence alignment (MSA) using the [vg toolkit](https://github.com/vgteam/vg).
The resulting graph is then transformed into a contig variation graph and node
abundances are computed by mapping the original read set to the graph.
```
python /path/to/vg-flow/scripts/build_graph_msga.py -f <forward_fastq> -r <reverse_fastq> -c <contig_fasta> -vg /path/to/vg/executable -t <num_threads>
```
The graph and the node abundances are stored in `contig_graph.final.gfa` and
`node_abundances.txt` and used as input for step 2.
Note: this is the most time-consuming step in the workflow, because it requires read
mapping to the contig variation graph. Allowing multithreading with `-t <num_threads>`
can speed-up this step significantly.

***Remark:*** VG-Flow also allows you to input a homemade variation graph using
the option --reuse_graph. Make sure to save your graph in vg format under the
name `contig_graph.norm.vg`. VG-Flow will use the vg toolkit to index your graph
and perform read alignment, followed by the rest of the VG-Flow workflow.

**Step 2: build haplotypes**

This is the most interesting and also the most challenging step. VG-flow solves a
linear program (derived from a min-cost flow problem) that yields contig abundance
estimates. Based on these abundances, a selection of candidate haplotypes is generated.

Next, we assign relative abundance rates to each of the candidate paths, while minimizing
the difference between node abundance and the sum of abundance rates of all
strains going through the node, summing over all nodes. This is formalized as
an optimization problem and solved using the LP solver in Gurobi.

This step is executed by running `vg-flow.py`. The user needs to specify the threshold
values `-m`, the minimal node abundance, and `-c`, the minimal strain abundance. We
recommend setting -m=0.005\*sequencing depth and -c=0.01\*sequencing depth. For example,
with an average sequencing depth of 20,000x, set -m=100 and -c=200.
```
python /path/to/vg-flow/scripts/vg-flow.py -m <node_ab> -c <strain_ab> node_abundance.txt contig_graph.final.gfa
```
The final haplotypes are written to `haps.final.fasta` and the genome variation graph to `genome_graph.gfa`.

The conda environment can be deactivated by typing `conda deactivate`.


### Example ###

The example directory contains a small example to test VG-Flow before applying it to
your own data. The example data consists of a small set of contigs `input.fasta`, and
two read files `forward.fastq` and `reverse.fastq`, corresponding to a sequencing depth
of 200x.

To run the example, enter the example directory, activate the conda environment, and
execute the following shell commands:
```
python ../scripts/build_graph_msga.py -f forward.fastq -r reverse.fastq -c input.fasta -vg /path/to/vg/executable -t 8
python ../scripts/vg-flow.py -m 1 -c 2 node_abundance.txt contig_graph.final.gfa
```
Note that we use `-m=1` and `-c=2` because of the total sequencing depth of 200x. Please
make sure to adjust these parameter settings to your own data as described in **step 2**.


### Contact ###

In case of any questions or issues, please contact Jasmijn Baaijens.
