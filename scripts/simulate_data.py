import sys
import random
from graph_tool.all import Graph
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description="Simulate a genome variation graph")
    parser.add_argument('-g', '--genome_size', dest='genome_size', type=int, required=True, help="Size of the genome")
    parser.add_argument('-m', '--mutation_rate', dest='mutation_rate', type=float, required=True, help="Mutation rate")
    parser.add_argument('-k', '--num_haps', dest='num_haps', type=int, required=True, help="Number of haplotypes, creates test instances from 2 to i")
    parser.add_argument('-f', '--file_location', dest='file_location', type=str, required=True, help="Location of the folder where the data should be stored in")
    parser.add_argument('-c', '--contig_length', dest='contig_length', type=int, default=100, help="Number of bases of a contig")
    args = parser.parse_args()

    seq = create_genome(args.genome_size)
    haps = create_haplotypes(seq, args.mutation_rate, args.num_haps)
    # weights = exponential_weights(args.num_haps)
    weights = random_exponential_weights(args.num_haps)
    g, contigs = make_graph(haps, weights, args.contig_length)
    g = contract_vertices(g)
    save_haps_to_file(haps, weights, f'{args.file_location}/true_haps_{args.genome_size}_{args.num_haps}.fasta')
    write_gfa(g, f'{args.file_location}/graph_{args.genome_size}_{args.num_haps}.gfa')
    create_abundance_file(g, f'{args.file_location}/_abundance_{args.genome_size}_{args.num_haps}.txt')


def exponential_weights(size_array) -> list:
    # create exponential weights
    return [2**i for i in range(size_array)]


def random_exponential_weights(size_array) -> list:
    # create random weights
    return np.random.exponential(2**size_array, size_array).astype(int).tolist()


def create_abundance_file(graph: Graph, filename:str):
        with open(filename, 'w') as f:
            for v in graph.vertices():
                f.write(f"{int(v)}:{graph.vp.weight[v]}\n")


def create_haplotypes(genome: str, mutationrate: float, num_haps: int) -> np.array:
    haplotypes = []
    bases = "ACTG"

    for i in range(num_haps):
        haplotype = list(genome)
        for j in range(len(haplotype)):
            if random.random() < mutationrate:  # Mutate with probability mutationrate
                haplotype[j] = random.choice([b for b in bases if b != genome[i]])  # Choose a different base
        haplotypes.append("".join(haplotype))
    return np.array(haplotypes)


def create_genome(length: int) -> str:
    return ''.join(random.choices("ACTG", k=length))


def make_graph(haplotypes: np.ndarray, weights: np.ndarray, contig_length: int) -> tuple[Graph, list]:
    """
    Returns graph in graph-tool format.
    """
    g = Graph(directed=True) 
    g.vp.seq = g.new_vertex_property('string')
    g.vp.contigs = g.new_vertex_property('vector<string>')
    g.vp.weight = g.new_vertex_property('int', val=0)
    g.ep.ori = g.new_edge_property('string')


    last_vertices = []
    contig_dict = {}
    for i in range(len(haplotypes[0])):
        used_bases = []
        corresponding_vertices = []
        current_vertices = []
        for j in range(len(haplotypes)):
            if haplotypes[j][i] in used_bases: # node already exists
                index = used_bases.index(haplotypes[j][i])
                vertex = corresponding_vertices[index]
                g.vp.weight[vertex] += weights[j]
            else: # create new node
                vertex = g.add_vertex()
                g.vp.seq[vertex] = haplotypes[j][i]
                g.vp.weight[vertex] = weights[j]
                g.vp.contigs[vertex] = []
                used_bases.append(haplotypes[j][i])
                corresponding_vertices.append(vertex)
            if i != 0: # create edges
                if not g.edge(last_vertices[j], vertex):
                    e = g.add_edge(last_vertices[j], vertex)
                    g.ep.ori[e] = '++'  # Orientation of edge
            if i % contig_length == 0:
                contig = haplotypes[j][i::]
                if len(contig) > contig_length + 20:
                    contig = contig[:contig_length + 19]
                if vertex in contig_dict:
                    contig_dict[vertex].append(contig)
                else:
                    contig_dict[vertex] = [contig]
            current_vertices.append(vertex)
        last_vertices = current_vertices

    # add contigs to vertices
    contig_number = 0
    contig_numbers_list = []
    for v, contigs in contig_dict.items():
        for seq in contigs:
            assert g.vp.seq[v] == seq[0]
            g.vp.contigs[v].append(str(contig_number))
            seq = seq[1:]
            next_node = v
            while len(seq) > 0:
                new_node = [u  for _, u in next_node.out_edges() if g.vp.seq[u] == seq[0]][0]
                g.vp.contigs[new_node].append(str(contig_number))
                seq = seq[1:]
                next_node = new_node
            contig_numbers_list.append(str(contig_number))
            contig_number += 1
    return g, contig_numbers_list


def contract_vertices(g: Graph) -> Graph:
    """
    Contract vertices that form a simple path.
    """
    edges_to_contract = list(g.edges())
    vertices_to_remove = []
    for s,t in edges_to_contract:
        if (s.out_degree() == 1 and t.in_degree() == 1) and g.vp.weight[s] == g.vp.weight[t]:
            g.vp.seq[t] = g.vp.seq[s] + g.vp.seq[t]
            g.vp.contigs[t] = list(set(g.vp.contigs[s]) | set(g.vp.contigs[t]))
            for e in s.in_edges():
                new_start_vertex = e.source()
                e =  g.add_edge(new_start_vertex, t)
                g.ep.ori[e] = '++'  # Orientation of edge
            vertices_to_remove.append(s)
    g.remove_vertex(reversed(sorted(vertices_to_remove)), fast=False)
    g = g.copy()  # Create a copy to reindex vertices
    return g


def write_gfa(g, gfa_file, paths=None):
    """
    Write graph-tool graph to a gfa format storing nodes, edges, and paths.
    Returns nothing.
    """
    contigs = {} # dict mapping contigs to lists of nodes
    with open(gfa_file, 'w') as f:
        f.write("H\tVN:Z:1.0") # header line GFA 1.0
        for v in g.vertices():
            f.write("\nS\t{}\t{}".format(int(v), g.vp.seq[v]))
            for w in v.out_neighbors():
                e = g.edge(v, w)
                ori = g.ep.ori[e]
                assert len(ori) == 2
                f.write("\nL\t{}\t{}\t{}\t{}\t0M".format(int(v), ori[0],
                                                         int(w), ori[1]))
            for k in g.vp.contigs[v]:
                if k in contigs.keys():
                    contigs[k].append(v)
                else:
                    contigs[k] = [v]

        if paths != None:
            if len(paths) > 0:
                # write paths from path list
                path_info = paths
            else:
                # write paths based on graph traversal
                path_info = contigs
            for k, node_list in path_info.items():
                if len(node_list) == 0:
                    print("Skipping length 0 path")
                    continue
                path = ""
                # for v in contigs[k]:
                for v in node_list:
                    path += str(int(v)) + '+,'
                path = path.rstrip(',')
                f.write("\nP\t{}\t{}".format(k, path))
    return

def save_haps_to_file(haps, weights, filename):
    # save haps to file
    with open(filename, "w") as f:
        for i in range(len(haps)):
            f.write(f'>path{i} {weights[i]}x freq={round(weights[i] / sum(weights), 3)}\n{haps[i]}\n')


if __name__ == '__main__':
    sys.exit(main())