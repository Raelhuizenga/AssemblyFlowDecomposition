import sys
import random
from graph_tool.all import Graph
import numpy as np

def main():
    seq = create_genome(4)
    print(seq)
    haps = create_haplotypes(seq, 0.2, 3)
    print(haps)
    weights = [1, 2, 4]
    g, contigs = make_graph(haps, weights)
    print(contigs)
    g = contract_vertices(g)
    paths = get_contig_paths_from_graph(g, contigs)
    print(paths)
    print(str(g))


def create_haplotypes(genome: str, mutationrate: float, num_haps: int) -> np.array:
    haplotypes = [genome]
    bases = "ACTG"
    random.seed(43)

    for i in range(num_haps - 1):
        haplotype = list(genome)
        for j in range(len(haplotype)):
            if random.random() < mutationrate:  # Mutate with probability mutationrate
                haplotype[j] = random.choice([b for b in bases if b != genome[i]])  # Choose a different base
        haplotypes.append("".join(haplotype))
    return np.array(haplotypes)


def create_genome(length: int) -> str:
    random.seed(43)
    return ''.join(random.choices("ACTG", k=length))


def make_graph(haplotypes: np.ndarray, weights: np.ndarray) -> Graph:
    """
    Returns graph in graph-tool format.
    """
    g = Graph(directed=True) #Define a graph with its properties
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
                    g.ep.ori[e] = '+'  # Orientation of edge
            if i % 8 == 0:
                contig = haplotypes[j][i::]
                if len(contig) > 13:
                    contig = contig[:12]
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

def get_contig_paths_from_graph(g: Graph, contig_list: list) -> list:
    contigs = []
    for contig in contig_list:
        last_vertex = None
        contig_path = []
        for v in g.vertices():
            if contig in g.vp.contigs[v]:
                if last_vertex:
                    contig_path.append((int(last_vertex), int(v)))
                last_vertex = v
        contigs.append(contig_path)
    return contigs


def contract_vertices(g: Graph) -> Graph:
    edges_to_contract = list(g.edges())
    vertices_to_remove = []
    for s,t in edges_to_contract:
        if (s.out_degree() == 1 and t.in_degree() == 1) and g.vp.weight[s] == g.vp.weight[t]:
            g.vp.seq[t] = g.vp.seq[s] + g.vp.seq[t]
            g.vp.contigs[t] = list(set(g.vp.contigs[s]) | set(g.vp.contigs[t]))
            for e in s.in_edges():
                new_start_vertex = e.source()
                g.add_edge(new_start_vertex, t)
            vertices_to_remove.append(s)
    g.remove_vertex(reversed(sorted(vertices_to_remove)), fast=False)
    g = g.copy()  # Create a copy to reindex vertices

    # Reindex vertices
    # g.vp.new_index = g.new_vertex_property("int")  # Temporary property to store new indices
    # for new_idx, v in enumerate(g.vertices()):
    #     g.vp.new_index[v] = new_idx

    # g.vertex_index = g.vp.new_index  # Overwrite the default index mapping

    return g


if __name__ == '__main__':
    sys.exit(main())