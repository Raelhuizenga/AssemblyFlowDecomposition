import networkx as nx
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert a GFA file to a networkx graph")
    parser.add_argument('-g', '--gfa_file', dest='gfa_file', type=str, required=True, help="GFA file")
    parser.add_argument('-a', '--abundance_file', dest='abundance_file', type=str, required=True, help="Abundance file")
    args = parser.parse_args()
    network, subpaths = convert_gfa_to_networkx(args.gfa_file, args.abundance_file)
    network = add_source_and_sink_node(network)
    

def convert_gfa_to_networkx(gfa_file, abundance_file):
    network = nx.DiGraph()
    subpaths = []
    with open(abundance_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(':')
            network.add_node(line[0], flow=line[1])
        
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'L':
                network.add_edge(line[1], line[3])
            if line[0] == 'P':
                path_nodes = [node.strip('+') for node in line[2].split(',')]
                subpaths.append(path_nodes)
    return network, subpaths

def add_source_and_sink_node(graph):
    graph.add_node('s')
    graph.add_node('t')
    for node in graph.nodes():
        if graph.in_degree(node) == 0 and node != 's':
            graph.add_edge('s', node)
        if graph.out_degree(node) == 0 and node != 't':
            graph.add_edge(node, 't')
    return graph         


if __name__ == "__main__":
    main()