import networkx as nx
import argparse
import flowpaths as fp

def main():
    parser = argparse.ArgumentParser(description="Convert a GFA file to a networkx graph")
    parser.add_argument('-n', '--graph_number', dest='graph_number', type=int, required=True, help="graph number")
    # parser.add_argument('-g', '--gfa_file', dest='gfa_file', type=str, required=True, help="GFA file")
    # parser.add_argument('-a', '--abundance_file', dest='abundance_file', type=str, required=True, help="Abundance file")
    args = parser.parse_args()
    # network, subpaths = convert_gfa_to_networkx(args.gfa_file, args.abundance_file)
    # network = add_source_and_sink_node(network)
    graphs_files = [
        # ("output/vg-flow/data/graph_100_3.gfa", "output/vg-flow/data/abundance_100_3.txt"),
        ("simulated_data/gridsearch/graph_1000_6.gfa", "simulated_data/gridsearch/abundances_1000_6.txt"),
        ("simulated_data/gridsearch/graph_1000_7.gfa", "simulated_data/gridsearch/abundances_1000_7.txt"),
        ("simulated_data/gridsearch/graph_1000_8.gfa", "simulated_data/gridsearch/abundances_1000_8.txt"),
        ("simulated_data/gridsearch/graph_10000_6.gfa", "simulated_data/gridsearch/abundances_10000_6.txt"),
        ("simulated_data/gridsearch/graph_10000_7.gfa", "simulated_data/gridsearch/abundances_10000_7.txt"),
        ("simulated_data/gridsearch/graph_10000_8.gfa", "simulated_data/gridsearch/abundances_10000_8.txt"),
        ("simulated_data/gridsearch/graph_100000_6.gfa", "simulated_data/gridsearch/abundances_100000_6.txt"),
        ("simulated_data/gridsearch/graph_100000_7.gfa", "simulated_data/gridsearch/abundances_100000_7.txt"),
        ("simulated_data/gridsearch/graph_100000_8.gfa", "simulated_data/gridsearch/abundances_100000_8.txt"),
    ]

    # for graph_file in graphs_files:
    #     process_graph(graph_file[0], graph_file[1])
    process_graph(graphs_files[args.graph_number][0], graphs_files[args.graph_number][1])


    
def process_graph(graph_file: str, abundance_file: str):

    graph, subpath_constraints_nodes = convert_gfa_to_networkx(graph_file, abundance_file)

    # We create a node expanded graph, where the weights are taken from the attribute "flow"
    neGraph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")

    # We transform the constraints into constraints in the node expanded graph
    ne_subpath_constraints_nodes = neGraph.get_expanded_subpath_constraints_nodes(subpath_constraints_nodes)

    # Fastest optimizations for these inputs
    optimization_options = {
        "optimize_with_safe_paths": False,
        "optimize_with_safe_sequences": False,
        "optimize_with_safe_zero_edges": False,
        "optimize_with_greedy": False,
        "external_solver": "gurobi",
        "threads": 8,
    }

    # # Greedy optimization allowed to kick in
    # optimization_options = {
    #     "optimize_with_safe_paths": False,
    #     "optimize_with_safe_sequences": True,
    #     "optimize_with_safe_zero_edges": True,
    #     "optimize_with_greedy": True,
    #     "external_solver": "highs",
    #     "threads": 4,
    # }

    # We create the solver and solve the model
    mfd_model = fp.MinFlowDecomp(
        neGraph, 
        # weight_type=int,
        # weight_type=float, # default is float
        flow_attr="flow",
        edges_to_ignore=neGraph.edges_to_ignore,
        subpath_constraints_coverage=1,
        subpath_constraints=ne_subpath_constraints_nodes,
        optimization_options=optimization_options,
        )
    mfd_model.solve()
    process_solution(graph, neGraph, mfd_model)

def process_solution(graph: nx.DiGraph, neGraph: fp.NodeExpandedDiGraph, model: fp.MinFlowDecomp):

    if model.is_solved():
        solution = model.get_solution()
        expanded_paths = solution["paths"]
        paths = neGraph.get_condensed_paths(expanded_paths)

        # Printing the paths and weights
        print("Paths:", paths)
        print("Weights:", solution["weights"])
        
        # Double checking that the model returned a valid flow decomposition
        print("model.is_valid_solution()", model.is_valid_solution())
        
        # printing some running times
        print(model.solve_statistics)

        # # Basic drawing of the solution paths
        # fp.utils.graphutils.draw_solution_basic(graph, flow_attr="flow", paths = paths, weights = solution["weights"], id = model.G.graph["id"])
    else:
        print("Model could not be solved.")

def convert_gfa_to_networkx(gfa_file, abundance_file):
    network = nx.DiGraph()
    network.graph["id"] = gfa_file
    subpaths = []
    with open(abundance_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(':')
            network.add_node(line[0], flow=float(line[1]))
        
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'L':
                network.add_edge(line[1], line[3])
            if line[0] == 'P':
                path_nodes = [node.strip('+') for node in line[2].split(',')]
                subpaths.append(path_nodes)
    return network, subpaths     


if __name__ == "__main__":
    main()