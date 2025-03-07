import sys
import os
from graph_tool.all import Graph
from graph_functions import read_gfa, rev_comp, write_gfa
from collections import OrderedDict
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from graph_tool.draw import graph_draw, sfdp_layout
from graph_tool.topology import is_DAG
import time
import json
import argparse

usage = """
    This script is used to decompose a flow graph into a set of paths.
"""

def main():
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-t', '--threads', dest='threads', type=int, required=True, help="Number of threads to use, will use all available if 0")
    parser.add_argument('-g', '--graph', dest='graph', type=str, required=True, help="Input graph file")
    parser.add_argument('-a', '--abundance', dest='abundance', type=str, required=True, help="Input abundance file")
    parser.add_argument('-o', '--output_folder', dest='output_folder', type=str, required=True, help="Output folder")
    parser.add_argument('-i', '--input_folder', dest='input_folder', type=str, required=True, help="Input folder")
    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()

    test_with_read_from_file(args.input_folder, args.graph, args.abundance, args.output_folder, threads)


def test_with_read_from_file(data_folder: str, gfa_file: str, abundance_file: str, output_folder: str, threads):
    # read abundances from file
    abundance_list = read_abundances(f'{data_folder}/{abundance_file}')

    # read graph from file
    graph, paths = read_gfa(f'{data_folder}/{gfa_file}')[:2]

    # reformate paths
    paths_edges = create_contig_paths(paths)
 
    graph.vp.weight = graph.new_vertex_property('int', 0)
    for i, v in enumerate(list(graph.vertices())[:-2]):
        graph.vp.weight[v] = abundance_list[i]

    # count contigs and vertices
    nvert = graph.num_vertices()
    contig_IDs = set()
    for contigs in graph.vp.contigs:
        contig_IDs = contig_IDs.union(set(contigs))
    ncontigs = len(contig_IDs)
    average_degree = float(sum(graph.get_total_degrees(graph.get_vertices())) / graph.num_vertices())
    max_degree = int(max(graph.get_total_degrees(graph.get_vertices())))

    print("********* input graph *********")
    print("#vertices = {}".format(nvert))
    print("#edges = {}".format(graph.num_edges()))
    print("#contigs = {}".format(ncontigs))
    print("average degree graph = {}".format(average_degree))
    print("max degree graph = {}".format(max_degree))
    print("*******************************\n")

    graph = add_source_and_sink(graph)

    x, obj, g = induce_flow_on_graph(graph, threads=threads)

    print("Objective value: ", obj)
    
    assert(is_DAG(g))

    start_time = time.time()
    w_sol, sol_paths, z_sol = mfd_algorithm(graph, paths_edges, max_strains=30, threads=threads)
    sol_paths = [get_sequence_from_path(graph, path) for path in sol_paths]

    end_time = time.time()

    print("Execution Time: {:.2f} seconds".format(end_time - start_time))

    # save solution paths and weights to file
    with open(f"{output_folder}/solution_{gfa_file[:-3]}.json", "w") as f:
        json.dump({"weights": w_sol, "paths": sol_paths, "time": end_time-start_time, "num_vertices": graph.num_vertices(), "num_edges": graph.num_edges(), "average_degree": average_degree, "max_degree": max_degree}, f)

    with open(f"{output_folder}/solution_paths_{gfa_file[:-3]}.final.fasta", "w") as f:
        for i in range(len(w_sol)):
            f.write(f'>path{i} {w_sol[i]}x freq={round(w_sol[i] / sum(w_sol), 3)}\n{sol_paths[i]}\n')


def induce_flow_on_graph(graph: Graph, threads=1) -> Graph:
    # define new model
    m = gp.Model('lp')
    # define objective
    obj = gp.LinExpr()
    # edge variables
    x = m.addVars(list(graph.edges()), lb=0, vtype=GRB.CONTINUOUS, name='x')
    # add additional variables to implement absolute values
    y = m.addVars(list(graph.vertices()), lb=0, vtype=GRB.CONTINUOUS, name='y')

    # set constraints for flow conservation
    for v in graph.vertices():
        if graph.vp.seq[v] == 's' or graph.vp.seq[v] == 't':
            # source/sink node
            continue
        sum_flow = gp.LinExpr()
        for u in v.in_neighbors():
            sum_flow += x[graph.edge(u,v)]
        obj += y[v] #abs(abundance - sum_xv)
        # set constraints on y[v] to obtain absolute value
        m.addConstr(y[v] >= sum_flow - graph.vp.weight[v], "y_{}_-".format(v))
        m.addConstr(y[v] >=  graph.vp.weight[v] - sum_flow, "y_{}_+".format(v))
        
        for w in v.out_neighbors():
            sum_flow -= x[graph.edge(v,w)]
        m.addConstr(sum_flow == 0)

    # set objective and minimize
    m.setObjective(obj, GRB.MINIMIZE)
    print('\nObjective function ready, starting Gurobi optimization:\n')
    m.update()

    m.Params.LogToConsole = 0
    m.Params.Threads = threads
    m.Params.NumericFocus = 0
    m.Params.PoolSearchMode = 0
    m.Params.PoolSolutions = 10
    m.Params.Method = 4

    # solve ILP
    print("\n*** Running optimization to induce flow***\n")
    m.optimize()

    print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL:
        x_final = {edge: v.x for edge, v in zip(graph.edges(), m.getVars()) if 'x' in v.varName}
        objVal = m.objVal
        for v in graph.vertices():
            new_v_value = 0
            if graph.vp.seq[v] == 't':
                continue 
            for u in v.in_neighbors():
                new_v_value += x_final[graph.edge(u,v)]
            graph.vp.weight[v] = new_v_value
        return x_final, objVal, graph

    else:
        try:
            m.computeIIS()
            # Print the names of all of the constraints in the IIS set.
            print("IIS constraints:")
            for c in m.GetConstrs():
                if c.Get(GRB.IntAttr.IISConstr) > 0:
                    print(c.Get(GRB.StringAttr.ConstrName))
            # Print the names of all of the variables in the IIS set.
            print("IIS variables:")
            for v in m.GetVars():
                if v.Get(GRB.IntAttr.IISLB) > 0 or v.Get(GRB.IntAttr.IISUB) > 0:
                    print(v.Get(GRB.StringAttr.VarName))
            print("ERROR: Infeasible model.")
        except:
            #print(m.getAttr(GRB.Attr.UnbdRay, m.getVars()))
            print("ERROR: Unbounded model.")

        print('\nNo optimal solution found, exiting.')
        sys.exit(1)

def add_source_and_sink(graph) -> Graph:
    # add source and sink node
    source_node = graph.add_vertex()
    graph.vp.seq[source_node] = 's'
    graph.vp.weight[source_node] = 0
    sink_node = graph.add_vertex()
    graph.vp.seq[sink_node] = 't'
    graph.vp.weight[sink_node] = 0

     # add edges to source and sink node
    for i, v in enumerate(list(graph.vertices())[:-2]):
        if v.in_degree() == 0:
            e = graph.add_edge(source_node, v)
            graph.ep.ori[e] = '+'
        if v.out_degree() == 0 and v != source_node and v != sink_node:
            graph.add_edge(v, sink_node)
            graph.ep.ori[e] = '+'
    return graph


def create_contig_paths(paths) -> list:
    # @ToDo check if this function is right, does it take the right index if i > 10?
    formatted_paths = []

    for _, path in paths.items():
        
        # If pairwise tuples are needed, create them
        if len(path) >= 2:
            formatted_paths.append([(path[i][0], path[i+1][0]) if path[i][1] == '+' else (path[i+1][0], path[i][0]) for i in range(len(path)-1)])
        else:
            formatted_paths.append(path[0])

    return formatted_paths 


def read_abundances(abundance_file: str) -> list:
    # read node abundances from file
    abundance_list = []
    old_id = -1
    
    with open(abundance_file, 'r') as f:
        for line in f:
            [node_id, abundance] = line.rstrip().split(':')
            # make sure that nodes are ordered incrementally
            new_id = int(node_id)
            if old_id != -1:
                assert new_id == old_id + 1
            old_id = new_id
            # add node abundance to list
            abundance_list.append(float(abundance))
    return abundance_list


def get_sequence_from_path(graph, path):
    seq = ''
    for s,t in path[:-1]:
         ori = graph.ep.ori[graph.edge(s,t)]
         seq += get_seq(graph, t, ori)
    return seq


def get_seq(g, v, ori):
    """Extract nucleotide sequence from graph for given node and orientation"""
    seq = g.vp.seq[g.vertex(int(v))]
    if ori == "-":
        seq = rev_comp(seq)
    return seq


def get_source_and_sink(graph):
    sources = [v for v in graph.vertices() if v.in_degree() == 0]
    sinks = [v for v in graph.vertices() if v.out_degree() == 0]
    return sources[0], sinks[0]


def draw_graph(g):
    pos = sfdp_layout(g)

    # Create a new property map combining `seq` and `weight`
    text_prop = g.new_vertex_property(
        "string", 
         [f"{g.vp.seq[v]}, {g.vp.weight[v]}" for v in g.vertices()]
    )

    fontsize_prop = g.new_vertex_property("double")
    fontsize_prop.a = 20

    graph_draw(
        g,
        pos=pos,
        vertex_text=text_prop,
        vertex_font_size=fontsize_prop,
        output_size=(1000, 1000),
        bg_color=[1, 1, 1, 1],  # white background (RGBA format)
        output="flow_graph.png"
    )


def mfd_algorithm(graph: Graph, subpaths, max_strains=10, threads=0) -> tuple:
    for i in range(2, max_strains): 
        w_sol, paths, z_sol = fd_fixed_size(graph, subpaths, i, threads=threads)
        if w_sol != []:
            return w_sol, paths, z_sol
    return w_sol, paths, z_sol


def build_base_ilp_model(graph, subpaths, size, use_only_subpaths=False):
    max_flow_value = max(graph.vp.weight)

    # create extra sets
    T = [(u, v, k) for (u, v) in graph.get_edges() for k in range(size)]
    SC = list(range(size))
    R = [(k,s) for k in range(size) for s in range(len(subpaths))]

    # Create a new model
    model = gp.Model('MFD')
    model.setParam('LogToConsole', 0)


    # Create variables
    x = model.addVars(T, vtype=GRB.BINARY, name='x')                    # Binary variable, 1 if an edge (u,v) is in a path, 0 otherwise
    w = model.addVars(SC, vtype=GRB.INTEGER, name='w', lb=1)            # Weight of path i
    z = model.addVars(T, vtype=GRB.CONTINUOUS, name='z', lb=0)          # Linearized w_i * xᵤᵥᵢ
    r = model.addVars(R, vtype=GRB.BINARY,name='r')                     # Binary variable, 1 if subpath s in path k, 0 otherwise


    # flow conservation
    for k in range(size):
        # model.addConstr(sum(x[u, v, k] for (u, v) in graph.get_edges()) >= 1) # ensure at least one edge is in the path
        for v in graph.get_vertices():
            if graph.vp.seq[v] == 's':
                model.addConstr(sum(x[v, w2, k] for _, w2 in graph.get_out_edges(v)) == 1, name="One edge out of source")
            elif graph.vp.seq[v] == 't':
                model.addConstr(sum(x[u, v, k] for u, _ in graph.get_in_edges(v)) == 1, name="One edge into sink")
            else:
                model.addConstr(sum(x[v, w2, k] for _, w2 in graph.get_out_edges(v)) - sum(x[u, v, k] for u, _ in graph.get_in_edges(v)) == 0, name=f"Flow conservation {v} for path {k}")

    # flow balance
    for u in graph.get_vertices():
        weight = graph.vp.weight[u]
        if weight > 0:
            model.addConstr(weight == sum(z[u, v, k] for k in range(size) for v in graph.get_out_neighbors(u)), name=f"Flow balance {v}")

    # supbatph constraints
    for k in range(0,size):
        if use_only_subpaths:
            for s, subpath in enumerate(subpaths):
                model.addConstr(sum(x[u,v,k] for u,v in subpath) == len(subpath)*r[k,s], name=f"Strict subpath {s} in path {k}")
            for (u, v) in graph.get_edges():  # Iterate over all possible edges
                if graph.vp.seq[u] != 's' and graph.vp.seq[v] != 't':
                    model.addConstr(
                        x[u,v,k] <= sum(r[k,s] for s, subpath in enumerate(subpaths) if (u,v) in subpath),
                        name=f"Edge ({u},{v}) in path {k} covered by subpath"
                    )
        else:
            for s, subpath in enumerate(subpaths):
                model.addConstr(sum(x[u,v,k] for u,v in subpath) >= len(subpath)*r[k,s], name=f"Length subpath {s} in path {k}")
        
    if not use_only_subpaths:	
        model.addConstrs((sum(r[k,s] for k in range(size)) >= 1 for s in range(len(subpaths))), name="At least one path contains subpath")

    # linearization
    for (u, v) in graph.get_edges():
        for k in range(size):
            # linearization - x*w
            model.addConstr(z[u, v, k] <= max_flow_value * x[u, v, k], name=f"Linearization 1 {u} {v} {k}")
            model.addConstr(w[k] - (1 - x[u, v, k]) * max_flow_value <= z[u, v, k], name=f"Linearization 2 {u} {v} {k}")
            model.addConstr(z[u, v, k] <= w[k], name=f"Linearization 3 {u} {v} {k}, flow balance")
    # # Add symmetry breaking constraints        
    model.addConstrs(w[k] <= w[k+1] for k in range(size-1))
    model.setObjective(gp.quicksum(w[k] for k in SC), GRB.MINIMIZE)
    return model, x, w, z


def fd_fixed_size(graph, subpaths, size, threads = 0) -> tuple:
    w_sol = []
    paths = []
    z_sol = {}
    source, sink = get_source_and_sink(graph)

    # calculate a flow decomposition into size paths
    try:
        # Create a new model
        model, _, _, _ = build_base_ilp_model(graph, subpaths, size)
        model.setParam('Threads', threads) 

        # objective function
        model.optimize()

        if model.status == GRB.OPTIMAL:
            T = [(u, v, k) for (u, v) in graph.get_edges() for k in range(size)]
            w_sol = [0] * size
            paths = [list() for _ in range(size)]
            for k in range(size):
                w_sol[k] = round(model.getVarByName(f'w[{k}]').x)
            for (u, v, k) in T:
                z_var = model.getVarByName(f'z[{u},{v},{k}]')
                if z_var is not None:  # Store the flow values
                    z_sol[(u, v, k)] = round(z_var.x)  # Adjust precision as needed
                if round(model.getVarByName(f'x[{u},{v},{k}]').x) == 1:
                    paths[k].append((int(u), int(v)))
            for k in range(len(paths)):
                paths[k] = reorder_path(paths[k], int(source), int(sink))
        else:
            print(f'Model is infeasible for size {size}')
            model.computeIIS()
            # for c in model.getConstrs():
            #     if c.IISConstr:  # Identifies constraints in the IIS
            #         print(f"Infeasible constraint {c.constrName}: {model.getRow(c)} <= {c.RHS}")
            # for v in model.getVars():
            #     if v.IISLB:
            #         print(f"Lower bound causing infeasibility: {v.varName} >= {v.lb}")
            #     if v.IISUB:
            #         print(f"Upper bound causing infeasibility: {v.varName} <= {v.ub}")


    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return w_sol, paths, z_sol


def reorder_path(edges, start_vertex, end_vertex):
    path = []
    current_vertex = start_vertex
    edge_dict = {u: v for u, v in edges}  # Map source to target for quick traversal
    while current_vertex in edge_dict:
        next_vertex = edge_dict[current_vertex]
        path.append((current_vertex, next_vertex))
        if next_vertex == end_vertex:
            break
        current_vertex = next_vertex
    return path


if __name__ == '__main__':
    sys.exit(main())
