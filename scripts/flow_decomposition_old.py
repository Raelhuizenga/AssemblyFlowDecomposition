import sys
import numpy as np
from graph_tool.all import Graph
from graph_functions import read_gfa
from graph_tool.draw import graph_draw, sfdp_layout
from collections import OrderedDict
import gurobipy as gp
from gurobipy import GRB



def main():
    # Example call to build_flow_graph
    adj_out = {'0': ['4'], '1': ['5', '4'], '2': ['1'], '4': [], '5': []}
    contig_abundances = {'0': 26.0, '1': 31.0, '2': 30.0, '4': 27.0, '5': 30.0}
    # Build the graph
    g, node_dict = build_contig_flow_graph(adj_out, contig_abundances)

    abundances, paths, slack = mfd_algorithm(g)

    print(abundances, paths, slack)

    draw_graph(g)

def draw_graph(g):
    pos = sfdp_layout(g)
    graph_draw(
        g,
        pos=pos,
        vertex_text=g.vp.id,
        edge_text=g.ep.weight,
        output_size=(1000, 1000),
        vertex_fill_color=g.vp.color,
        bg_color=[1, 1, 1, 1],  # white background (RGBA format)
        output="flow_graph.png"
    )

def exact_flow_decomposition_fixed_size(graph: Graph, size: int):
    w_sol = []
    paths = []
    slack = []

    # calculate a flow decomposition into size paths
    try:
        # Create a new model
        model, _, _, _ = build_base_ilp_model(graph, size)

        # objective function
        model.optimize()

        if model.status == GRB.OPTIMAL:
            T = [(u, v, k) for (u, v) in graph.get_edges() for k in range(size)]
            w_sol = [0] * len(range(size))
            slack = [0] * len(range(size))
            paths = [list() for _ in range(size)]
            for k in range(size):
                w_sol[k] = round(model.getVarByName(f'w[{k}]').x)
                slack[k] = round(model.getVarByName(f'pho[{k}]').x)
            for (u, v, k) in T:
                if round(model.getVarByName(f'x[{u},{v},{k}]').x) == 1:
                    paths[k].append((int(u), int(v)))
            for k in range(len(paths)):
                paths[k] = reorder_path(paths[k])

    except gp.GurobiError as e:
        print(f'Error code {e.errno}: {e}', file=sys.stderr)

    except AttributeError:
        print('Encountered an attribute error', file=sys.stderr)

    return w_sol, paths, slack


def mfd_algorithm(graph: Graph, max_strains=10):
    tolerance = 1e-2
    old_w_sol, old_paths, old_slack = [], [], -1
    for i in range(2, max_strains): # Note: can be improved by looping over edge width instead of number of edges @TODO
        w_sol, paths, slack = exact_flow_decomposition_fixed_size(graph, i)
        print(w_sol)
        if w_sol != []:
            if abs(np.sum(slack) - old_slack) < tolerance:
                return old_w_sol, old_paths, old_slack # return when the slack does not change
            old_w_sol, old_paths, old_slack = w_sol, paths, np.sum(slack)
    return w_sol, paths

def build_base_ilp_model(graph : Graph, size : int):
    if not graph.ep.weight.a.size:
        raise ValueError("Edge weights are not properly defined.")
    max_flow_value = max(graph.ep.weight.a)
    M = 1e3

    # create extra sets
    T = [(u, v, k) for (u, v) in graph.get_edges() for k in range(size)]
    SC = list(range(size))

    # Create a new model
    model = gp.Model('MFD')
    model.setParam('LogToConsole', 0)
    model.setParam('Threads', 1) # use one for now (should be fixed with input parameter later)


    # Create variables
    x = model.addVars(T, vtype=GRB.BINARY, name='x')                    # Binary variable, 1 if an edge (u,v) is in a path, 0 otherwise
    w = model.addVars(SC, vtype=GRB.INTEGER, name='w', lb=0)            # Weight of path i
    pho = model.addVars(SC,vtype=GRB.INTEGER,name="pho",lb=0)           # Slack of path i
    z = model.addVars(T, vtype=GRB.CONTINUOUS, name='z', lb=0)          # Linearized w_i * xᵤᵥᵢ
    phi = model.addVars(T, vtype=GRB.CONTINUOUS, name='phi', lb=0)      # Linearized ρ_i * xᵤᵥᵢ


    # flow conservation
    for k in range(size):
        for v in graph.get_vertices():
            if graph.vp.id[v] == 's':
                model.addConstr(sum(x[v, w, k] for _, w in graph.get_out_edges(v)) == 1)
            elif graph.vp.id[v] == 't':
                model.addConstr(sum(x[u, v, k] for u, _ in graph.get_in_edges(v)) == 1)
            else:
                model.addConstr(sum(x[v, w, k] for _, w in graph.get_out_edges(v)) - sum(x[u, v, k] for u, _ in graph.get_in_edges(v)) == 0)

    # flow balance
    for e in graph.get_edges():
        weight = graph.ep.weight[e]
        if weight > 0:
            u, v = e
            model.addConstr(weight - sum(z[u, v, k] for k in range(size)) <= sum(phi[u, v, k] for k in range(size)))
            model.addConstr(weight - sum(z[u, v, k] for k in range(size)) >= - sum(phi[u, v, k] for k in range(size)))


    for (u, v) in graph.get_edges():
        for k in range(size):
            # linearization - x*w
            model.addConstr(z[u, v, k] <= max_flow_value * x[u, v, k])
            model.addConstr(w[k] - (1 - x[u, v, k]) * max_flow_value <= z[u, v, k])
            model.addConstr(z[u, v, k] <= w[k])

            # linearization - x*pho
            model.addConstr(phi[u, v, k] <= M * x[u, v, k])
            model.addConstr(pho[k] - (1 - x[u, v, k]) * M <= phi[u, v, k])
            model.addConstr(phi[u, v, k] <= pho[k])
    
    # Set the objective to minimize the sum of slacks
    model.setObjective(gp.quicksum(pho[k] for k in SC), GRB.MINIMIZE)

    return model, x, w, z


def build_contig_flow_graph(adj_out, contig_abundances):
    # build flow graph
    g = Graph(directed=True)
    g.vp.id = g.new_vertex_property('string')
    g.ep.weight = g.new_edge_property('int')
    g.vp.color = g.new_vertex_property("vector<double>")

    # Add global source and sink
    source_node = g.add_vertex()
    g.vp.id[source_node] = 's'
    g.vp.color[source_node] = [1, 0, 0, 1] 
    sink_node = g.add_vertex()
    g.vp.id[sink_node] = 't'
    g.vp.color[sink_node] = [1, 0, 0, 1] 

    # Add 2 nodes for every contig
    node_dict = {}
    for v in adj_out.keys():
        v_in = g.add_vertex()
        g.vp.color[v_in] = [0, 0, 1, 1] 
        g.vp.id[v_in] = str(int(v_in))  # Use implicit vertex index
        v_out = g.add_vertex()
        g.vp.color[v_out] = [0, 0, 1, 1] 
        g.vp.id[v_out] = str(int(v_out))  # Use implicit vertex index
        e = g.add_edge(v_in, v_out)  # contig-edge
        g.ep.weight[e] = contig_abundances[v]  # contig abundance
        node_dict[v] = (v_in, v_out)

    # Add edges for overlaps
    for v, adj_list in adj_out.items():
        v_in, v_out = node_dict[v]
        if len(adj_list) == 0:
            e = g.add_edge(v_out, sink_node)  # sink-edge
            g.ep.weight[e] = -1
        for w in adj_list:
            w_in, w_out = node_dict[w]
            e = g.add_edge(v_out, w_in)  # overlap-edge
            g.ep.weight[e] = -1

    # Add edges for source nodes
    for v in adj_out.keys():
        v_in, v_out = node_dict[v]
        if v_in.in_degree() == 0:
            e = g.add_edge(source_node, v_in)  # source-edge
            g.ep.weight[e] = -1

    return g, node_dict

def reorder_path(edges):
    path = []
    current_vertex = 0
    edge_dict = {u: v for u, v in edges}  # Map source to target for quick traversal
    while current_vertex in edge_dict:
        next_vertex = edge_dict[current_vertex]
        path.append((current_vertex, next_vertex))
        if next_vertex == 1:
            break
        current_vertex = next_vertex
    return path


def convert_flow_path_to_contig_path(graph: Graph, nodeDict: dict, paths: list):
    contig_paths = []
    for path in paths:
        contig_path = []
        for i in range(len(path)):
            v, w = path[i][0], path[i][1]
            if graph.vp.id[v] == 's' or graph.vp.id[w] == 't':
                continue
            for contig, (v_in, v_out) in nodeDict.items():
                if v_in == v:
                    contig_path.append(contig)
                    break
        contig_paths.append(contig_path)    
    return contig_paths


if __name__ == '__main__':
    sys.exit(main())