import unittest
from graph_tool.all import Graph
from scripts.flow_decomposition_old import draw_graph, build_contig_flow_graph, exact_flow_decomposition_fixed_size, convert_flow_path_to_contig_path
import numpy as np
from collections import Counter

class TestFlowDecomposition(unittest.TestCase):


    def test_exact_mfd_two_paths(self):
        adj_out = {'0': ['4'], '1': ['5'], '2': ['1'], '4': [], '5': []}
        contig_abundances = {'0': 25.0, '1': 29.0, '2': 29.0, '4': 25.0, '5': 29.0}
        true_abundances = np.array([25, 29])
        true_paths = [[(0, 6), (6, 7), (7, 4), (4, 5), (5, 10), (10, 11), (11, 1)], 
                      [(0, 2), (2, 3), (3, 8), (8, 9), (9, 1)]]
        graph, _ = build_contig_flow_graph(adj_out, contig_abundances)
        abundances, paths, slack = exact_flow_decomposition_fixed_size(graph, 2)
        self.assertTrue(np.array_equal(np.sort(true_abundances), np.sort(np.array(abundances))))
        self.assertEqual(Counter(map(tuple, true_paths)), Counter(map(tuple, paths)))
        self.assertEqual(np.sum(slack), 0)

    def test_inexact_mfd_two_paths(self):
        adj_out = {'0': ['4'], '1': ['5'], '2': ['1'], '4': [], '5': []}
        contig_abundances = {'0': 25.0, '1': 29.0, '2': 30.0, '4': 25.0, '5': 29.0}
        true_abundances = np.array([25, 29])
        true_paths = [[(0, 2), (2, 3), (3, 8), (8, 9), (9, 1)], 
                      [(0, 6), (6, 7), (7, 4), (4, 5), (5, 10), (10, 11), (11, 1)]]
        graph, _ = build_contig_flow_graph(adj_out, contig_abundances)
        abundances, paths, slack = exact_flow_decomposition_fixed_size(graph, 2)
        self.assertTrue(np.array_equal(np.sort(true_abundances), np.sort(np.array(abundances))))
        self.assertEqual(Counter(map(tuple, true_paths)), Counter(map(tuple, paths)))
        self.assertEqual(np.sum(slack), 1)

    def test_exact_mfd_three_paths(self):
        adj_out = {'0': ['4'], '1': ['5', '4'], '2': ['1'], '4': [], '5': []}
        contig_abundances = {'0': 25.0, '1': 31.0, '2': 31.0, '4': 27.0, '5': 29.0}
        true_abundances = np.array([25, 29, 2])
        true_paths = [[(0, 6), (6, 7), (7, 4), (4, 5), (5, 10), (10, 11), (11, 1)], 
                      [(0, 2), (2, 3), (3, 8), (8, 9), (9, 1)], 
                      [(0, 6), (6, 7), (7, 4), (4, 5), (5, 8), (8, 9), (9, 1)]]
        graph, _ = build_contig_flow_graph(adj_out, contig_abundances)
        abundances, paths, slack = exact_flow_decomposition_fixed_size(graph, 3)
        self.assertTrue(np.array_equal(np.sort(true_abundances), np.sort(np.array(abundances))))
        self.assertEqual(Counter(map(tuple, true_paths)), Counter(map(tuple, paths)))
        self.assertEqual(np.sum(slack), 0)

    def test_inexact_mfd_three_paths(self):
        adj_out = {'0': ['4'], '1': ['5', '4'], '2': ['1'], '4': [], '5': []}
        contig_abundances = {'0': 26.0, '1': 31.0, '2': 30.0, '4': 27.0, '5': 30.0}
        true_abundances = np.array([29, 26, 1])
        true_paths = [[(0, 2), (2, 3), (3, 8), (8, 9), (9, 1)], 
                      [(0, 6), (6, 7), (7, 4), (4, 5), (5, 10), (10, 11), (11, 1)], 
                      [(0, 6), (6, 7), (7, 4), (4, 5), (5, 8), (8, 9), (9, 1)]]
        graph, _ = build_contig_flow_graph(adj_out, contig_abundances)
        abundances, paths, slack = exact_flow_decomposition_fixed_size(graph, 3)
        self.assertTrue(np.array_equal(np.sort(true_abundances), np.sort(np.array(abundances))))
        self.assertEqual(Counter(map(tuple, true_paths)), Counter(map(tuple, paths)))
        self.assertEqual(np.sum(slack), 1)

    def test_convert_flow_path_to_contig_path(self):
        adj_out = {'0': ['4'], '1': ['5'], '2': ['1'], '4': [], '5': []}
        contig_abundances = {'0': 25.0, '1': 29.0, '2': 29.0, '4': 25.0, '5': 29.0}
        true_contig_paths = [['2', '1', '5'], ['0', '4']]
        graph, node_dict = build_contig_flow_graph(adj_out, contig_abundances)
        _, paths, _ = exact_flow_decomposition_fixed_size(graph, 2)
        contig_paths = convert_flow_path_to_contig_path(graph, node_dict, paths)
        self.assertEqual(contig_paths, true_contig_paths)

if __name__ == '__main__':
    unittest.main()