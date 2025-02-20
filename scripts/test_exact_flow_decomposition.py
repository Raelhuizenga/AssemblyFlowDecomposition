import unittest
from graph_tool.all import Graph
from exact_flow_decomposition import fd_fixed_size, add_source_and_sink, get_sequence_from_path, induce_flow_on_graph, draw_graph
import numpy as np
from collections import Counter
from simulate_data import make_graph

class TestFlowDecomposition(unittest.TestCase):


    def test_exact_mfd_two_paths(self):
        haplotypes = np.array(['ATG', 'ATG', 'AGG'])
        true_haplopyes = np.array(['ATG', 'AGG'])
        true_abundance = np.array([2, 1])
        g = make_graph(haplotypes)
        g = add_source_and_sink(g)
        subpaths = [[(0,1)], [(0,2)], [(1,3)], [(2,3)]]
        abundances, paths, _ = fd_fixed_size(g, subpaths,  2)
        paths = [get_sequence_from_path(g, path) for path in paths]
        self.assertTrue(np.array_equal(np.sort(true_abundance), np.sort(np.array(abundances))))
        self.assertTrue(np.array_equal(np.sort(true_haplopyes), np.sort(np.array(paths))))

    # def test_inexact_mfd_two_paths(self):
    #     haplotypes = np.array(['ATG', 'ATG', 'AGG'])
    #     true_haplopyes = np.array(['ATG', 'AGG'])
    #     true_abundance = np.array([2, 1])
    #     g = make_graph(haplotypes)
    #     first_vertex = g.vertex(0)
    #     g.vp.weight[first_vertex] += 1
    #     g = add_source_and_sink(g)
    #     draw_graph(g)
    #     _, _, g = induce_flow_on_graph(g)
      
    #     abundances, paths, _ = fd_fixed_size(g, [],  2)
    #     paths = [get_sequence_from_path(g, path) for path in paths]
    #     self.assertTrue(np.array_equal(np.sort(true_abundance), np.sort(np.array(abundances))))
    #     self.assertTrue(np.array_equal(np.sort(true_haplopyes), np.sort(np.array(paths))))

if __name__ == '__main__':
    unittest.main()