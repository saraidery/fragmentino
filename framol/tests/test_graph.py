import numpy as np
import pytest
import os


from framol import SimpleWeightedGraph


class TestGraph:
    def test_add_edge_1(self):

        g = SimpleWeightedGraph()
        g.add_vertex("1")
        g.add_vertex("2")
        g.add_vertex("3")
        g.add_vertex("4")

        g.add_edge(0, 1, 0.2)
        g.add_edge(2, 0, 0.5)
        g.add_edge(3, 2, 0.1)
        g.add_edge(3, 1, 0.3)

        edges = [[0, 1], [0, 2], [2, 3], [1, 3]]
        assert np.allclose(edges, g.edges)
        weights = [0.2, 0.5, 0.1, 0.3]
        assert np.allclose(weights, g.weights)

    def test_add_illegal_edge_1(self):


        g = SimpleWeightedGraph()
        g.add_vertex("1")
        g.add_vertex("2")

        with pytest.raises(ValueError, match="Cannot add edge between non-existing vertices"):
            g.add_edge(2, 0, 0.5)

    def test_add_illegal_edge_2(self):


        g = SimpleWeightedGraph()
        g.add_vertex("1")
        g.add_vertex("2")

        with pytest.raises(ValueError, match="Cannot add edge for a single vertex"):
            g.add_edge(0, 0, 0.5)

