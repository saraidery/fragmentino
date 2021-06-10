#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
import numpy as np
import pytest
import os


from fragmentino import SimpleWeightedGraph


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

        with pytest.raises(
            ValueError, match="Cannot add edge between non-existing vertices"
        ):
            g.add_edge(2, 0, 0.5)

    def test_add_illegal_edge_2(self):
        g = SimpleWeightedGraph()
        g.add_vertex("1")
        g.add_vertex("2")

        with pytest.raises(ValueError, match="Cannot add edge for a single vertex"):
            g.add_edge(0, 0, 0.5)

    def test_n_vertices(self):
        g = SimpleWeightedGraph()
        g.add_vertex("1")
        g.add_vertex("2")
        g.add_vertex("3")
        g.add_vertex("4")

        g.add_edge(0, 1, 0.2)
        g.add_edge(2, 0, 0.5)
        g.add_edge(3, 2, 0.1)
        g.add_edge(3, 1, 0.3)

        assert g.n_vertices == 4

    def test_n_edges(self):
        g = SimpleWeightedGraph()
        g.add_vertices(["1", "2", "3", "4"])

        g.add_edge(0, 1, 0.2)
        g.add_edge(2, 0, 0.5)
        g.add_edge(3, 2, 0.1)

        assert g.n_edges == 3

    def test_size(self):
        g = SimpleWeightedGraph()
        g.add_vertex("1")
        g.add_vertex("2")
        g.add_vertex("3")
        g.add_vertex("4")

        g.add_edge(0, 1, 0.2)
        g.add_edge(2, 0, 0.5)
        g.add_edge(3, 2, 0.1)
        g.add_edge(3, 1, 0.3)

        assert g.size == 4
