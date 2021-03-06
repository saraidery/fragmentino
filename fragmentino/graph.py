#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
from copy import deepcopy


class SimpleWeightedGraph:
    """Simple weighted graph class

    Attributes
    ----------
    vertices : list

    edges : numpy.ndarray

    weights : numpy.ndarray

    """

    def __init__(self):
        self.vertices = []
        self.weights = []
        self.edges = []

    def add_vertex(self, vertex):
        """Adds a vertex to the graph
        Parameters
        ----------
        vertex
        """
        self.vertices.append(vertex)

    def add_vertices(self, vertex_list):
        """Adds a list of vertices to the current set of vertices in the graph

        Parameters
        ----------
        vertex_list : list
            list of vertices

        """
        self.vertices.extend(vertex_list)

    def add_edge(self, v1, v2, weight):
        """Adds a weighted edge between two existing vertices (v1 and v2) in the graph

        Parameters
        ----------
        v1 : int
            Index of first vertex
        v2 : int
            Index of second vertex
        weight : float
            Weight of edge
        """
        if v1 >= len(self.vertices) or v2 >= len(self.vertices):
            raise ValueError("Cannot add edge between non-existing vertices")

        if v1 < v2:
            edge = [v1, v2]
        elif v1 > v2:
            edge = [v2, v1]
        else:
            raise ValueError("Cannot add edge for a single vertex")

        self.edges.append(edge)
        self.weights.append(weight)

    @property
    def n_vertices(self):
        """The number of vertices"""
        return len(self.vertices)

    @property
    def n_edges(self):
        """The number of edges"""
        return len(self.edges)

    @property
    def size(self):
        """The size of the graph, defined as the number of edges"""
        return self.n_edges


class ContractableWeightedGraph(SimpleWeightedGraph):
    """Weighted graph that can be contracted by merging vertices
    that are connected by edges with small weights

    Attributes
    ----------
    vertices : list
        list of any type that implements a ``merge`` and a ``size`` method, that
        respectively can merge one instance with another of the same type and
        has some meaningfull size property

    edges : list, numpy.ndarray
        stores pairs of vertex indices representing the edges

    weights : list, numpy.ndarray
        stores the edge weights

    max_vertex_size : int
        Maximal size of vertex

    """

    def __init__(self, max_vertex_size):
        self.vertices = []
        self.weights = []
        self.edges = []
        self._max_vertex_size = max_vertex_size

    def contract_by_smallest_weight(self):
        """
        Contract edges (merge vertices) until no vertices
        can be merged without exceeding the maximal vertex size

        Always try to merge the vertices with the smallest edge weight

        Warning
        -------
        Converts weights and edges to numpy.ndarray and sorts them according
        to ascending weight

        """
        self.weights = np.array(self.weights)
        self.edges = np.array(self.edges)

        self._sort_edges_by_weight()
        edge_index = self._determine_next_graph_contraction()

        while edge_index != -1:
            self._graph_contraction(edge_index)
            edge_index = self._determine_next_graph_contraction()

    def _sort_edges_by_weight(self):
        """Sort edges by weights"""
        self.edges = self.edges[self.weights.argsort()]
        self.weights.sort()

    def _determine_next_graph_contraction(self):
        """Determine next edge contraction

        Returns
        -------
        int
            index for next edge to contract
        """
        for edge_index, edge in enumerate(self.edges):
            if self._can_contract_edge(edge):
                return edge_index
        return -1

    def _can_contract_edge(self, edge):
        """Can contract edge

        Checks if the passed edge can be contracted
        without resulting in a new vertex that exceeds the maximum vertex size

        Parameters
        ----------
        edge : numpy.ndarray

        Returns
        -------
        bool
            True if contraction is possible
        """
        v1, v2 = edge
        new_v_size = self.vertices[v1].size + self.vertices[v2].size
        return new_v_size <= self._max_vertex_size

    def _graph_contraction(self, edge_index):
        """Graph contraction

        Merges the two vertices of the edge (given by edge_index) and
        contracts the edge.

        """
        v1, v2 = self.edges[edge_index]

        self._delete_edge(edge_index)
        self._merge_vertices(v1, v2)
        self._remove_duplicate_edges()

    def _merge_vertices(self, v1, v2):
        """Update vertices"""
        v1_copy = deepcopy(self.vertices[v1])
        v2_copy = deepcopy(self.vertices[v2])

        del self.vertices[v2]
        del self.vertices[v1]

        v1_copy.merge(v2_copy)
        self.vertices.append(v1_copy)

        self._update_vertex_indices_in_edges(v1, v2)

    def _delete_edge(self, edge_index):
        """Delete edge"""
        self.edges = np.delete(self.edges, edge_index, axis=0)
        self.weights = np.delete(self.weights, edge_index, axis=0)

    def _remove_duplicate_edges(self):
        """Remove duplicate edges"""

        self.edges, indices = np.unique(self.edges, axis=0, return_index=True)
        self.weights = self.weights[np.array(indices)]
        self._sort_edges_by_weight()

    def _update_vertex_indices_in_edges(self, v1, v2):
        """Update vertex indices in edges

        Parameters
        ----------
        v1 : int
            index of first vertex that was merged
        v2 : int
            index of second vertex that was merged
        """

        # Rules for updating vertex index (v), when vertex v1 and v2 have been merged
        def _update_vertex_index(v, v1, v2):
            if (v > v1 and v < v2) or (v > v2 and v < v1):
                return v - 1
            elif v > v2 and v > v1:
                return v - 2
            elif v == v2 or v == v1:
                return len(self.vertices) - 1
            else:
                return v

        for edge in self.edges:
            edge[0] = _update_vertex_index(edge[0], v1, v2)
            edge[1] = _update_vertex_index(edge[1], v1, v2)

        self.edges = np.sort(self.edges, axis=1)

    def swap_vertices(self, v1, v2):
        """
        Swaps the order of two vertices
        """
        self.vertices[v1], self.vertices[v2] = self.vertices[v2], self.vertices[v1]

        # Rules for updating vertex index (v), when vertex v1 and v2 have been swapped
        def _update_vertex_index(v, v1, v2):
            if v == v1:
                v = v2
            elif v == v2:
                v = v1
            return v

        for edge in self.edges:
            edge[0] = _update_vertex_index(edge[0], v1, v2)
            edge[1] = _update_vertex_index(edge[1], v1, v2)

        self.edges = np.sort(self.edges, axis=1)
