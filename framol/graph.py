import numpy as np
from copy import deepcopy


class SimpleWeightedGraph:
    """Simple weighted graph class

    Attributes
    ----------
    vertices : list
        list of any type
    max_vertex_size : int
        Maximal size of vertex

    """

    def __init__(self, vertices):
        self.vertices = vertices
        self.weights = []
        self.edges = []

    def add_edge(self, v1, v2, weight):

        if v1 >= len(self.vertices) or v2 >= len(self.vertices):
            raise ValueError("Cannot add edge between non-existing vertices")

        self.edges.append([v1, v2])
        self.weights.append(weight)


class WeightedGraph(SimpleWeightedGraph):
    """Weighted graph class

    Attributes
    ----------
    vertices : list
        list of any type that implements a merge and a size method, that
        respectively can merge one instance with another of the same type and
        has some meaningfull size property
    max_vertex_size : int
        Maximal size of vertex

    """

    def __init__(self, vertices, max_vertex_size):
        self.vertices = vertices
        self.weights = []
        self.edges = []
        self._max_vertex_size = max_vertex_size

    def contract_by_smallest_weight(self):
        """Contract by smallest weights

        Contract edges (merge vertices) until no vertices
        can be merged without exceeding the maximal vertex size

        """

        self._sort_edges_by_weight()
        edge_index = self._determine_next_edge_contraction()

        while edge_index != -1:

            self._edge_contraction(edge_index)
            edge_index = self._determine_next_edge_contraction()


    def _sort_edges_by_weight(self):
        """Sort edges by weights"""

        self.weights = np.array(self.weights)
        self.edges = np.array(self.edges)

        self.edges = self.edges[self.weights.argsort()]
        self.weights.sort()

    def _can_merge_vertices(self, edge):
        """Can merge vertices

        Parameters
        ----------

        Returns
        -------
        bool
            True if a merge of vertex v1 and v2 will not result in a vertex that exceeds the maximal size
        """
        v1, v2 = edge
        return (self.vertices[v1].size + self.vertices[v2].size) < self._max_vertex_size

    def _determine_next_edge_contraction(self):
        """Determine next edge contraction

        Returns
        -------
        int
            index for next edge to contract
        """

        for edge_index, edge in enumerate(self.edges):
            if self._can_merge_vertices(edge):
                return edge_index
        return -1

    def _edge_contraction(self, edge_index):
        """Edge contraction"""

        self._update_vertices(edge_index)
        self._update_edges(edge_index)

    def _update_vertices(self, edge_index):
        """Update vertices inplace"""

        v1, v2 = self.edges[edge_index]

        # Copy the vertices that are to be merged
        v1_copy = deepcopy(self.vertices[v1])
        v2_copy = deepcopy(self.vertices[v2])

        # Delete vertices that are to be merged
        if v1 > v2:
            del self.vertices[v1]
            del self.vertices[v2]
        if v2 > v1:
            del self.vertices[v2]
            del self.vertices[v1]

        # Merge vertices and add to list
        v1_copy.merge(v2_copy)
        self.vertices.append(v1_copy)


    def _update_edges(self, edge_index):
        """Update edges inplace"""

        v1, v2 = self.edges[edge_index]

        self._contract_edge(edge_index)
        self._remove_duplicate_edges(v1, v2)

        for edge in self.edges:
            self._update_vertex_indices_in_edges(edge, v1, v2)


    def _contract_edge(self, edge_index):
        self.edges = np.delete(self.edges, edge_index, axis=0)
        self.weights = np.delete(self.weights, edge_index, axis=0)

    def _remove_duplicate_edges(self, v1, v2):
        found_indices = []
        delete_indices = []

        def _update_found_and_delete(j, v, found_indices, delete_indices):
            if v in found_indices:
                delete_indices.append(j)
            else:
                found_indices.append(v)

        for j, edge in enumerate(self.edges):

            if edge[0] == v1 or edge[0] == v2:
                _update_found_and_delete(j, edge[1], found_indices, delete_indices)

            elif edge[1] == v1 or edge[1] == v2:
                _update_found_and_delete(j, edge[0], found_indices, delete_indices)

        self.edges = np.delete(self.edges, delete_indices, axis=0)
        self.weights = np.delete(self.weights, delete_indices, axis=0)

    def _update_vertex_indices_in_edges(self, edge, v1, v2):
        """Update vertex indices in edges inplace"""

        def _update_vertex_index(v, v1, v2):
            if (v > v1 and v < v2) or (v > v2 and v < v1):
                return v - 1
            elif v > v2 and v > v1:
                return v - 2
            elif v == v2 or v == v1:
                return len(self.vertices) - 1
            else:
                return v

        edge[0] = _update_vertex_index(edge[0], v1, v2)
        edge[1] = _update_vertex_index(edge[1], v1, v2)
