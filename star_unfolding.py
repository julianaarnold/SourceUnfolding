from basic_unfolding import BasicUnfolding

import numpy as np
from pygeodesic import geodesic

from unfolding_utils import *



class StarUnfolding(BasicUnfolding):
    def __init__(self, verts, faces):
        super().__init__(verts, faces)

    def execute(self):
        self.cut_star_unfolding()
        self.unfold()
    
    def cut_star_unfolding(self):
        geoalg = geodesic.PyGeodesicAlgorithmExact(self.vertices, self.faces)

        added_vertices = np.zeros((0,3))

        source_id = 0

        paths_indices = []

        for target_id in range(1, len(self.vertices)):
            _, path = geoalg.geodesicDistance(source_id, target_id)

            print(path)

            assert len(path) >= 2

            # store indices of vertices used for cut path
            path_indices = []
            path_indices.append(target_id)

            for i in range(len(path)-2):
                path_indices.append(len(self.vertices)+len(added_vertices)+i)

            path_indices.append(source_id)
            paths_indices.append(path_indices)

            added_vertices = np.append(added_vertices, path[1:-1], axis = 0)


        for v in added_vertices:

            self.vertices = np.append(self.vertices, [v], axis = 0)

            cut_edge = find_cut_edge_vertex_ids(self.vertices, self.faces, v)
            cut_faces = find_faces_shared_by_cut_edge(cut_edge, self.faces)
            self.faces = cut_faces_in_two(self.faces, cut_faces, cut_edge, len(self.vertices) - 1)

        #now, find the faces next to the cut lines
        self.faces_to_separate = []
        for path in paths_indices:
            for i in range(1, len(path)):
                v_a = path[i-1]
                v_b = path[i]

                self.faces_to_separate.append(find_faces_shared_by_cut_edge([v_a, v_b], self.faces))

        