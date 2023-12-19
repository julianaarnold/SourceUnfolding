from basic_unfolding import BasicUnfolding

import numpy as np
from pygeodesic import geodesic

from unfolding_utils import *

from visualization import *
import matplotlib.pyplot as plt

class StarUnfolding(BasicUnfolding):
    def __init__(self, verts, faces, point):
        self.point = np.array(point)
        super().__init__(verts, faces)

    def execute(self):
        #i, verts, faces = insert_point_into_mesh(self.vertices, self.faces, self.point)
        #assert i != -1
        self.source_face_id = 3#i
        #print("Source face id: ", self.source_face_id)
        #self.vertices = verts
        #self.faces = faces

        print("vertices: ", self.vertices, "faces: ", self.faces)

        self.cut_star_unfolding()
        self.unfold()
    
    def cut_star_unfolding(self):
        geoalg = geodesic.PyGeodesicAlgorithmExact(self.vertices, self.faces)

        added_vertices = np.zeros((0,3))
        
        paths_indices = []

        for target_id in range(len(self.vertices)):
            if target_id == self.source_face_id:
                continue

            _, path = geoalg.geodesicDistance(self.source_face_id, target_id)

            print(path)

            assert len(path) >= 2

            # store indices of vertices used for cut path
            path_indices = []
            path_indices.append(target_id)

            for i in range(len(path)-2):
                path_indices.append(len(self.vertices)+len(added_vertices)+i)

            path_indices.append(self.source_face_id)
            paths_indices.append(path_indices)

            added_vertices = np.append(added_vertices, path[1:-1], axis = 0)


        for v in added_vertices:

            self.vertices = np.append(self.vertices, [v], axis = 0)

            print("adding vertex: ", v, "with id: ", len(self.vertices) - 1)

            cut_edge = find_cut_edge_vertex_ids(self.vertices, self.faces, v)
            print("cut edge: ", cut_edge)
            cut_faces = find_faces_shared_by_cut_edge(cut_edge, self.faces)
            print("faces before cut: ", self.faces)
            self.faces = cut_faces_in_two(self.faces, cut_faces, cut_edge, len(self.vertices) - 1)
            print("faces after cut: ", self.faces)

            if len(self.vertices) == 12:
                self.unfold()
                plot_polygons(self.unfolded_polygons.values())
                plot_cut_edges(self.unfolded_polygons, self.faces, [[len(self.faces)-3, len(self.faces)-4], [len(self.faces)-1, len(self.faces)-2]])
                plt.show()

        #now, find the faces next to the cut lines
        self.faces_to_separate = []
        for path in paths_indices:
            print("path: ", path)
            for i in range(1, len(path)):
                v_a = path[i-1]
                v_b = path[i]

                print("v_a: ", v_a, "v_b: ", v_b)
                self.faces_to_separate.append(find_faces_shared_by_cut_edge([v_a, v_b], self.faces))

        