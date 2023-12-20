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
        i, verts, faces = insert_point_into_mesh(self.vertices, self.faces, self.point)
        assert i != -1
        self.source_point_id = i
        self.vertices = verts
        self.faces = faces

        self.cut_star_unfolding()
        self.unfold()
    
    def cut_star_unfolding(self):

        added_vertices = np.zeros((0,3))
        
        paths_indices = []

        self.face_mapping = {}

        num_vertices = len(self.vertices)
        for target_id in range(num_vertices):
            if target_id == self.source_point_id:
                continue
  
            geoalg = geodesic.PyGeodesicAlgorithmExact(self.vertices, self.faces)
            _, path = geoalg.geodesicDistance(self.source_point_id, target_id)

            assert len(path) >= 2

            # store indices of vertices used for cut path
            path_indices = []
            path_indices.append(target_id)

            for v in path[1:-1]:
                self.vertices = np.append(self.vertices, [v], axis = 0)
                cut_edge = find_cut_edge_vertex_ids(self.vertices, self.faces, v)
                cut_faces = find_faces_shared_by_cut_edge(cut_edge, self.faces)
                self.faces, mapping = cut_faces_in_two(self.faces, cut_faces, cut_edge, len(self.vertices) - 1)#

                self.face_mapping.update(mapping)

                path_indices.append(len(self.vertices) - 1)

            path_indices.append(self.source_point_id)
            paths_indices.append(path_indices)

            added_vertices = np.append(added_vertices, path[1:-1], axis = 0)
        

        #now, find the faces next to the cut lines
        self.faces_to_separate = []
        for path in paths_indices:
            for i in range(1, len(path)):
                v_a = path[i-1]
                v_b = path[i]

                self.faces_to_separate.append(find_faces_shared_by_cut_edge([v_a, v_b], self.faces))

    def find_source_point_images(self):
        source_point_images = []

        for face_id, face in enumerate(self.faces):
            if self.source_point_id in face:
                source_point_images.append(self.unfolded_polygons[face_id][np.where(face == self.source_point_id)[0][0]])

        return np.array(source_point_images)             