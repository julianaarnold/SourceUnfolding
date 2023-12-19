from unfolding_utils import *

class BasicUnfolding:
    def __init__(self, vertices, faces):
        self.original_vertices = vertices
        self.original_faces = faces
        self.vertices = vertices[:]
        self.faces = faces[:]
        self.faces_to_separate = {}

        self.execute()

    def execute(self):
        self.unfold()

    def unfold(self):
        self.unfolded_polygons = {}  # resulting polygons, represented as lists of 2D coordinates

        source_face_id = 0
        mesh_graph = graph_from_mesh(self.faces)

        mesh_graph = delete_cut_line_edges(mesh_graph, self.faces_to_separate)
        parent_dict = nx.dfs_predecessors(mesh_graph, source=source_face_id)
        parent_dict[source_face_id] = None  # add the source face, as networkX is not doing this by default


        source_face_normal = igl.per_face_normals(self.vertices, self.faces, np.ones((1, 3)))[0]
        projection_to_2d = get_2d_projection(source_face_normal)

        for face_id, parent_face_id in parent_dict.items():
            print("face_id: ", face_id)
            # retrieve the coordinates of current face
            face_coordinates = [self.vertices[vertex_id] for vertex_id in self.faces[face_id]]

            # iterate over all parents and apply unfolding rotations accordingly
            selected_face_id = face_id
            selected_parent_face_id = parent_face_id

            while selected_face_id != source_face_id:
                # get edge between selected face and parent as tuple of two vertex_ids
                hinge_edge = find_common_edge(self.faces, selected_parent_face_id, selected_face_id)
                
                offset = self.vertices[hinge_edge[0]]
                rotation_angle = dihedral_angle(self.vertices, self.faces,  selected_parent_face_id, selected_face_id)
                rotation_axis = self.vertices[hinge_edge[0]] - self.vertices[hinge_edge[1]]
                rotation_matrix = get_rotation_matrix(rotation_axis, rotation_angle)

                for i in range(3):
                    face_coordinates[i] = rotation_matrix.dot(face_coordinates[i] - offset) + offset

                # traverse up the tree
                selected_face_id = selected_parent_face_id
                selected_parent_face_id = parent_dict[selected_parent_face_id]

            # project 3D coordinates into the 2D plane that
            for i in range(3):
                face_coordinates[i] = projection_to_2d.dot(face_coordinates[i])  # TODO apply 2D projection to each vertex

            self.unfolded_polygons[face_id] = face_coordinates