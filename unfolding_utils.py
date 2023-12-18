import igl
from meshplot import plot
import numpy as np
import networkx as nx
import drawsvg as draw
from scipy.spatial.transform import Rotation

def is_vertex_on_line_segment(v, a, b):
  if np.linalg.norm(np.cross(b-a, v-a)) > 0.01:
      return False

  AVdotAB = (v - a).dot(b - a)
  ABdotAB = (b - a).dot(b - a)

  return (AVdotAB > 0) and (AVdotAB < ABdotAB)

def find_cut_edge_vertex_ids(vertices, faces, new_vert):
  edges = igl.edges(faces)

  for e in edges:
    if is_vertex_on_line_segment(new_vert, vertices[e[0]], vertices[e[1]]):
      return e

  return None

def find_faces_shared_by_cut_edge(cut_edge, faces):
  shared_faces = []

  for i in range(len(faces)):
    f = faces[i]
    if cut_edge[0] in f and cut_edge[1] in f:
      shared_faces.append(i)

  assert len(shared_faces) == 2

  return np.array(shared_faces)

def get_cut_faces(faces, face_id, cut_edge_vertices, new_vert_id):
  face_to_cut = faces[face_id]

  first_half = np.copy(face_to_cut)
  second_half = np.copy(face_to_cut)

  first_half[np.nonzero(first_half == cut_edge_vertices[0])[0][0]] = new_vert_id
  second_half[np.nonzero(second_half == cut_edge_vertices[1])[0][0]] = new_vert_id

  return np.array([first_half, second_half])

def cut_faces_in_two(faces, shared_face_ids, cut_edge_vertices, new_vert_id):
  cut_a = get_cut_faces(faces, shared_face_ids[0], cut_edge_vertices, new_vert_id)
  cut_b = get_cut_faces(faces, shared_face_ids[1], cut_edge_vertices, new_vert_id)
  faces = np.delete(faces, list(shared_face_ids), 0)

  faces = np.append(faces, cut_a, axis=0)
  faces = np.append(faces, cut_b, axis=0)

  return faces

def add_edges_from_mesh(graph, faces):
    face_adjacency_matrix, _ = igl.triangle_triangle_adjacency(faces)
    for face_id, _ in enumerate(faces):
        for ajd_face_id in face_adjacency_matrix[face_id]:
            graph.add_edge(face_id,ajd_face_id)

# adds nodes to graph G for every face in mesh
def add_nodes_from_mesh(graph, faces): [graph.add_node(face_id) for face_id, face in enumerate(faces)]

# create networkx graph from given mesh
def graph_from_mesh(faces):
    graph = nx.Graph()
    add_nodes_from_mesh(graph, faces)
    add_edges_from_mesh(graph, faces)
    return graph

# returns a rotation matrix from a (unnormalized) axis and an angle
def get_rotation_matrix(axis, angle):
    if np.linalg.norm(axis) == 0:
        return np.eye(3)
    
    return Rotation.from_rotvec(axis/np.linalg.norm(axis) * angle).as_matrix()

# returns a matrix that maps 3D space onto a 2D plane (the orientation of which is specified by 'face_normal').
def get_2d_projection(face_normal):
    xy_plane_normal = np.array([0,0,1])  # aka 'the z-axis'
    rotation_axis = np.cross(face_normal, xy_plane_normal)
    angle = np.arccos(np.clip(np.dot(xy_plane_normal, face_normal), -1.0, 1.0))

    discard_z_matrix = np.array([
        [1, 0, 0],
        [0, 1, 0]
    ])

    rotation_matrix = get_rotation_matrix(rotation_axis, angle)
    return discard_z_matrix.dot(rotation_matrix)


def draw_polygons(polygons):
    # generate svg visualization
    drawing = None
    drawing = draw.Drawing(1000, 300, origin='center')
    for polygon in polygons:
        # polygon = [coords[0:2] for coords in polygon]
        scale_factor = 1
        drawing.append(draw.Lines(*np.array(polygon).flatten()*scale_factor,
                                  close=True, fill='#eeee00', stroke='#000', stroke_width=.1))

    drawing.rasterize()
    return drawing

# find commone edge of two adjacent faces
def find_common_edge(faces, face_id_a, face_id_b):
    # make sure that the resulting vertex ids are clockwise wrt. source face
    face_vertex_array_a = faces[face_id_a]
    face_vertex_array_b = faces[face_id_b]

    for i in range(3):
        if face_vertex_array_a[i] in face_vertex_array_b and face_vertex_array_a[(i+1) % 3] in face_vertex_array_b:
            return (face_vertex_array_a[i], face_vertex_array_a[(i+1) % 3])

    return None

def get_face_normal(vertices, faces, face_id):
    face_normals = igl.per_face_normals(vertices, faces, np.ones((1, 3)))
    return face_normals[face_id]

# get angle between the normals of two faces
def dihedral_angle(vertices, faces, face_a_id, face_b_id):
    face_normals = igl.per_face_normals(vertices, faces, np.ones((1, 3)))
    return np.arccos(np.clip(np.dot(face_normals[face_a_id], face_normals[face_b_id]), -1.0, 1.0))

def delete_cut_line_edges(graph, faces_to_separate):
  for pair in faces_to_separate:
    graph.remove_edge(pair[0], pair[1])

  return graph

def get_unfolding_transform(source_face_id, face_to_unfold_id, faces, vertices):
    hinge_edge = find_common_edge(faces, source_face_id, face_to_unfold_id)

    offset = vertices[hinge_edge[0]]
    rotation_angle = dihedral_angle(vertices, faces,  source_face_id, face_to_unfold_id)
    rotation_axis = vertices[hinge_edge[0]] - vertices[hinge_edge[1]]
    rotation_matrix = get_rotation_matrix(rotation_axis, rotation_angle)

    return offset, rotation_matrix

def create_unfolding_transform_matrix(offset, rotation_matrix):
    # get the matrix that computes (p - offset) * rotation_matrix + offset as 4x4 matrix

    transform_matrix = np.eye(4)
    transform_matrix[:3, :3] = rotation_matrix
    transform_matrix[:3, 3] = offset + np.dot(rotation_matrix, -offset)
    
    return transform_matrix

def apply_4x4_matrix_to_3d_point(matrix, point):
    # apply 4x4 matrix to 3d points
    point = np.append(point, 1)
    point = np.dot(matrix, point)
    point = point[:3]
    return point


def get_line_intersection(line1, line2):
    print(line1, line2)
    a1 = np.array(line1[0])
    a2 = np.array(line1[1])
    b1 = np.array(line2[0])
    b2 = np.array(line2[1])

    # see https://stackoverflow.com/a/42727584
    s = np.vstack([a1,a2,b1,b2])        # s for stacked
    h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
    l1 = np.cross(h[0], h[1])           # get first line
    l2 = np.cross(h[2], h[3])           # get second line
    x, y, z = np.cross(l1, l2)          # point of intersection
    if z == 0:                          # lines are parallel
        return (float('inf'), float('inf'))
    return np.array([x/z, y/z])


def unfold(vertices, faces, faces_to_separate=[]):
    polygons = []  # resulting polygons, represented as lists of 2D coordinates

    source_face_id = 0
    mesh_graph = graph_from_mesh(faces)
    mesh_graph = delete_cut_line_edges(mesh_graph, faces_to_separate)
    parent_dict = nx.dfs_predecessors(mesh_graph, source=source_face_id)
    parent_dict[source_face_id] = None  # add the source face, as networkX is not doing this by default


    source_face_normal = igl.per_face_normals(vertices, faces, np.ones((1, 3)))[0]
    projection_to_2d = get_2d_projection(source_face_normal)

    for face_id, parent_face_id in parent_dict.items():
        # retrieve the coordinates of current face
        face_coordinates = [vertices[vertex_id] for vertex_id in faces[face_id]]

        # iterate over all parents and apply unfolding rotations accordingly
        selected_face_id = face_id
        selected_parent_face_id = parent_face_id

        while selected_face_id != source_face_id:
            # get edge between selected face and parent as tuple of two vertex_ids
            hinge_edge = find_common_edge(faces, selected_parent_face_id, selected_face_id)
            
            offset = vertices[hinge_edge[0]]
            rotation_angle = dihedral_angle(vertices, faces,  selected_parent_face_id, selected_face_id)
            rotation_axis = vertices[hinge_edge[0]] - vertices[hinge_edge[1]]
            rotation_matrix = get_rotation_matrix(rotation_axis, rotation_angle)

            for i in range(3):
                face_coordinates[i] = rotation_matrix.dot(face_coordinates[i] - offset) + offset

            # traverse up the tree
            selected_face_id = selected_parent_face_id
            selected_parent_face_id = parent_dict[selected_parent_face_id]

        # project 3D coordinates into the 2D plane that
        for i in range(3):
            face_coordinates[i] = projection_to_2d.dot(face_coordinates[i])  # TODO apply 2D projection to each vertex

        polygons.append(face_coordinates)

    return polygons