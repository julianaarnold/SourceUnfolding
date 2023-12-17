import igl
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString

from scipy.spatial import Voronoi, voronoi_plot_2d

from unfolding_utils import *
from visualization import *


EPSILON = 0.00001

def unfold_recursive(current_face, target_face, graph, vertices, faces, blocked_faces):

    if current_face == target_face:
        return [[current_face]]
    

    blocked_faces.append(current_face)
    face_neighbors = list(graph.neighbors(current_face))

    paths = []

    
    
    for neighbor in face_neighbors:
        if neighbor in blocked_faces:
            continue

        neighbor_paths = unfold_recursive(neighbor, target_face, graph, vertices, faces, blocked_faces)

        for path in neighbor_paths:
            paths.append([current_face] + path)

    blocked_faces.pop()

    return paths

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

def clamp_to_edge(edge, point):
    # clamp a point to an edge
    edge_vector = edge[1] - edge[0]
    point_vector = point - edge[0]

    if np.dot(edge_vector, point_vector) < 0:
        return edge[0]

    if np.dot(edge_vector, point_vector) > np.dot(edge_vector, edge_vector):
        return edge[1]

    return point

def is_projection_behind_point(a, proj, source_point):
    # check if a point is behind a projection
    return np.dot(a - source_point, proj - a) > -EPSILON

    
def get_image_on_edge(receiving_edge, edge_to_project, source_point):

    # find the intersection of the receiving edge and the line between the source point and the edge to project
    first_intersection = get_line_intersection([receiving_edge[0], receiving_edge[1]], [edge_to_project[0], source_point])
    second_intersection = get_line_intersection([receiving_edge[0], receiving_edge[1]], [edge_to_project[1], source_point])

    # if there is no intersection, return None
    if first_intersection[0] == float('inf') or second_intersection[0] == float('inf'):
        return None
    
    # ensure that the intersections are behind the edge to project
    if not is_projection_behind_point(edge_to_project[0], first_intersection, source_point):
        return None
    
    if not is_projection_behind_point(edge_to_project[1], second_intersection, source_point):
        return None
    
    # limit the intersection to the edge
    first_intersection = clamp_to_edge(receiving_edge, first_intersection)
    second_intersection = clamp_to_edge(receiving_edge, second_intersection)

    if np.linalg.norm(first_intersection - second_intersection) < EPSILON:
        return None
    
    return [first_intersection, second_intersection]



def unfoldPathCandidate(vertices, faces, path, source_point, projection_to_2d):

    crossed_edges = []
    crossed_polygons = []

    crossed_polygons.append([projection_to_2d.dot(vertices[v]) for v in faces[path[0]]])

    

    print(path)

    current_transform = np.eye(4)
    for i in range(0, len(path)-1):

        crossed_edge_points = [vertices[i] for i in find_common_edge(faces, path[i], path[i+1])]

        offset, rot = get_unfolding_transform(path[i], path[i+1], faces, vertices)
        next_transform = create_unfolding_transform_matrix(offset, rot)
        current_transform = np.dot(current_transform, next_transform)

        transformed_crossed_edge = [projection_to_2d.dot(apply_4x4_matrix_to_3d_point(current_transform, p)) for p in crossed_edge_points]

        crossed_edges.append(transformed_crossed_edge)

        crossed_polygons.append([projection_to_2d.dot(apply_4x4_matrix_to_3d_point(current_transform, vertices[v])) for v in faces[path[i+1]]])

    unfolded_source_point = projection_to_2d.dot(apply_4x4_matrix_to_3d_point(current_transform, source_point))

    
    crossable_segment = None
    for crossed_edge in crossed_edges:
        if crossable_segment is None:
            crossable_segment = crossed_edge
        else:
            crossable_segment = get_image_on_edge(crossable_segment, crossed_edge, unfolded_source_point)

            if crossable_segment is None:
                return None, None, None

    #plt.plot([unfolded_source_point[0]], [unfolded_source_point[1]], 'ro')
    #plt.plot([crossable_segment[0][0], crossable_segment[1][0]], [crossable_segment[0][1], crossable_segment[1][1]], 'r-')
    #plot_polygons(crossed_polygons)
            
    # is the crossable segment the full edge of the face?
    full_edge = True
    if np.linalg.norm(crossable_segment[0] - crossed_edges[0][0]) > EPSILON and np.linalg.norm(crossable_segment[0] - crossed_edges[0][1]) > EPSILON:
        full_edge = False

    if np.linalg.norm(crossable_segment[1] - crossed_edges[0][0]) > EPSILON and np.linalg.norm(crossable_segment[1] - crossed_edges[0][1]) > EPSILON:
        full_edge = False

    return unfolded_source_point, crossable_segment, full_edge

def recursive_calculate_voronoi_tesselation(vertices, faces, face_id, path_data, projection_to_2d, depth = 0, points_to_consider = []):
    print(depth, projection_to_2d)

    if depth == len(path_data):
        """if len(points_to_consider) > 2:
            points_to_consider = np.array(points_to_consider)
            print(points_to_consider)
            vor = Voronoi(points_to_consider)
            voronoi_plot_2d(vor)
        else:
            p1 = points_to_consider[0]
            p2 = points_to_consider[1]

            p = (p1 + p2) / 2
            v = p2 - p1
            v = np.array([-v[1], v[0]])
            v = v * 200  # make the vector long enough to not limit anyting

            mid_segment = np.array([p - v, p + v])
            
            plt.plot(mid_segment[:,0],mid_segment[:,1], ls='--', color='red', linewidth=1.0)"""

        return 1

    p, seg, full_edge = path_data[depth]

    if full_edge:
        return recursive_calculate_voronoi_tesselation(vertices, faces, face_id, path_data, projection_to_2d, depth + 1, points_to_consider + [p])
    else:
        tesselation_with_p = recursive_calculate_voronoi_tesselation(vertices, faces, face_id, path_data, projection_to_2d, depth + 1, points_to_consider + [p])
        tesselation_without_p = recursive_calculate_voronoi_tesselation(vertices, faces, face_id, path_data, projection_to_2d, depth + 1, points_to_consider)

        # do funky stuff here
        merged_tesselation = None

        return tesselation_with_p + tesselation_without_p
        


def apply_source_unfolding(vertices, faces, source_point, source_face_id):
    graph = graph_from_mesh(faces)

    for i in range(len(faces)):
        if i == source_face_id:
            continue

        path_candidates = unfold_recursive(i, source_face_id, graph, vertices, faces, [])
        
        face_normal = igl.per_face_normals(vertices, faces, np.ones((1, 3)))[i]
        projection_to_2d = get_2d_projection(face_normal)

        valid_path_data = []
        for path in path_candidates:
            p, seg, full_edge = unfoldPathCandidate(vertices, faces, path, source_point, projection_to_2d)
            
            if p is not None:
                valid_path_data.append((p, seg, full_edge))

        if len(valid_path_data) == 0:
            print("No valid path found")
            continue

        print("Num_Branches", recursive_calculate_voronoi_tesselation(vertices, faces, i, valid_path_data, projection_to_2d))

        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plot_polygons([[projection_to_2d.dot(vertices[v]) for v in faces[i]]])

        """
        for j in range(len(valid_path_data)):
            for k in range(j):
                p1, seg1 = valid_path_data[j]
                p2, seg2 = valid_path_data[k]

                # construct midpoint line between the two points
                p = (p1 + p2) / 2
                v = p2 - p1
                v = np.array([-v[1], v[0]])
                v = v * 200  # make the vector long enough to not limit anyting

                mid_segment = np.array([p - v, p + v])

                

                # find line segment, where both paths can reach the midpoint line
                mid_segment = get_image_on_edge(seg1, mid_segment, p1)
                if mid_segment is None:
                    continue
                mid_segment = get_image_on_edge(seg2, mid_segment, p2)
                if mid_segment is None:
                    continue

                mid_segment = np.array(mid_segment)
                    
                
                
                
                # find the intersection with the current face
                face = [vertices[v] for v in faces[i]]
                face = Polygon(face)
                mid_segment = LineString(mid_segment)
                intersection = np.array(face.intersection(mid_segment).coords)
                print(intersection)
                # plot the intersection line
                if face.intersects(mid_segment):
                    plt.plot(intersection[:,0],intersection[:,1], ls='--', color='red', linewidth=1.0)
        """
        if len(valid_path_data) >= 2:
            source_points = [p for p, seg, _ in valid_path_data]
            source_points = np.array(source_points)
            vor = Voronoi(source_points)
            voronoi_plot_2d(vor)

            
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        plot_polygons([[projection_to_2d.dot(vertices[v]) for v in faces[i]]])

                




        

    return vertices, faces, []