import sys
import igl
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
from scipy.spatial import Voronoi, voronoi_plot_2d

from pygeodesic import geodesic

from basic_unfolding import BasicUnfolding
from star_unfolding import StarUnfolding

from unfolding_utils import *
from visualization import *

class SourceUnfolding(BasicUnfolding):
    def __init__(self, vertices, faces, source_point, show_intermediate_results=False, report_errors=False):
        self.source_point = np.array(source_point)
        self.debug = show_intermediate_results
        self.report_errors = report_errors
        super().__init__(vertices, faces)

    def execute(self):
        i, verts, faces = insert_point_into_mesh(self.vertices, self.faces, self.source_point)

        assert i != -1
        self.source_face_id = i
        self.vertices = verts
        self.faces = faces

        self.cut_source_unfolding()
        self.unfold()

    def cut_source_unfolding(self):
        self.find_cut_locus()

        self.faces_to_separate = []

        cut_paths = []
        for segment in self.cut_locus:
            edge_verts = []
            for p in segment:
                i, verts, faces = insert_point_into_mesh(self.vertices, self.faces, p)

                assert i != -1
                self.vertices = verts
                self.faces = faces

                edge_verts.append(i)

            # there could be some edge between the two vertices
            geoalg = geodesic.PyGeodesicAlgorithmExact(self.vertices, self.faces)
            _, path = geoalg.geodesicDistance(edge_verts[0], edge_verts[1])

            if len(path) <= 1:
                continue

            if len(path) > 2:
                path_indices = []
                path_indices.append(edge_verts[1])
                # cut the edge
                for v in path[1:-1]:
                    i, self.vertices, self.faces = insert_point_into_mesh(self.vertices, self.faces, v)
                    path_indices.append(len(self.vertices) - 1)

                path_indices.append(edge_verts[0])
                edge_verts = path_indices

            cut_paths.append(edge_verts)

        for path in cut_paths:
            # find faces shared by the two vertices
            for i in range(len(path) - 1):
                try:
                    cut_faces = find_faces_shared_by_cut_edge([path[i], path[i+1]], self.faces)
                    self.faces_to_separate.append(cut_faces)
                except:
                    if (self.debug):
                        print(path)

    def compute_voronoi_lines(self, vor):
        # copied from scipy.spatial.voronoi_plot_2d

        center = vor.points.mean(axis=0)
        ptp_bound = vor.points.ptp(axis=0)

        finite_segments = []
        infinite_segments = []
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            if np.all(simplex >= 0):
                finite_segments.append(vor.vertices[simplex])
            else:
                i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

                t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal

                midpoint = vor.points[pointidx].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                if (vor.furthest_site):
                    direction = -direction
                far_point = vor.vertices[i] + direction * ptp_bound.max()

                infinite_segments.append([vor.vertices[i], far_point])

        return finite_segments, infinite_segments
    
    def _debug_plot_point_images_with_voronoi(self, source_images, star_unfolding, finite_segments, infinite_segments):
        if not self.debug:
            return
        
        plt.plot(source_images[:, 0], source_images[:, 1], 'ro')
        plot_polygons(star_unfolding.unfolded_polygons.values())

        # plot finite segments as solid lines and infinite segments as dashed lines
        for segment in np.array(finite_segments):
            plt.plot(segment[:, 0], segment[:, 1], 'r-')

        for segment in np.array(infinite_segments):
            plt.plot(segment[:, 0], segment[:, 1], 'r--')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def _debug_plot_cut_locus_3d(self):
        if not self.debug:
            return
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for segment in self.cut_locus:
            ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], 'r-', linewidth=5.0)

        plt.gca().set_aspect('equal', adjustable='box')

        # plot the original mesh
        plot_mesh(self.original_vertices, self.original_faces, ax)
        plt.show()

    
    def limit_to_star_unfolding(self, segments, star_unfolding):
        faces_to_unfolded = {}

        for key, value in star_unfolding.unfolded_polygons.items():
            original_face = star_unfolding.face_mapping[key]

            if not original_face in faces_to_unfolded.keys():
                faces_to_unfolded[original_face] = []

            faces_to_unfolded[original_face].append((key, value))

        merged_faces = {}
        for key, value in faces_to_unfolded.items():
            faces_to_merge = value[:]

            while len(faces_to_merge) > 0:
                # keep track of the face id to be able to refold the mesh later
                face_id, merged_face = faces_to_merge.pop()

                # iterate through all faces to merge in reverse order
                for i in range(len(faces_to_merge) - 1, -1, -1):
                    face = faces_to_merge[i][1]

                    merged_face_candidate = try_merge_polygons_2d(merged_face, face)

                    if merged_face_candidate is not None:
                        merged_face = merged_face_candidate
                        faces_to_merge.pop(i)

                if not key in merged_faces.keys():
                    merged_faces[key] = []

                merged_faces[key].append((face_id, merged_face))

        intersected_segments = {}
        for segment in segments:
            for key, faces in merged_faces.items():
                for face_id, face in faces:
                    p = Polygon(face).buffer(EPSILON/100)
                    line = LineString(segment)

                    def same_line(line1, line2):
                        return (np.linalg.norm(line1[0] - line2[0]) < EPSILON and np.linalg.norm(line1[1] - line2[1]) 
                            < EPSILON or np.linalg.norm(line1[0] - line2[1]) < EPSILON and np.linalg.norm(line1[1] - line2[0]) < EPSILON)

                    try:
                        if p.intersects(line):
                            intersected_geometry = p.intersection(line)
                            intersected_line = np.array(intersected_geometry.coords)
                            
                            if len(intersected_line) == 1:
                                continue

                            # check if the intersection is a point
                            if intersected_geometry.length < EPSILON:
                                continue

                            # check for duplicate lines
                            duplicate = False
                            for key1, value in intersected_segments.items():
                                for seg in value:
                                    if same_line(seg, intersected_line):
                                        duplicate = True
                                        break

                                if duplicate:
                                    break

                            if duplicate:
                                continue

                            if not face_id in intersected_segments.keys():
                                intersected_segments[face_id] = []

                            intersected_segments[face_id].append(intersected_line)
                    except Exception as e:
                        # do not deal with this exception for now
                        if self.report_errors:
                            print(e)
                            print("exception occured for face: ", face)
                            print("segment: ", segment)
                        

        if self.debug:
            plt.gca().set_aspect('equal', adjustable='box')
            explode_polygons_with_intersected_segments(merged_faces, intersected_segments)

        return intersected_segments


    def find_cut_locus(self):
        star_unfolding = StarUnfolding(self.vertices, self.faces, self.source_point)

        source_images = star_unfolding.find_source_point_images()

        # plot voronoi diagram of source images
        vor = Voronoi(source_images)
        
        finite_segments, infinite_segments = self.compute_voronoi_lines(vor)
        segments = np.array(finite_segments + infinite_segments)

        self._debug_plot_point_images_with_voronoi(source_images, star_unfolding, finite_segments, infinite_segments)

        intersected_segments = self.limit_to_star_unfolding(segments, star_unfolding)
        
        # project the segments to 3D
        projected_segments = []
        for key, value in intersected_segments.items():
            # to fold, use inverse of applied transformation from unfolding
            folding_matrix = star_unfolding.applied_transformations[key]
            folding_matrix = np.linalg.inv(folding_matrix)

            for segment in value:
                projected_segment = []
                for point in segment:
                    projected_segment.append(apply_4x4_matrix_to_3d_point(folding_matrix, np.append(point, 0)))

                projected_segments.append(projected_segment)

        self.cut_locus = np.array(projected_segments)

        self._debug_plot_cut_locus_3d()
