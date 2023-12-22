import igl
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
from scipy.spatial import Voronoi, voronoi_plot_2d

from basic_unfolding import BasicUnfolding
from star_unfolding import StarUnfolding

from unfolding_utils import *
from visualization import *

class SourceUnfolding(BasicUnfolding):
    def __init__(self, vertices, faces, source_point):
        self.source_point = np.array(source_point)
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
    
    def limit_to_star_unfolding(self, segments, star_unfolding):
        faces_to_unfolded = {}

        print(star_unfolding.face_mapping)

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

        for key, faces in merged_faces.items():
            plot_polygons([face[1] for face in faces])


        intersected_segments = {}
        for segment in segments:
            for key, faces in merged_faces.items():
                for face_id, face in faces:
                    p = Polygon(face)
                    line = LineString(segment)

                    def same_line(line1, line2):
                        return (np.linalg.norm(line1[0] - line2[0]) < 0.0001 and np.linalg.norm(line1[1] - line2[1]) 
                            < 0.0001 or np.linalg.norm(line1[0] - line2[1]) < 0.0001 and np.linalg.norm(line1[1] - line2[0]) < 0.0001)

                    try:
                        if p.intersects(line):
                            intersected_line = np.array(p.intersection(line).coords)

                            if len(intersected_line) == 1:
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

                            """# plot intersected line and polygon it was intersected with
                            plt.plot(intersected_line[:, 0], intersected_line[:, 1], 'r-')
                            plot_polygons([star_unfolding.unfolded_polygons[key]])"""
                    except Exception as e:
                        # do not deal with this exception for now
                        print(e)
                        print("exception occured for face: ", face)
                        print("segment: ", segment)
                        
                        



        return intersected_segments


    def find_cut_locus(self):
        star_unfolding = StarUnfolding(self.vertices, self.faces, self.source_point)

        source_images = star_unfolding.find_source_point_images()

        plt.plot(source_images[:, 0], source_images[:, 1], 'ro')
        plot_polygons(star_unfolding.unfolded_polygons.values())

        # plot voronoi diagram of source images
        vor = Voronoi(source_images)
        
        finite_segments, infinite_segments = self.compute_voronoi_lines(vor)
        segments = np.array(finite_segments + infinite_segments)

        for segment in segments:
            plt.plot(segment[:, 0], segment[:, 1], 'r-')

        plt.show()

        intersected_segments = self.limit_to_star_unfolding(segments, star_unfolding)

        print("intersected segments: ", intersected_segments)

        # plot lines
        for key, value in intersected_segments.items():
            for segment in value:
                plt.plot(segment[:, 0], segment[:, 1], 'r-')

        # set axis limits
        plt.xlim(-6, 6)
        plt.ylim(-6, 6)

        plt.show()

        
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
        

        print("applied transformations: ", star_unfolding.applied_transformations)

        # project unfolded polygons back to 3D
        projected_polygons = []
        for key, value in star_unfolding.unfolded_polygons.items():
            # to fold, use inverse of applied transformation from unfolding
            folding_matrix = star_unfolding.applied_transformations[key]
            folding_matrix = np.linalg.inv(folding_matrix)

            projected_polygon = []
            for point in value:
                projected_polygon.append(apply_4x4_matrix_to_3d_point(folding_matrix, np.append(point, 0)))

            projected_polygons.append(projected_polygon)

        projected_polygons = np.array(projected_polygons)

        # plot the projected segments in 3D
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for segment in self.cut_locus:
            ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], 'r-', linewidth=2.0)
        
        ax.set_xlim(-6, 6)
        ax.set_ylim(-6, 6)
        ax.set_zlim(-6, 6)

        # plot the original mesh
        plot_mesh(self.original_vertices, self.original_faces, ax)

        #plot_polygons_3d(projected_polygons, ax)
        
        plt.show()

        self.cut_locus 