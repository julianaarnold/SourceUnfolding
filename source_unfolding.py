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
        pass

    def find_cut_locus(self):
        star_unfolding = StarUnfolding(self.vertices, self.faces, self.source_point)

        source_images = star_unfolding.find_source_point_images()

        plt.plot(source_images[:, 0], source_images[:, 1], 'ro')
        plot_polygons(star_unfolding.unfolded_polygons.values())

        # plot voronoi diagram of source images
        vor = Voronoi(source_images)
        
        finite_segments, infinite_segments = self.compute_voronoi_lines(vor)
        segments = np.append(np.array(finite_segments), np.array(infinite_segments), axis = 0)

        # plot lines
        for segment in segments:
            plt.plot(segment[:, 0], segment[:, 1], 'r-')

        # set axis limits
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)

        print(finite_segments)
        print(infinite_segments)

        plt.show()