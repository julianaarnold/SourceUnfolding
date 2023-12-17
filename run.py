import igl
import numpy as np

from unfolding_utils import *
from visualization import *
from source_unfolding import *

raw_vertices, raw_faces = igl.read_triangle_mesh("./meshes/cube.stl")
vertices, faces = igl.remove_duplicates(raw_vertices, raw_faces, 0.00001)
vertices, faces, faces_to_separate = apply_source_unfolding(vertices, faces, np.array([0, 0, 1]), 0)


polygons = unfold(vertices, faces, faces_to_separate=faces_to_separate)
plot_polygons(polygons)