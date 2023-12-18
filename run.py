import igl
import numpy as np

from unfolding_utils import *
from visualization import *
from source_unfolding import *

from basic_unfolding import BasicUnfolding

raw_vertices, raw_faces = igl.read_triangle_mesh("./meshes/cube.stl")
vertices, faces = igl.remove_duplicates(raw_vertices, raw_faces, 0.00001)

unfolding = BasicUnfolding(vertices, faces)
unfolding.unfold()

plot_polygons(unfolding.unfolded_polygons)