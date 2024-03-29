import igl
import numpy as np

from unfolding_utils import *
from visualization import *

from basic_unfolding import BasicUnfolding
from star_unfolding import StarUnfolding
from source_unfolding import SourceUnfolding

raw_vertices, raw_faces = igl.read_triangle_mesh("./meshes/pyramid.stl")
vertices, faces = igl.remove_duplicates(raw_vertices, raw_faces, EPSILON)

unfolding = SourceUnfolding(vertices, faces, [0.5, 0.0, 0.0], show_intermediate_results=True, report_errors=True)

plot_polygons(unfolding.unfolded_polygons.values())
plot_cut_edges(unfolding.unfolded_polygons, unfolding.faces, unfolding.faces_to_separate)

plt.gca().set_aspect('equal', adjustable='box')
plt.show()