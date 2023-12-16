import igl
from meshplot import plot
import numpy as np
import pygeodesic.geodesic as geodesic
import networkx as nx
import drawsvg as draw

raw_vertices, raw_faces = igl.read_triangle_mesh("./meshes/cube.stl")
vertices, faces = igl.remove_duplicates(raw_vertices, raw_faces, 0.00001)

print(vertices, faces)