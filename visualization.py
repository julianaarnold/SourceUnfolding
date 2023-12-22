import numpy as np
import matplotlib.pyplot as plt

def plot_polygons(polygons):
    for polygon in polygons:
        polygon = np.array(polygon)

        plt.fill(polygon[:, 0], polygon[:, 1], 'b', alpha=0.5)

        for i in range(len(polygon)):
            plt.plot([polygon[i][0], polygon[(i+1) % len(polygon)][0]], [polygon[i][1], polygon[(i+1) % len(polygon)][1]], ls='-', color='black', linewidth=1.0)


def plot_cut_edges(polygons, faces, faces_to_separate):
    # find the shared vertices foreach [f1, f2] in faces_to_separate

    edges_to_plot = []

    for f1, f2 in faces_to_separate:
        shared_vertices = np.intersect1d(faces[f1], faces[f2])

        # find indices of shared vertices in f1 and f2
        v1_f1 = np.where(faces[f1] == shared_vertices[0])[0][0]
        v2_f1 = np.where(faces[f1] == shared_vertices[1])[0][0]

        v1_f2 = np.where(faces[f2] == shared_vertices[0])[0][0]
        v2_f2 = np.where(faces[f2] == shared_vertices[1])[0][0]

        edges_to_plot.append([polygons[f1][v1_f1], polygons[f1][v2_f1]])
        edges_to_plot.append([polygons[f2][v1_f2], polygons[f2][v2_f2]])

    for edge in edges_to_plot:
        plt.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], ls='-', color='red', linewidth=2.0)

    

def plot_polygons_3d(polygons, ax):
    for polygon in polygons:
        polygon = np.array(polygon)

        # plot lines in 3D
        for i in range(len(polygon)):
            ax.plot([polygon[i][0], polygon[(i+1) % len(polygon)][0]], [polygon[i][1], polygon[(i+1) % len(polygon)][1]], [polygon[i][2], polygon[(i+1) % len(polygon)][2]], ls='-', color='black', linewidth=1.0)



def plot_mesh(vertices, faces, ax):
    polygons = [vertices[face] for face in faces]
    plot_polygons_3d(polygons, ax)