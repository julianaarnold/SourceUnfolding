import numpy as np
import matplotlib.pyplot as plt

def plot_polygons(polygons):
    for polygon in polygons:
        polygon = np.array(polygon)

        plt.fill(polygon[:, 0], polygon[:, 1], 'b', alpha=0.5)

        for i in range(len(polygon)):
            plt.plot([polygon[i][0], polygon[(i+1) % len(polygon)][0]], [polygon[i][1], polygon[(i+1) % len(polygon)][1]], ls='-', color='black', linewidth=1.0)
    
    plt.show()