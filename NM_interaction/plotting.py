import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import numpy as np
import math
from .geometryUtils import arrange_rectangular_bars, rotate_points_3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_rectangular_section(ax, canvas, b, h, rebar_coords, event=None):
    """Plot a rectangular section with reinforcement bars."""
    ax.clear()
    rect = Rectangle((0, 0), b, h, facecolor="blue", edgecolor="black", alpha=0.5)
    ax.add_patch(rect)
    ax.set_xlim(0, 1.1 * max(b, h))
    ax.set_ylim(0, 1.1 * max(b, h))
    ax.scatter(rebar_coords[:, 0], rebar_coords[:, 1], color="red")
    ax.grid(True)
    ax.set_axis_on()
    canvas.draw()

def plot_circular_section(ax, canvas, diameter, radial_num_bars, cover, event=None):
    """Plot a circular section with reinforcement bars."""
    ax.clear()
    radius = diameter / 2
    circle = Circle([0, 0], radius, color="blue", alpha=0.5)
    ax.add_patch(circle)
    ax.set_xlim(-1.1 * radius, 1.1 * radius)
    ax.set_ylim(-1.1 * radius, 1.1 * radius)

    if radial_num_bars > 0 and cover > 0:
        theta = 2 * math.pi / radial_num_bars
        bar_radius = radius - cover
        x = [bar_radius * math.sin(theta * i) for i in range(radial_num_bars)]
        y = [bar_radius * math.cos(theta * i) for i in range(radial_num_bars)]
        ax.scatter(x, y, color="red")

    # Set equal aspect ratio to ensure the circle is round
    ax.set_aspect("equal", "box")
    ax.grid(True)
    ax.set_axis_on()
    canvas.draw()

def plot_arbitrary_section(ax, canvas, coord_list, bar_list, event=None):
    """Plot an arbitrary section with reinforcement bars."""
    ax.clear()

    if coord_list:
        vertices_x_vals, vertices_y_vals = zip(*coord_list)
        x_vals_closed = list(vertices_x_vals) + [vertices_x_vals[0]]
        y_vals_closed = list(vertices_y_vals) + [vertices_y_vals[0]]
        ax.plot(vertices_x_vals, vertices_y_vals, linestyle="-", marker="x", color="green")
        ax.fill(x_vals_closed, y_vals_closed, "lightblue", alpha=0.5)

    if bar_list:
        bar_x_vals, bar_y_vals = zip(*bar_list)
        ax.scatter(bar_x_vals, bar_y_vals, marker="o", color="red")

    ax.grid(True)
    ax.set_axis_on()
    canvas.draw()

def plot3Dpoints(section, help=False):
    if help:
        print('''
            Plots 3D interaction diagram for a given column 
            Inputs: section - section object from LSLtools
                    
              Outputs: 3D plot of interaction diagram
                  ''')
    if section.shape == "circular":
        # due to symmetry, only need to plot one half of the interaction diagram and rotate about centroid
        centroid = [0,0,0]
        # Number of symmetrical rotations
        n_rotations = 20
        # retrieve base [M,N] points for circular column from functions in LSL tools
        N_rds_maj_pos, M_rds_maj_pos, N_rds_maj_neg, M_rds_maj_neg, N_rds_min_pos, M_rds_min_pos, N_rds_min_neg, M_rds_min_neg = section.determine_NM_envelope_points()
        base_points = np.column_stack((np.zeros(len(M_rds_maj_pos)),M_rds_maj_pos, N_rds_maj_pos))
        # Store all rotated sets of points
        NM_points_3d = []
        for i in range(n_rotations):
            angle = i * (360 / n_rotations)  # Rotate by equal increments
            rotated_points = rotate_points_3d(base_points, angle)
            NM_points_3d.append(rotated_points)

        # Convert to a NumPy array for easy indexing
        NM_points_3d = np.array(NM_points_3d)
        # Create surface faces (connect corresponding points in successive rotations)
        faces = []
        for i in range(n_rotations):
            p1 = NM_points_3d[i]
            p2 = NM_points_3d[(i + 1) % n_rotations]  # Wrap around for last set
    
        for j in range(len(base_points) - 1):
            # Each quad face is made of 4 points: (p1[j], p1[j+1], p2[j+1], p2[j])
            faces.append([p1[j], p1[j + 1], p2[j + 1], p2[j]])

        # Plot the surface
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Create the 3D surface
        ax.add_collection3d(Poly3DCollection(faces, alpha=0.6, edgecolor='k', facecolor='cyan'))

        # Axes settings
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("3D Rotationally Symmetric Surface Around Z-Axis")
        ax.set_xlim(-6, 6)
        ax.set_ylim(-6, 6)
        ax.set_zlim(-6, 12)
        ax.view_init(elev=30, azim=45)  # Adjust viewing angle

        plt.show()