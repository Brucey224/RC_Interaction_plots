import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import numpy as np
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .resistanceCalcs import determine_envelope_value_major_axis_negative, determine_envelope_value_major_axis_positive, determine_envelope_value_minor_axis_negative, determine_envelope_value_minor_axis_positive
from .geometryUtils import arrange_rectangular_bars, arrange_circular_bars, rotate_points_3d

def plot_rectangular_section(ax, canvas, b_input, h_input, cover_input, link_dia_input, bar_dia_input, n_x_input, n_y_input, event=None):
    ## Plot a rectangular section with reinforcement bars
    try:
        b = float(b_input.get())
        h = float(h_input.get())
        cover = float(cover_input.get())
        link_dia = float(link_dia_input.get())
        bar_dia = float(bar_dia_input.get())
        n_x = int(n_x_input.get())
        n_y = int(n_y_input.get())
        ax.clear()
        rebar_coords = arrange_rectangular_bars(b, h, cover, link_dia, bar_dia, n_x, n_y)
        rect = Rectangle((0, 0), b, h, facecolor="blue", edgecolor="black", alpha=0.5)
        ax.add_patch(rect)
        ax.set_xlim(0, 1.1 * max(b, h))
        ax.set_ylim(0, 1.1 * max(b, h))
        ax.scatter(rebar_coords[:, 0], rebar_coords[:, 1], color="red")
        ax.grid(True)
        ax.set_axis_on()
        canvas.draw()
    except ValueError as e:
        print(f"Invalid Input: {e}")

def plot_circular_section(ax, canvas, diameter, radial_num_bars, cover, event=None):
    ## Plot a circular section with reinforcement bars
    ax.clear()
    radius = diameter / 2
    circle = Circle([0, 0], radius, color="blue", alpha=0.5)
    ax.add_patch(circle)
    ax.set_xlim(-1.1 * radius, 1.1 * radius)
    ax.set_ylim(-1.1 * radius, 1.1 * radius)

    rebar_coords = arrange_circular_bars(diameter, cover, radial_num_bars)
    ax.scatter(rebar_coords[:, 0], rebar_coords[:, 1], color="red")

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

def plot_major_axis_failure_envelope(ax, canvas, column, N_Ed, M_Edy, My_02):
    #Collect user inputs for processing
    N_Rd_positive_list = []
    N_Rd_negative_list = []
    M_Rdy_positive_list = []
    M_Rdy_negative_list = []
    Npl_Rd = (column.section.A * column.concrete_properties.f_cd + len(column.reinforcement.arrangement) * math.pi * column.reinforcement.bar_diameter**2 / 4 * column.reinforcement.f_yd)*1e-3

    if column.section.shape == "rectangular" or column.section.shape == "arbitrary":
        y_limit = column.section.h
    elif column.section.shape == "circular":
        y_limit = column.section.diameter
    for y in range(1, int(y_limit*5),1):
        #add axial force and major axis moment to list for plotting
        N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, y)
        N_Rd_positive_list.append(min(Npl_Rd,N_Rd))
        M_Rdy_positive_list.append(M_Rdy)
        
    # PLOT NEGATIVE MOMENT SIDE OF DIAGRAM
    for y in range(1, int(y_limit*5), 1):
        # add axial force and major axis moment to list for plotting
        N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_negative(column, 0.8, y)
        N_Rd_negative_list.append(min(Npl_Rd,N_Rd))
        M_Rdy_negative_list.append(M_Rdy)
        
    N_ratio = N_Ed / Npl_Rd

    if N_ratio <= 0.1:
        a = 1.0
    elif N_ratio < 0.7:
        a = 1.0 + (N_ratio - 0.1) * 0.5 / 0.6
    elif N_ratio <= 1.0:
        a = 1.5 + (N_ratio - 0.7) * 0.5 / 0.3

    ax.clear()
    ax.plot(M_Rdy_positive_list, N_Rd_positive_list, label = 'Interaction envelope - major axis - positive moment', color = '#006D62')
    ax.plot(M_Rdy_negative_list, N_Rd_negative_list, label = 'Interaction envelope - major axis - negative moment', linestyle = '--', color = '#802628')
    ax.scatter(My_02,N_Ed, label = 'First Order design action effects', color = '#88BBC2')
    ax.scatter(M_Edy,N_Ed, label = 'Second Order design action effects', color = '#C2B658')
    ax.axhline(0, color='black', linewidth=0.5)  # Horizontal line at y = 0
    ax.axvline(0, color='black', linewidth=0.5)  # Vertical line at x = 0
    ax.grid(True, which='major', linestyle='-', linewidth=0.4, color='gray', alpha=0.5)
    ax.spines['left'].set_position(('data', 0)) # Move the left spine to x=0
    ax.spines['right'].set_color('none') # Hide the right spine
    ax.spines['bottom'].set_position(('data', 0)) # Move the bottom spine to y=0
    ax.spines['top'].set_color('none') # Hide the top spine
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', direction='inout', length=6)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('M [kNm]', loc='right')  # Move the x-axis label to the right
    ax.set_ylabel('N [kN]', loc='top')
    canvas.draw()

def plot_minor_axis_failure_envelope(ax, canvas, column, N_Ed, M_Edz, Mz_02):

    Npl_Rd = (column.section.A * column.concrete_properties.f_cd + len(column.reinforcement.arrangement) * math.pi * column.reinforcement.bar_diameter**2 / 4 * column.reinforcement.f_yd)*1e-3
    # PLOT MAJOR AXIS INTERACTION DIAGRAM
    N_Rd_positive_list = []
    N_Rd_negative_list = []
    M_Rdz_positive_list = []
    M_Rdz_negative_list = []
    if column.section.shape == "rectangular" or column.section.shape == "arbitrary":
        x_limit = column.section.b
    elif column.section.shape == "circular":
        x_limit = column.section.diameter
    for x in range(1, int(x_limit*5) ,1):
        N_Rd, M_Rdz, steel_stresses, steel_strains = determine_envelope_value_minor_axis_positive(column, 0.8, x)
        N_Rd_positive_list.append(min(Npl_Rd,N_Rd))
        M_Rdz_positive_list.append(M_Rdz)

    for x in range(1, int(x_limit*5) ,1):
        # add axial force and major axis moment to list for plotting
        N_Rd, M_Rdz, steel_stresses, steel_strains = determine_envelope_value_minor_axis_negative(column, 0.8, x)
        N_Rd_negative_list.append(min(Npl_Rd,N_Rd))
        M_Rdz_negative_list.append(M_Rdz)

    ax.clear()
    ax.plot(M_Rdz_positive_list, N_Rd_positive_list, label = 'Interaction envelope - major axis - positive moment', color = '#006D62')
    ax.plot(M_Rdz_negative_list, N_Rd_negative_list, label = 'Interaction envelope - major axis - negative moment', linestyle = '--', color = '#802628')
    ax.scatter(Mz_02,N_Ed, label = 'First Order design action effects', color = '#88BBC2')
    ax.scatter(M_Edz,N_Ed, label = 'Second Order design action effects', color = '#C2B658')
    ax.axhline(0, color='black', linewidth=0.5)  # Horizontal line at y = 0
    ax.axvline(0, color='black', linewidth=0.5)  # Vertical line at x = 0
    ax.grid(True, which='major', linestyle='-', linewidth=0.4, color='gray', alpha=0.5)
    ax.spines['left'].set_position(('data', 0)) # Move the left spine to x=0
    ax.spines['right'].set_color('none') # Hide the right spine
    ax.spines['bottom'].set_position(('data', 0)) # Move the bottom spine to y=0
    ax.spines['top'].set_color('none') # Hide the top spine
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', direction='inout', length=6)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('M [kNm]', loc='right')  # Move the x-axis label to the right
    ax.set_ylabel('N [kN]', loc='top')
    canvas.draw()

def plot3Dpoints(shape,  help=False):
    if help:
        print('''
            Plots 3D interaction diagram for a given column 
            Inputs: section - section object from LSLtools
                    
              Outputs: 3D plot of interaction diagram
                  ''')
    if shape == "circular":
        # due to symmetry, only need to plot one half of the interaction diagram and rotate about centroid
        centroid = [0,0,0]
        # Number of symmetrical rotations
        n_rotations = 20
        # retrieve base [M,N] points for circular column from functions in LSL tools
        N_rds_maj_pos, M_rds_maj_pos, = determine_envelope_value_major_axis_positive()
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