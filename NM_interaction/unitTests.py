from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from .classes import Column, Concrete_Section, Reinforcement, Concrete_Material
from .resistanceCalcs import determine_envelope_value_major_axis_positive, compute_plastic_axial_capacity
from .utils import get_concrete_properties
import numpy as np
from .geometryUtils import integrate_part_of_circle

def test_positive_rectangular_major_plot():
    print("Testing positive interaction diagram for major axis rectangular section")
    print("400x400 rectangular section, 12 bars")
    concrete_properties = get_concrete_properties("C30/37")
    h=400
    b=400
    concrete_section = Concrete_Section(shape='rectangular', h=h, b=b)
    reinforcement = Reinforcement(f_yk=500, E_s=200, bar_diameter=16, link_diameter=10, cover = 30, shape="rectangular", b = 400, h = 400, num_of_rows_of_rebar = 4, num_of_cols_of_rebar = 4, gamma_s=1.15)
    column = Column(concrete_section, reinforcement, concrete_properties, L_eff_y=3, L_eff_z=3)
    T_pl_Rd, N_pl_Rd = compute_plastic_axial_capacity(column)
    print(f"Plastic axial capacity (Npl_Rd): {N_pl_Rd} ") 
    fig = plt.figure()
    x=[]
    N = []
    M = []
    N.append(T_pl_Rd)
    M.append(0)
    for NA in range(1,4*h,5):
        N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, NA)
        x.append(NA)
        N.append(N_Rd)
        M.append(M_Rdy)
    N.append(N_pl_Rd)
    M.append(0)
    # Plot the points
    plt.scatter(M, N, color='red', marker='x')

    # Add faint gridlines to the plot
    plt.grid(visible=True, which='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.7)

    # Add solid black lines for the x=0 and y=0 axes
    plt.axhline(0, color='black', linewidth=1.5, linestyle='-')
    plt.axvline(0, color='black', linewidth=1.5, linestyle='-')

    # Add labels and legend
    plt.xlabel('M')
    plt.ylabel('N')
    plt.title('N, M Points')
    plt.show()

def test_positive_circular_major_plot():
    diameter = 400
    cover = 30
    radial_num_bars = 12
    concrete_properties = get_concrete_properties("C30/37")
    concrete_section = Concrete_Section(shape="circular",  diameter = diameter)
    reinforcement = Reinforcement(f_yk=500, E_s=200, bar_diameter=16, link_diameter=10, cover = 30, shape="circular", diameter = diameter, radial_number = radial_num_bars, gamma_s=1.15)
    column = Column(concrete_section, reinforcement, concrete_properties, L_eff_y=3, L_eff_z=3)
    T_pl_Rd, N_pl_Rd = compute_plastic_axial_capacity(column)
    x=[]
    N = []
    M = []
    N.append(T_pl_Rd)
    M.append(0)
    for NA in range(1,10*diameter):
        steel_contribution_N, steel_contribution_M, concrete_contribution_N, concrete_contribution_M,axial_contribution_M, N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, NA)
        x.append(NA)
        N.append(N_Rd)
        M.append(M_Rdy)
    N.append(N_pl_Rd)
    M.append(0)

    for i in range(len(x)):
        plt.annotate(f'x={x[i]}', (M[i], N[i]), textcoords="offset points", xytext=(5, 5), ha='center', fontsize=8)
    # Plot the points
    plt.scatter(M, N, color='red', marker='x')

    # Add faint gridlines to the plot
    plt.grid(visible=True, which='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.7)

    # Add solid black lines for the x=0 and y=0 axes
    plt.axhline(0, color='black', linewidth=1.5, linestyle='-')
    plt.axvline(0, color='black', linewidth=1.5, linestyle='-')

    # Add labels and legend
    plt.xlabel('M')
    plt.ylabel('N')
    plt.title('N, M Points')
    plt.show()


def plot_M_contributions():
    diameter = 400
    cover = 30
    radial_num_bars = 12
    concrete_properties = get_concrete_properties("C30/37")
    concrete_section = Concrete_Section(shape="circular",  diameter = diameter)
    reinforcement = Reinforcement(f_yk=500, E_s=200, bar_diameter=16, link_diameter=10, cover = 30, shape="circular", diameter = diameter, radial_number = radial_num_bars, gamma_s=1.15)
    column = Column(concrete_section, reinforcement, concrete_properties, L_eff_y=3, L_eff_z=3)
    T_pl_Rd, N_pl_Rd = compute_plastic_axial_capacity(column)
    fig = plt.figure()
    x=[]
    steel_M = []
    concrete_M = []
    axial_M = []
    for NA in range(1,800,10):
        print(f"neutral axis: {NA}")
        steel_contribution_N, steel_contribution_M, concrete_contribution_N, concrete_contribution_M,axial_contribution_M, N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, NA)
        x.append(NA)
        steel_M.append(steel_contribution_M)
        #segment_area, centroid_y, lever_arm_concrete_y = integrate_part_of_circle(NA, diameter, 0.8)
        #print(f"x: {NA}, segment area: {segment_area}, centroid_y: {centroid_y}, lever arm: {lever_arm_concrete_y}")
        concrete_M.append(concrete_contribution_M)
        axial_M.append(axial_contribution_M)
    print("Plotting Moment resistance contributions from concrete, steel and out of balance axial force")

    # Plot x vs M with contributions from different components
    fig, ax = plt.subplots()

    # Bar width and offsets
    bar_width = 0.2
    x_positions = np.arange(len(x))  # Create positions for the x bins

    # Plot bars next to each other
    plt.bar(x_positions - bar_width, steel_M, bar_width, label='Steel Moments', color='red')
    plt.bar(x_positions, concrete_M, bar_width, label='Concrete Moments', color='blue')
    plt.bar(x_positions + bar_width, axial_M, bar_width, label='Axial Moments', color='green')

    # Add labels and legend
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Add a horizontal line at y=0
    plt.xticks(x_positions, x)  # Set x-ticks to the original x values
    plt.xlabel('x')
    plt.ylabel('M')
    plt.title('x vs M with Component Contributions (Side-by-Side Bars)')
    plt.legend()

    plt.show()

def plot_N_contributions():
    diameter = 400
    cover = 30
    radial_num_bars = 12
    concrete_properties = get_concrete_properties("C30/37")
    concrete_section = Concrete_Section(shape="circular",  diameter = diameter)
    reinforcement = Reinforcement(f_yk=500, E_s=200, bar_diameter=16, link_diameter=10, cover = 30, shape="circular", diameter = diameter, radial_number = radial_num_bars, gamma_s=1.15)
    column = Column(concrete_section, reinforcement, concrete_properties, L_eff_y=3, L_eff_z=3)
    T_pl_Rd, N_pl_Rd = compute_plastic_axial_capacity(column)
    fig = plt.figure()
    x=[]
    N = []
    steel_N = []
    concrete_N = []
    for NA in range(1,800,10):
        steel_contribution_N, steel_contribution_M, concrete_contribution_N, concrete_contribution_M,axial_contribution_M, N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, NA)
        x.append(NA)
        steel_N.append(steel_contribution_N)
        concrete_N.append(concrete_contribution_N)
        N.append(steel_contribution_N+concrete_contribution_N)
    print("Plotting Moment resistance contributions from concrete, steel and out of balance axial force")

    # Plot x vs M with contributions from different components
    fig, ax = plt.subplots()

    # Bar width and offsets
    bar_width = 0.2
    x_positions = np.arange(len(x))  # Create positions for the x bins

    # Plot bars next to each other
    plt.bar(x_positions - bar_width, steel_N, bar_width, label='Sum of Steel axial forces', color='red')
    plt.bar(x_positions, concrete_N, bar_width, label='Concrete axial forces', color='blue')
    plt.bar(x_positions + bar_width, N, bar_width, label='Resultant axial', color='green')

    # Add labels and legend
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Add a horizontal line at y=0
    plt.xticks(x_positions, x)  # Set x-ticks to the original x values
    plt.xlabel('x')
    plt.ylabel('M')
    plt.title('x vs M with Component Contributions (Side-by-Side Bars)')
    plt.legend()

    plt.show()

def review_steel_strains():
    concrete_properties = get_concrete_properties("C30/37")
    h=400
    b=400
    concrete_section = Concrete_Section(shape='rectangular', h=h, b=b)
    reinforcement = Reinforcement(f_yk=500, E_s=200, bar_diameter=16, link_diameter=10, cover = 30, shape="rectangular", b = 400, h = 400, num_of_rows_of_rebar = 4, num_of_cols_of_rebar = 4, gamma_s=1.15)
    column = Column(concrete_section, reinforcement, concrete_properties, L_eff_y=3, L_eff_z=3)
    fig = plt.figure()
    x=[]
    strains = []
    bar_coords = column.reinforcement.arrangement
    concrete_N = []
    print(f"rebar coords: {bar_coords}")
    for NA in range(1,40):
        print(f"neutral axis: {NA}")
        steel_contribution_N, steel_contribution_M, concrete_contribution_N, concrete_contribution_M,axial_contribution_M, N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, NA)
        print(f"Rebar Strain: {steel_strains}")
        print(f"Rebar Stress: {steel_stresses}")
        print(f"Sum of steel forces: {0.001*sum(steel_stresses)*3.14*16**2/4} kN")
        x.append(NA)
        strains.append(steel_strains)




