from NM_interaction.utils import get_concrete_properties
from NM_interaction.plotting import plot_rectangular_section, plot_circular_section, plot_arbitrary_section, plot_major_axis_failure_envelope, plot_minor_axis_failure_envelope
from NM_interaction.resistanceCalcs import determine_envelope_value_major_axis_negative, determine_envelope_value_major_axis_positive, determine_envelope_value_minor_axis_negative, determine_envelope_value_minor_axis_positive, compute_plastic_axial_capacity
from NM_interaction.classes import Column, Concrete_Section, Reinforcement, Concrete_Material
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from NM_interaction.geometryUtils import compute_area_and_centroid, clip_polygon_at_y, clip_polygon_at_x, integrate_part_of_circle
import numpy as np

concrete_properties = get_concrete_properties("C30/37")
concrete_section = Concrete_Section(shape="arbitrary", vertices = [[-200,0], [200,0], [500,300], [250,550], [0,300], [-250,500], [-500,300], [-200,0]])
reinforcement = Reinforcement(f_yk = 500, E_s=200, bar_diameter = 16, link_diameter=10, shape="arbitrary", bar_list = [[150,50], [300,200], [450,300], [250,500], [50,300], [-50,300], [-250,500], [-450,300], [-150,50]], gamma_s=1.15)
column = Column(concrete_section, reinforcement, concrete_properties, L_eff_y=4, L_eff_z=4)

centroid =compute_area_and_centroid(column.concrete_section.polygon)[1]

print(f"centroid(y): {centroid[1]}")

fig, ax = plt.subplots()

T_pl_Rd, Npl_Rd = compute_plastic_axial_capacity(column)
    #Collect user inputs for processing
N_Rd_positive_list = []
N_Rd_negative_list = []
M_Rdy_positive_list = []
M_Rdy_negative_list = []
N_concrete = []
M_concrete = []
N_steel = []
M_steel = []
M_axial = []

for y in range(1, 50000, 5):
    #add axial force and major axis moment to list for plotting
    N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_positive(column, 0.8, y)
    N_Rd_positive_list.append(min(Npl_Rd,N_Rd))
    M_Rdy_positive_list.append(M_Rdy)
    N_Rd, M_Rdy, steel_stresses, steel_strains = determine_envelope_value_major_axis_negative(column, 0.8, y)
    N_Rd_negative_list.append(min(Npl_Rd,N_Rd))
    M_Rdy_negative_list.append(M_Rdy)

ax.clear()
ax.plot(M_Rdy_positive_list, N_Rd_positive_list, marker = "x", color = 'black')
ax.plot(M_Rdy_negative_list, N_Rd_negative_list, linestyle = '--', marker = "x", color = 'black')
ax.axhline(0, color='black', linewidth=0.5)  # Horizontal line at y = 0
ax.axvline(0, color='black', linewidth=0.5)  # Vertical line at x = 0
# Add faint gridlines to the plot
ax.grid(visible=True, which='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
ax.spines['left'].set_position(('data', 0)) # Move the left spine to x=0
ax.spines['right'].set_color('none') # Hide the right spine
ax.spines['bottom'].set_position(('data', 0)) # Move the bottom spine to y=0
ax.spines['top'].set_color('none') # Hide the top spine
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.tick_params(axis='y', direction='inout', length=6)
ax.tick_params(axis='x', rotation=90)
ax.set_title("Major Axis interaction Envelope")
ax.set_xlabel('M [kNm]', loc='right')  # Move the x-axis label to the right
ax.set_ylabel('N [kN]', loc='top')

plt.show()

