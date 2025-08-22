import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize 
from shapely.geometry import Polygon, LineString, GeometryCollection, Point, box
from shapely.ops import split
import pandas as pd
import math
from scipy.integrate import quad
from matplotlib import font_manager
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Circle, Rectangle
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
#import forallpeople

#forallpeople.environment('structural', top_level=True)

class Section():
    #Section class
    def __init__(self, shape= None, h=None, b=None, diameter = None, vertices = None):
        self.shape = shape
        if shape == 'rectangular':
            self.b = float(b)#*mm
            self.h = float(h)#*mm
            A = b*h
            self.A = A
            I_yy = b*h**3/12
            I_zz = h*b**3/12
            self.I_yy = I_yy
            self.I_zz = I_zz
            self.i_y = (I_yy / A)**0.5
            self.i_z = (I_zz / A)**0.5
        elif shape == 'circular':
            self.diameter = float(diameter)#*mm
            A = math.pi*diameter**2/4
            self.A = A
            I_yy = math.pi*diameter**4/64
            self.I_yy = I_yy 
            I_zz = math.pi*diameter**4/64
            self.I_zz = I_zz
            self.i_y = (I_yy / A)**0.5
            self.i_z = (I_zz / A)**0.5
        elif shape == 'arbitrary':
            self.vertices = vertices
            self.polygon = Polygon(vertices)
            self.left_of_section = min(self.polygon.exterior.coords.xy[0])#*mm
            self.right_of_section = max(self.polygon.exterior.coords.xy[0])#*mm
            self.top_of_section = max(self.polygon.exterior.coords.xy[1])#*mm
            self.bottom_of_section = min(self.polygon.exterior.coords.xy[1])#*mm
            self.h = (self.top_of_section - self.bottom_of_section)
            self.b = (self.right_of_section - self.left_of_section)
            self.A = self.polygon.area#*mm**2
            self.I_zz = abs(compute_second_moment_area(self.polygon) [0]) #mm**4
            self.I_yy = abs(compute_second_moment_area(self.polygon) [1]) #mm**4
            self.i_y = (self.I_yy / self.A)**0.5
            self.i_z = (self.I_zz / self.A)**0.5


class Column():
    #Column class
    def __init__(self, section, reinforcement, concrete_properties, L_eff_y, L_eff_z):
        self.section = section
        self.reinforcement = reinforcement
        self.concrete_properties = concrete_properties
        self.L_eff_y = float(L_eff_y)#*m
        self.L_eff_z = float(L_eff_z)#*m
        if section.shape == 'rectangular':
            self.dy = (section.h - reinforcement.cover - reinforcement.link_diameter - reinforcement.bar_diameter/2)#*mm
            self.d_2 = (reinforcement.cover + reinforcement.link_diameter + reinforcement.bar_diameter/2)#*mm
            self.dz = (section.b - reinforcement.cover - reinforcement.link_diameter - reinforcement.bar_diameter/2)#*mm
        elif section.shape == 'circular':
            self.d = (section.diameter - reinforcement.cover - reinforcement.link_diameter - reinforcement.bar_diameter/2)#*mm
            self.d_2 = (reinforcement.cover + reinforcement.link_diameter + reinforcement.bar_diameter/2)#*mm
        elif section.shape == 'arbitrary':
            self.dy = self.section.h - min(coord[1] for coord in self.reinforcement.arrangement)
            self.dy_2 = min(coord[1] for coord in self.reinforcement.arrangement)
            self.dz = self.section.b - min(coord[0] for coord in self.reinforcement.arrangement)
            self.dz_2 = min(coord[0] for coord in self.reinforcement.arrangement)
            # reinforcement for arbitrary sections are defined by the user
            

class Concrete_Material():
    def __init__(self, concrete_grade, f_ck, f_ck_cube, f_cm, f_ctm, f_ctk_0_05, f_ctk_0_95, E_cm, eps_c1, eps_c3, eps_cu1, eps_cu2, eps_cu3, eta, gamma_c = 1.5, alpha_cc = 0.85):
        self.grade = str(concrete_grade)
        self.f_ck = float(f_ck)#*MPa
        self.f_ck_cube = float(f_ck_cube)#*MPa
        self.f_cd = float(alpha_cc*f_ck/gamma_c)#*MPa
        self.f_cm = float(f_cm)#*MPa
        self.f_ctm = float(f_ctm)#*MPa
        self.f_ctk_0_05 = float(f_ctk_0_05)#*MPa
        
        self.f_ctk_0_95 = float(f_ctk_0_95)#*MPa
        self.E_cm = float(E_cm)#*GPa
        self.eps_c1 = float(eps_c1)
        self.eps_c3 = float(eps_c3)
        self.eps_cu1 = float(eps_cu1)
        self.eps_cu2 = float(eps_cu2)
        self.eps_cu3 = float(eps_cu3)
        self.eta = float(eta)

class Reinforcement():
    def __init__(self, f_yk, E_s, bar_diameter, link_diameter, cover = None, shape=None, b =None, h = None, diameter = None, num_of_rows_of_rebar = None, num_of_cols_of_rebar = None, radial_number = None, rebar_coords = None, gamma_s=1.15):
        self.shape = shape
        self.f_yk = float(f_yk)#*MPa
        self.E_s = float(E_s)#*GPa
        self.bar_diameter = float(bar_diameter)#*mm
        self.link_diameter = float(link_diameter)#*mm
        self.f_yd = float(f_yk / gamma_s)#*MPa
        self.eps_ud = float((f_yk / gamma_s) /E_s)
        self.cover = cover#*mm
        A_s = 0
        if self.shape == 'rectangular':
            rebar_coords = []
            width_x = b - 2*cover - 2*link_diameter - bar_diameter
            height_y = h - 2*cover - 2*link_diameter - bar_diameter
            spacing_x = width_x / (num_of_cols_of_rebar-1)
            spacing_y = height_y / (num_of_rows_of_rebar - 1)
            #create rebar in top + bottom layer
            for i in range(num_of_cols_of_rebar):
                x1 = (cover + link_diameter + bar_diameter / 2 + i*spacing_x)#*mm
                y1 = (cover + link_diameter + bar_diameter/2)#*mm
                rebar_coords.append([x1, y1])
                x2 = x1
                y2 = h - cover - link_diameter - bar_diameter/2
                rebar_coords.append([x2,y2])
                A_s += 2* math.pi * (bar_diameter/2)**2
                print(f'A_s: {A_s}')
            # create side bars
            for j in range(num_of_rows_of_rebar-2):
                x1 = cover + link_diameter + bar_diameter/2
                y1 = cover+link_diameter + bar_diameter/2 + spacing_y*(j+1)
                rebar_coords.append([x1, y1])
                x2 = b - cover - link_diameter - bar_diameter/2
                y2 = y1
                rebar_coords.append([x2,y2])
                A_s += 2* math.pi * (bar_diameter/2)**2
                print(f'A_s: {A_s}')

        if self.shape == 'circular':
            rebar_coords = []
            r = diameter / 2 - cover - link_diameter - bar_diameter/2
            radial_spacing = 2*math.pi / radial_number
            for i in range(radial_number):
                x = r*math.cos(radial_spacing*i) + diameter/2
                y = r*math.sin(radial_spacing*i) + diameter/2
                rebar_coords.append([x,y])
                A_s += math.pi * (bar_diameter/2)**2

        if self.shape == 'arbitrary':
            A_s = len(rebar_coords)*math.pi * (bar_diameter/2)**2

        self.arrangement = rebar_coords
        self.A_s = A_s

class InputApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RC Column Designer")
        
        # Radio Button Variable
        self.shape_var = tk.StringVar(value="rectangular")
        
        self.coord_list = []
        self.bar_list = []

        # Create the main layout
        self.create_main_layout()

    def create_main_layout(self):
        # Options input for concrete grade
        options = ["C12/15","C16/20","C20/25","C25/30","C30/37","C35/45","C40/50","C45/55","C50/60","C55/67","C60/75","C70/85","C80/95","C90/105"]
        #Create row for concrete choice
        concrete_choice_label = ttk.Label(self.root, text="Concrete Grade:")
        concrete_choice_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        
        # Create the Combobox
        self.concrete_dropdown = ttk.Combobox(root, values=options)
        self.concrete_dropdown.current()  # Set the default option (index 0)
        self.concrete_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        self.concrete_dropdown.bind("<<ComboboxSelected>>", self.dropdown_selected)

        f_yk_label = ttk.Label(self.root, text="Reinforcement yield strength, f_yk:")
        f_yk_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")
        self.f_yk_entry = ttk.Entry(self.root)
        self.f_yk_entry.grid(row=1, column=1, padx=10, pady=5, sticky="ew")
        self.f_yk_entry.insert(0,"500")
        ttk.Label(self.root, text="MPa").grid(row=1, column=2, padx=10, pady=5, sticky="w")

        E_s_label = ttk.Label(self.root, text="Reinforcing steel Young's Modulus, E_s:")
        E_s_label.grid(row=2, column=0, padx=10, pady=5, sticky="e")
        self.E_s_entry = ttk.Entry(self.root)
        self.E_s_entry.grid(row=2, column=1, padx=10, pady=5, sticky="ew")
        self.E_s_entry.insert(0,"210")
        ttk.Label(self.root, text="GPa").grid(row=2, column=2, padx=10, pady=5, sticky="w")

        bar_dia_label = ttk.Label(self.root, text="Main bar Diameter:")
        bar_dia_label.grid(row=3, column=0, padx=10, pady=5, sticky="e")
        self.bar_dia_entry = ttk.Entry(self.root)
        self.bar_dia_entry.grid(row=3, column=1, padx=10, pady=5, sticky="ew")
        self.bar_dia_entry.insert(0,"25")
        ttk.Label(self.root, text="mm").grid(row=3, column=2, padx=10, pady=5, sticky="w")

        link_dia_label = ttk.Label(self.root, text="Link Diameter:")
        link_dia_label.grid(row=4, column=0, padx=10, pady=5, sticky="e")
        self.link_dia_entry = ttk.Entry(self.root)
        self.link_dia_entry.grid(row=4, column=1, padx=10, pady=5, sticky="ew")
        self.link_dia_entry.insert(0, "10")
        ttk.Label(self.root, text="mm").grid(row=4, column=2, padx=10, pady=5, sticky="w")

        # Section shape label
        shape_label = ttk.Label(self.root, text="Section shape:")
        shape_label.grid(row=5, column=0, padx=10, pady=10, sticky="w")

        # --- Add Right-Hand Plots below the buttons ---
        # Create a frame for the two graphs below the buttons
        
        self.lengthloads_frame = ttk.Frame(self.root)
        self.lengthloads_frame.grid(row=0, column=6, rowspan=5, columnspan=5, padx=0, pady=0, sticky="nsew")

        self.L_Effy_label = ttk.Label(self.lengthloads_frame, text="L_eff,y:")
        self.L_Effy_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        self.L_Effy_entry = ttk.Entry(self.lengthloads_frame)
        self.L_Effy_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")
        self.L_Effy_units = ttk.Label(self.lengthloads_frame, text="m")
        self.L_Effy_units.grid(row=0, column=2, padx=10, pady=5, sticky="w")

        self.L_effz_label = ttk.Label(self.lengthloads_frame, text="L_eff,z:")
        self.L_effz_label.grid(row=0, column=3, padx=10, pady=5, sticky="e")
        self.L_effz_entry = ttk.Entry(self.lengthloads_frame)
        self.L_effz_entry.grid(row=0, column=4, padx=10, pady=5, sticky="ew")
        self.L_effz_units = ttk.Label(self.lengthloads_frame, text="m")
        self.L_effz_units.grid(row=0, column=5, padx=10, pady=5, sticky="w")

        self.N_Ed_label = ttk.Label(self.lengthloads_frame, text="N_Ed:")
        self.N_Ed_label.grid(row=2, column=0, padx=10, pady=5, sticky="e")
        self.N_Ed_entry = ttk.Entry(self.lengthloads_frame)
        self.N_Ed_entry.grid(row=2, column=1, padx=10, pady=5, sticky="ew")
        self.N_Ed_units = ttk.Label(self.lengthloads_frame, text="kN")
        self.N_Ed_units.grid(row=2, column=2, padx=10, pady=5, sticky="w")

        self.First_order_moments_label = ttk.Label(self.lengthloads_frame, text="First Order Moments:")
        self.First_order_moments_label.grid(row=3, column=0, padx=10, pady=5, sticky="w")

        self.M_y_top_label = ttk.Label(self.lengthloads_frame, text="M_y [top]:")
        self.M_y_top_label.grid(row=4, column=0, padx=10, pady=5, sticky="e")
        self.M_y_top_entry = ttk.Entry(self.lengthloads_frame)
        self.M_y_top_entry.grid(row=4, column=1, padx=10, pady=5, sticky="ew")
        self.M_y_top_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_y_top_units.grid(row=4, column=2, padx=10, pady=5, sticky="w")

        self.M_y_bottom_label = ttk.Label(self.lengthloads_frame, text="M_y [bottom]:")
        self.M_y_bottom_label.grid(row=5, column=0, padx=10, pady=5, sticky="e")
        self.M_y_bottom_entry = ttk.Entry(self.lengthloads_frame)
        self.M_y_bottom_entry.grid(row=5, column=1, padx=10, pady=5, sticky="ew")
        self.M_y_bottom_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_y_bottom_units.grid(row=5, column=2, padx=10, pady=5, sticky="w")

        self.M_z_top_label = ttk.Label(self.lengthloads_frame, text="M_z [top]:")
        self.M_z_top_label.grid(row=4, column=3, padx=10, pady=5, sticky="e")
        self.M_z_top_entry = ttk.Entry(self.lengthloads_frame)
        self.M_z_top_entry.grid(row=4, column=4, padx=10, pady=5, sticky="ew")
        self.M_z_top_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_z_top_units.grid(row=4, column=5, padx=10, pady=5, sticky="w")

        self.M_z_bottom_label = ttk.Label(self.lengthloads_frame, text="M_z [botom]:")
        self.M_z_bottom_label.grid(row=5, column=3, padx=10, pady=5, sticky="e")
        self.M_z_bottom_entry = ttk.Entry(self.lengthloads_frame)
        self.M_z_bottom_entry.grid(row=5, column=4, padx=10, pady=5, sticky="ew")
        self.M_z_bottom_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_z_bottom_units.grid(row=5, column=5, padx=10, pady=5, sticky="w")

        self.plot_major_axis_button = ttk.Button(self.root, text="Plot Major Axis Failure Envelope", command = self.plot_major_axis_failure_envelope)
        self.plot_major_axis_button.grid(row=5, column=6, padx=10, pady=5, sticky="ew")

        self.plot_minor_axis_button = ttk.Button(self.root, text="Plot Minor Axis Failure Envelope", command = self.plot_minor_axis_failure_envelope)
        self.plot_minor_axis_button.grid(row=5, column=8, padx=10, pady=5, sticky="ew")

        self.plot_frame = ttk.Frame(self.root)
        self.plot_frame.grid(row=6, column=6, columnspan=4, rowspan=10, padx=10, pady=10, sticky="nsew")

        # Add first plot frame
        self.plot1_frame = ttk.Frame(self.plot_frame)
        self.plot1_frame.grid(row=0, column=0, padx=5, pady=5)
        self.fig1, self.ax1 = plt.subplots(figsize=(4, 4))
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.plot1_frame)
        self.canvas1.get_tk_widget().pack()
        self.ax1.set_title("Major Axis Plot")
        self.ax1.set_xlabel("M [kNm]")
        self.ax1.set_ylabel("N [kN]")

        # Add second plot frame next to the first
        self.plot2_frame = ttk.Frame(self.plot_frame)
        self.plot2_frame.grid(row=0, column=2, padx=5, pady=5)
        self.fig2, self.ax2 = plt.subplots(figsize=(4, 4))
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.plot2_frame)
        self.canvas2.get_tk_widget().pack()
        self.ax2.set_title("Minor Axis Plot")
        self.ax2.set_xlabel("M [kNm]")
        self.ax2.set_ylabel("N [kN]")

        # Radio buttons for shape selection
        self.radio_frame = ttk.Frame(self.root)
        self.radio_frame.grid(row=5, column=1, padx=10, pady=10)
        
        shapes = ["rectangular", "circular", "arbitrary"]
        for shape in shapes:
            ttk.Radiobutton(self.radio_frame, text=shape.capitalize(), value=shape, variable=self.shape_var, command=self.update_input_fields).pack(side="left", padx=5)

        # Create frame for dynamic input fields
        self.input_frame = ttk.Frame(self.root)
        self.input_frame.grid(row=6, column=0, columnspan = 3, padx=5,  pady=10)

        # Initialize input fields for the default selection (rectangular)
        self.update_input_fields()

    def dropdown_selected(self, event):
        print(f"Selected option: {self.concrete_dropdown.get()}")

    def update_input_fields(self):
        # Clear any existing widgets in the input_frame
        for widget in self.input_frame.winfo_children():
            widget.destroy()

        # Display input fields based on the selected shape
        shape = self.shape_var.get()
        
        if shape == "rectangular":
            self.create_rectangular_inputs()
        elif shape == "circular":
            self.create_circular_inputs()
        elif shape == "arbitrary":
            self.create_arbitrary_inputs()

    def create_rectangular_inputs(self):
        # Rectangular input fields
        ttk.Label(self.input_frame, text="Section width, b:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.b_input = ttk.Entry(self.input_frame)
        self.b_input.insert(0, "500")
        self.b_input.grid(row=0, column=2, padx=5, pady=5)
        self.b_input.bind("<Return>", self.plot_rectangular_section)
        self.b_input.bind("<FocusOut>", self.plot_rectangular_section)
        ttk.Label(self.input_frame, text="mm").grid(row=0, column=3, padx=5, pady=5, sticky="w")
        
        ttk.Label(self.input_frame, text="Section height, h:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.h_input = ttk.Entry(self.input_frame)
        self.h_input.insert(0, "500")
        self.h_input.grid(row=1, column=2, padx=5, pady=5)
        self.h_input.bind("<Return>", self.plot_rectangular_section)
        self.h_input.bind("<FocusOut>", self.plot_rectangular_section)
        ttk.Label(self.input_frame, text="mm").grid(row=1, column=3, padx=5, pady=5, sticky="w")
        
        ttk.Label(self.input_frame, text="Cover:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.cover_input = ttk.Entry(self.input_frame)
        self.cover_input.insert(0, "25")
        self.cover_input.grid(row=2, column=2, padx=5, pady=5)
        self.cover_input.bind("<Return>",self.plot_rectangular_section)
        self.cover_input.bind("<FocusOut>",self.plot_rectangular_section)
        ttk.Label(self.input_frame, text="mm").grid(row=2, column=3, padx=5, pady=5, sticky="w")

        ttk.Label(self.input_frame, text="Number of bars distributed along the width of the column\nn_x:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.n_x_input = ttk.Entry(self.input_frame)
        self.n_x_input.grid(row=3, column=2, padx=5, pady=5)
        self.n_x_input.bind("<Return>",self.plot_rectangular_section)
        self.n_x_input.bind("<FocusOut>",self.plot_rectangular_section)
        ttk.Label(self.input_frame, text="#").grid(row=3, column=3, padx=5, pady=5, sticky="w")

        ttk.Label(self.input_frame, text="Number of bars distributed along the height of the column\n n_y:").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.n_y_input = ttk.Entry(self.input_frame)
        self.n_y_input.grid(row=4, column=2, padx=5, pady=5)
        self.n_y_input.bind("<Return>",self.plot_rectangular_section)
        self.n_y_input.bind("<FocusOut>",self.plot_rectangular_section)
        ttk.Label(self.input_frame, text="#").grid(row=4, column=3, padx=5, pady=5, sticky="w")

        self.fig0, self.ax0 = plt.subplots(figsize=(3,3))
        self.canvas0 = FigureCanvasTkAgg(self.fig0, master=self.input_frame)
        self.canvas0.get_tk_widget().grid(row=5, column=0, columnspan=4, padx=5, pady=5)
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)
        #self.ax0.set_axis_on()
        self.ax0.grid(True)
        self.canvas0.draw()

    def create_circular_inputs(self):
        # Circular input fields
        ttk.Label(self.input_frame, text="Diameter, d:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.diameter_entry = ttk.Entry(self.input_frame)
        self.diameter_entry.grid(row=0, column=1, padx=5, pady=5)
        self.diameter_entry.bind("<Return>", self.plot_circular_section)  # Update plot on pressing Enter
        self.diameter_entry.bind("<FocusOut>", self.plot_circular_section)  # Update plot on leaving the field
        
        ttk.Label(self.input_frame, text="Number of bars:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.radial_num_bars_input = ttk.Entry(self.input_frame)
        self.radial_num_bars_input.grid(row=1, column=1, padx=5, pady=5)
        self.radial_num_bars_input.bind("<Return>", self.plot_circular_section)
        self.radial_num_bars_input.bind("<FocusOut>",self.plot_circular_section)
        
        ttk.Label(self.input_frame, text="Cover:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.cover_input = ttk.Entry(self.input_frame)
        self.cover_input.grid(row=2, column=1, padx=5, pady=5)
        self.cover_input.bind("<Return>", self.plot_circular_section)
        self.cover_input.bind("<FocusOut>", self.plot_circular_section)

        self.fig0, self.ax0 = plt.subplots(figsize=(3,3))
        self.canvas0 = FigureCanvasTkAgg(self.fig0, master=self.input_frame)
        self.canvas0.get_tk_widget().grid(row=3, column=0, columnspan=4, padx=5, pady=5)
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)
        #self.ax0.set_axis_on()
        self.ax0.grid(True)
        self.canvas0.draw()

    def plot_rectangular_section(self, event=None):
        self.ax0.clear()
        rect = Rectangle((0,0), float(self.b_input.get()), float(self.h_input.get()), facecolor="blue", edgecolor="black", alpha=0.5)
        self.ax0.add_patch(rect)
        self.ax0.set_xlim(0, 1.1*max(float(self.b_input.get()), float(self.h_input.get())))
        self.ax0.set_ylim(0, 1.1*max(float(self.b_input.get()), float(self.h_input.get())))
        rebar_coords = self.arrange_rectangular_bars()
        self.ax0.scatter(rebar_coords[:,0], rebar_coords[:,1], color = "red")
        # Draw the updated canvas
        self.ax0.grid(True)
        self.ax0.set_axis_on()
        self.canvas0.draw()
        
    def arrange_rectangular_bars(self):
        rebar_coords = []
        cover = float(self.cover_input.get())
        b = float(self.b_input.get())
        h = float(self.h_input.get())
        link_dia = float(self.link_dia_entry.get())
        phi = float(self.bar_dia_entry.get())
        n_x = int(self.n_x_input.get())
        n_y = int(self.n_y_input.get())
        spacing_x = float((b - 2*cover - 2*link_dia - phi) / (n_x - 1))
        spacing_y = float((h - 2*cover - 2*link_dia - phi) / (n_y - 1))
        
        #create rebar in top + bottom layer
        for i in range(n_x):
            x1 = cover + link_dia + phi/2 + i*spacing_x
            y1 = cover + link_dia + phi/2
            rebar_coords.append([x1, y1])
            x2 = x1
            y2 = h - cover - link_dia - phi/2
            rebar_coords.append([x2,y2])
        # create side bars
        for j in range(n_y-2):
            x1 = cover + link_dia + phi/2
            y1 = cover + link_dia + phi/2 + spacing_y * (j+1)
            rebar_coords.append([x1, y1])
            x2 = b - cover - link_dia - phi/2
            y2 = y1
            rebar_coords.append([x2,y2])
        return np.array(rebar_coords)

    def plot_circular_section(self, event=None):
        # Create a filled circle
        radius = float(self.diameter_entry.get())/2
        circle = Circle([0,0], radius, color='blue', alpha=0.5)
        self.ax0.clear()
        # Add the circle to the axis
        self.ax0.add_patch(circle)
        
        # Set limits so the circle is properly shown
        self.ax0.set_xlim(- 1.1*float(self.diameter_entry.get())/2 , 1.1*float(self.diameter_entry.get())/2 )
        self.ax0.set_ylim( - 1.1*float(self.diameter_entry.get())/2, 1.1*float(self.diameter_entry.get())/2 )
        
        if self.radial_num_bars_input.get() and float(self.radial_num_bars_input.get()) > 0 and self.diameter_entry.get() and float(self.diameter_entry.get()) > 0 and self.cover_input.get() and float(self.cover_input.get()) > 0:
            theta = 2*math.pi / float(self.radial_num_bars_input.get())
            bar_radius = float(self.diameter_entry.get())/2 - float(self.cover_input.get())
            x = []
            y = []
            for bar in range(int(self.radial_num_bars_input.get())): 
                x.append(bar_radius*math.sin(theta*bar))
                y.append(bar_radius*math.cos(theta*bar))
            self.ax0.scatter(x,y, color = "red")    

        # Set equal aspect ratio to ensure the circle is round
        self.ax0.set_aspect('equal', 'box')
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)

        # Draw the updated canvas
        self.ax0.grid(True)
        self.ax0.set_axis_on()
        self.canvas0.draw()

    def create_arbitrary_inputs(self):
        # Arbitrary input fields
        ttk.Label(self.input_frame, text="Input vertices for section boundary [mm]:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.x_input = ttk.Entry(self.input_frame, width=10)
        self.x_input.grid(row=0, column=1, padx=5, pady=5)
        self.y_input = ttk.Entry(self.input_frame, width=10)
        self.y_input.grid(row=0, column=2, padx=5, pady=5)
        self.add_vertex_button = ttk.Button(self.input_frame, text="Add boundary point", command=self.add_vertex)
        self.add_vertex_button.grid(row=0, column=3, padx=5, pady=5)
        self.close_boundary_button = ttk.Button(self.input_frame, text="Close section boundary", command=self.close_boundary)
        self.close_boundary_button.grid(row=0, column=4, padx=5, pady=5)

        self.coord_display = tk.Text(self.input_frame, height=5, width=40, state="disabled")
        self.coord_display.grid(row=1, column=0, rowspan = 2, columnspan=4, padx=5, pady=5)
        # Add buttons for clearing and deleting last point
        self.delete_vertex_button = ttk.Button(self.input_frame, text="Delete Last Point", command=self.delete_last_point)
        self.delete_vertex_button.grid(row=1, column=3, padx=5, pady=5)
        self.clear_vertices_button = ttk.Button(self.input_frame, text="Clear Points", command=self.clear_vertices)
        self.clear_vertices_button.grid(row=2, column=3, padx=5, pady=5)

        ttk.Label(self.input_frame, text="Input [x,y] coordinates for reinforcement [mm]:").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.bar_x_input = ttk.Entry(self.input_frame, width=10)
        self.bar_x_input.grid(row=4, column=1, padx=5, pady=5)
        self.bar_y_input = ttk.Entry(self.input_frame, width=10)
        self.bar_y_input.grid(row=4, column=2, padx=5, pady=5)
        self.add_bar_button = ttk.Button(self.input_frame, text="Add reinforcement bar", command=self.add_bar)
        self.add_bar_button.grid(row=4, column=3, padx=5, pady=5)

        self.bar_display = tk.Text(self.input_frame, height=5, width=40, state="disabled")
        self.bar_display.grid(row=5, column=0, rowspan = 2, columnspan=4, padx=5, pady=5)
        self.delete_bar_button = ttk.Button(self.input_frame, text="Delete Last bar Input", command=self.delete_last_bar)
        self.delete_bar_button.grid(row=5, column=3, padx=5, pady=5)
        self.clear_bars_button = ttk.Button(self.input_frame, text="Clear all bars", command=self.clear_bars)
        self.clear_bars_button.grid(row=6, column=3, padx=5, pady=5)

        self.fig0, self.ax0 = plt.subplots(figsize=(3,3))
        self.canvas0 = FigureCanvasTkAgg(self.fig0, master=self.input_frame)
        self.canvas0.get_tk_widget().grid(row=7, column=0, columnspan=4, padx=5, pady=5)
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)
        self.ax0.set_axis_on()
        self.ax0.grid(True)
        self.canvas0.draw()

    def add_bar(self):
        try:
            x = float(self.bar_x_input.get())
            y = float(self.bar_y_input.get())
            self.bar_list.append([x, y])

            self.update_bar_display()
            self.update_section_plot()

            self.bar_x_input.delete(0, tk.END)
            self.bar_y_input.delete(0, tk.END)
        
        except ValueError:
            # Handle invalid inputs
            print("Invalid input. Please enter numerical values.")

    def add_vertex(self):
        # Add coordinates to the list and display them
        try:
            x = float(self.x_input.get())
            y = float(self.y_input.get())
            self.coord_list.append([x, y])
            
            # Update the list of coordiantes
            self.update_coord_display()
            # Update the displayed geometry
            self.update_section_plot()
            
            # Clear the input fields
            self.x_input.delete(0, tk.END)
            self.y_input.delete(0, tk.END)
        except ValueError:
            # Handle invalid inputs
            print("Invalid input. Please enter numerical values.")

    def update_coord_display(self):
        # Enable the text widget to update it
        self.coord_display.config(state="normal")
        self.coord_display.delete(1.0, tk.END)
        
        # Display all coordinates
        for coord in self.coord_list:
            self.coord_display.insert(tk.END, f"[{coord[0]}, {coord[1]}]\n")
        
        # Disable the text widget to prevent editing
        self.coord_display.config(state="disabled")

    def close_boundary(self):
        try:
            self.coord_list.append((self.coord_list[0]))
            # Update the list of coordiantes
            self.update_coord_display()

            # Update the displayed geometry
            self.update_section_plot()
        except:
            print("unable to close boundary")

    def update_bar_display(self):
        self.bar_display.config(state="normal")
        self.bar_display.delete(1.0, tk.END)

        for coord in self.bar_list:
            self.bar_display.insert(tk.END, f"[{coord[0]}, {coord[1]}]\n")

        # Disable the text widget to prevent editing
        self.coord_display.config(state="disabled")

    def update_section_plot(self):
        # Clear the previous plot
        self.ax0.clear()

        # Unpack the list of coordinates into x and y values
        if self.coord_list:
            vertices_x_vals, vertices_y_vals = zip(*self.coord_list)
            # Draw vectors (lines) between consecutive points
            self.ax0.plot(vertices_x_vals, vertices_y_vals, linestyle='-', marker='x', color='green')
            # Close the boundary by connecting the last point to the first
            x_vals_closed = list(vertices_x_vals) + [vertices_x_vals[0]]
            y_vals_closed = list(vertices_y_vals) + [vertices_y_vals[0]]

            # Fill the polygon formed by the boundary points and shade it
            self.ax0.fill(x_vals_closed, y_vals_closed, 'lightblue', alpha=0.5)

        if self.bar_list:
            bar_x_vals, bar_y_vals = zip(*self.bar_list)
            self.ax0.scatter(bar_x_vals, bar_y_vals, marker='o', color='red')
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)
        # Draw the updated canvas
        self.canvas0.draw()

    def delete_last_point(self):
        """Delete the last entered point from the list and update the graph."""
        if self.coord_list:
            self.coord_list.pop()
            self.update_coord_display()
            self.update_section_plot()

    def clear_vertices(self):
        """Clear all points from the list and update the graph."""
        self.coord_list.clear()
        self.update_coord_display()
        self.update_section_plot()

    def delete_last_bar(self):
        """Delete the last entered point from the list and update the graph."""
        if self.coord_list:
            self.bar_list.pop()
            self.update_bar_display()
            self.update_section_plot()

    def clear_bars(self):
        """Clear all points from the list and update the graph."""
        self.bar_list.clear()
        self.update_bar_display()
        self.update_section_plot()

    def plot_major_axis_failure_envelope(self, shape = None, h = None, b = None, diameter = None, vertices = None):
        #Collect user inputs for processing
        N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column = self.collect_user_input()
        My_01, My_02, Mz_01, Mz_02, M_Edy, M_Edz, slenderness_y, slenderness_z, slenderness_ratio_y, slenderness_ratio_z = check_slenderness(column, N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom)
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

        self.ax1.clear()
        self.ax1.plot(M_Rdy_positive_list, N_Rd_positive_list, label = 'Interaction envelope - major axis - positive moment', color = '#006D62')
        self.ax1.plot(M_Rdy_negative_list, N_Rd_negative_list, label = 'Interaction envelope - major axis - negative moment', linestyle = '--', color = '#802628')
        self.ax1.scatter(My_02,N_Ed, label = 'First Order design action effects', color = '#88BBC2')
        self.ax1.scatter(M_Edy,N_Ed, label = 'Second Order design action effects', color = '#C2B658')
        self.ax1.axhline(0, color='black', linewidth=0.5)  # Horizontal line at y = 0
        self.ax1.axvline(0, color='black', linewidth=0.5)  # Vertical line at x = 0
        self.ax1.grid(True, which='major', linestyle='-', linewidth=0.4, color='gray', alpha=0.5)
        self.ax1.spines['left'].set_position(('data', 0)) # Move the left spine to x=0
        self.ax1.spines['right'].set_color('none') # Hide the right spine
        self.ax1.spines['bottom'].set_position(('data', 0)) # Move the bottom spine to y=0
        self.ax1.spines['top'].set_color('none') # Hide the top spine
        self.ax1.xaxis.set_ticks_position('bottom')
        self.ax1.yaxis.set_ticks_position('left')
        self.ax1.tick_params(axis='y', direction='inout', length=6)
        self.ax1.tick_params(axis='x', rotation=90)
        self.ax1.set_xlabel('M [kNm]', loc='right')  # Move the x-axis label to the right
        self.ax1.set_ylabel('N [kN]', loc='top')
        self.canvas1.draw()  
        
    def plot_minor_axis_failure_envelope(self):
        N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column = self.collect_user_input()
        My_01, My_02, Mz_01, Mz_02, M_Edy, M_Edz, slenderness_y, slenderness_z, slenderness_ratio_y, slenderness_ratio_z = check_slenderness(column, N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom)

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

        self.ax2.clear()
        self.ax2.plot(M_Rdz_positive_list, N_Rd_positive_list, label = 'Interaction envelope - major axis - positive moment', color = '#006D62')
        self.ax2.plot(M_Rdz_negative_list, N_Rd_negative_list, label = 'Interaction envelope - major axis - negative moment', linestyle = '--', color = '#802628')
        self.ax2.scatter(Mz_02,N_Ed, label = 'First Order design action effects', color = '#88BBC2')
        self.ax2.scatter(M_Edz,N_Ed, label = 'Second Order design action effects', color = '#C2B658')
        self.ax2.axhline(0, color='black', linewidth=0.5)  # Horizontal line at y = 0
        self.ax2.axvline(0, color='black', linewidth=0.5)  # Vertical line at x = 0
        self.ax2.grid(True, which='major', linestyle='-', linewidth=0.4, color='gray', alpha=0.5)
        self.ax2.spines['left'].set_position(('data', 0)) # Move the left spine to x=0
        self.ax2.spines['right'].set_color('none') # Hide the right spine
        self.ax2.spines['bottom'].set_position(('data', 0)) # Move the bottom spine to y=0
        self.ax2.spines['top'].set_color('none') # Hide the top spine
        self.ax2.xaxis.set_ticks_position('bottom')
        self.ax2.yaxis.set_ticks_position('left')
        self.ax2.tick_params(axis='y', direction='inout', length=6)
        self.ax2.tick_params(axis='x', rotation=90)
        self.ax2.set_xlabel('M [kNm]', loc='right')  # Move the x-axis label to the right
        self.ax2.set_ylabel('N [kN]', loc='top')
        self.canvas2.draw()
        
    def collect_user_input(self):
        shape = self.shape_var.get()
        N_Ed = float(self.N_Ed_entry.get())
        M_y_top = float(self.M_y_top_entry.get()) 
        M_y_bottom = float(self.M_y_bottom_entry.get())
        M_z_top = float(self.M_z_top_entry.get())
        M_z_bottom = float(self.M_z_bottom_entry.get())
        concrete_grade = self.concrete_dropdown.get()
        concrete_properties = get_concrete_properties(concrete_grade)
        L_eff_y = float(self.L_Effy_entry.get())
        L_eff_z = float(self.L_effz_entry.get())
        f_yk = float(self.f_yk_entry.get())
        E_s = float(self.E_s_entry.get())
        bar_dia = float(self.bar_dia_entry.get())
        link_dia = float(self.link_dia_entry.get())
        if shape == "rectangular":
            h = float(self.h_input.get())
            b = float(self.b_input.get())
            n_x = int(self.n_x_input.get())
            n_y = int(self.n_y_input.get())
            section = Section(shape = shape, h=h, b=b)
            cover = float(self.cover_input.get())
            rebar = Reinforcement(f_yk, E_s, bar_dia, link_dia, cover = cover, shape=shape, b = b, h = h, num_of_rows_of_rebar = n_y, num_of_cols_of_rebar = n_x )
        if shape == "circular":
            diameter = float(self.diameter_entry.get())
            radial_number_bars = int(self.radial_num_bars_input.get())
            cover = float(self.cover_input.get())
            rebar = Reinforcement(f_yk, E_s, bar_dia, link_dia, cover = cover, shape=shape, diameter = diameter, radial_number = radial_number_bars )
            section = Section(shape=shape, diameter = diameter)
        if shape == "arbitrary":
            vertices = self.coord_list
            rebar_coords = self.bar_list
            section = Section(shape = shape, vertices = vertices)
            rebar = Reinforcement(f_yk, E_s, bar_dia, link_dia, shape=shape, rebar_coords = rebar_coords)
        column = Column(section, rebar, concrete_properties, L_eff_y, L_eff_z)
        return N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column

def get_concrete_properties(concrete_grade, gamma_c = 1.5, alpha_cc=0.85):
    #extract material properties from database
    df = pd.read_csv('concrete_grade_DB.csv')
    f_ck = (df.loc[df["Material Property"] == 'f_ck [MPa]', concrete_grade].values[0]) #MPa
    f_ck_cube = (df.loc[df["Material Property"] == 'f_ck_cube [MPa]', concrete_grade].values[0]) #MPa
    f_cm = (df.loc[df["Material Property"] == 'f_cm [MPa]', concrete_grade].values[0]) #MPa
    f_ctm = (df.loc[df["Material Property"] == 'f_ctm [MPa]', concrete_grade].values[0]) # MPa
    f_ctk_0_05 = (df.loc[df["Material Property"] == 'f_ctk_0_05 [MPa]', concrete_grade].values[0]) #MPa
    f_ctk_0_95 = (df.loc[df["Material Property"] == 'f_ctk_0_95 [MPa]', concrete_grade].values[0]) #MPa
    E_cm = (df.loc[df["Material Property"] == 'E_cm [GPa]', concrete_grade].values[0]) #GPa
    eps_c1 = (df.loc[df["Material Property"] == 'strain_c1 [%]', concrete_grade].values[0])*0.001
    eps_cu1 = (df.loc[df["Material Property"] == 'strain_cu1 [%]', concrete_grade].values[0])*0.001
    eps_cu2 = (df.loc[df["Material Property"] == 'strain_cu2 [%]', concrete_grade].values[0])*0.001
    eps_c3 = (df.loc[df["Material Property"] == 'strain_c3 [%]', concrete_grade].values[0])*0.001
    eps_cu3 = (df.loc[df["Material Property"] == 'strain_cu3 [%]', concrete_grade].values[0])*0.001
    if f_ck <50:
        eta = 1.0
    else:
        eta = 1.5
    concrete_properties = Concrete_Material(concrete_grade, f_ck, f_ck_cube, f_cm, f_ctm, f_ctk_0_05, f_ctk_0_95, E_cm, eps_c1, eps_c3, eps_cu1, eps_cu2, eps_cu3, eta, gamma_c, alpha_cc)
    return concrete_properties

def moment_capacity(x, A_s, h, b, d, d_2, E_s, f_yd, eps_cu2, eps_c3, f_cd, lambd):
    d_c = min(lambd*x, h)
    if x <= h:
        in_section = True
    else:
        in_section = False
    #Define concrete strain according to neutral axis location
    if in_section == True:
        concrete_strain = eps_cu2
        eps_sc = concrete_strain * (1 - d_2/x)
        eps_st = concrete_strain * (1 - d/x)
    else: 
        concrete_strain = eps_c3 * x / (x - h/2)
        eps_sc = concrete_strain * (x - d_2) / (x - h/2)
        eps_st = concrete_strain * (x - d) / (x - h/2)
    sigma_sc = max(min(eps_sc * E_s, f_yd),- f_yd)
    sigma_st = max(min(eps_st * E_s, f_yd), - f_yd)
    N = (f_cd*b*d_c + A_s/2 * sigma_sc + A_s / 2 * sigma_st)/10**3
    My_Rd = (f_cd*b*d_c*(x - d_c/2) + A_s/2 *sigma_sc*(x - d_2) - A_s/2*(d - x)*sigma_st)/ 10**6 + N*(h/2 - x) / 10**3
    
    return -My_Rd

def determine_max_moment(A_s, h, b, d, d_2, E_s, f_yd, eps_cu2, eps_c3, f_cd, lambd):
    result = scipy.optimize.minimize(moment_capacity, x0 = 200,bounds = [(1,h)], args=(A_s, h, b, d, d_2, E_s, f_yd, eps_cu2, eps_c3, f_cd, lambd))  # Minimize the negative of the function to find maximum
    # Extract the maximum value and corresponding argument
    max_value = -result.fun  # The result.fun attribute contains the minimum value
    arg_max = result.x  # The result.x attribute contains the value of x at the maximum

def integrate_part_of_circle(neutral_axis, diameter, lambd):
    r = diameter / 2 
    distance_from_centre = diameter / 2 - lambd * neutral_axis
    if distance_from_centre < diameter/2:
        # Calculate the area using the integral
        result, error = quad(lambda x, r: np.sqrt(r**2 - distance_from_centre**2), distance_from_centre, r, args=(r,))
        area = 2 * result
        # Calculate the angle for the centroid calculation
        theta = 2 * np.arccos(distance_from_centre/ r)
        area = r**2/2 *(theta - np.sin(theta))
        # Calculate the centroid y-coordinate (only y because x is 0 / not applicable)
        centroid = diameter/2 - (4 * r * (np.sin(theta / 2) ** 3)) / (3 * (theta - np.sin(theta)))
        lever_arm = neutral_axis - centroid
    else:
        area = math.pi*diameter**2/4
        lever_arm = neutral_axis - diameter/2
        centroid_y = diameter/2
    return area, centroid, lever_arm

def compute_second_moment_area(polygon: Polygon):
    c_x, c_y = polygon.centroid.x, polygon.centroid.y
    coords = list(polygon.exterior.coords)
    Ix = 0.0
    Iy = 0.0
    Ixy = 0.0

    for i in range(len(coords) - 1):
        x0, y0 = coords[i]
        x1, y1 = coords[i + 1]

        Ix += 1/12* (y0 - y1)*(x1+x0-2*c_x)*(x1**2+x0**2+2*c_x*(c_x-x0-x1))        
        Iy += 1/12* (x0 - x1)*(y1+y0-2*c_y)*(y1**2+y0**2+2*c_y*(c_y-y0-y1))
    return Ix, Iy

# Function to clip polygon at a certain y value
def clip_polygon_at_y(polygon, y_value):
    # Create a horizontal line at y_value
    line = LineString([(-1e10, y_value), (1e10, y_value)])
    # Split the polygon by the line
    result = split(polygon, line)
    #Return multiple geoemtry objects
    return result.geoms

# Function to clip polygon at a certain x value
def clip_polygon_at_x(polygon, x_value):
    # Create a horizontal line at x_value
    line = LineString([(x_value, -1e10), (x_value, 1e10)])
    # Split the polygon by the line
    result = split(polygon, line)
    #Return multiple geoemtry objects
    return result.geoms

def compute_area_and_centroid(polygon):
    # Get the coordinates of the polygon
    x, y = polygon.exterior.coords.xy
    # Convert to numpy arrays for easy calculations
    x = np.array(x)
    y = np.array(y)
    # Area Calculation
    A = 0.5 * np.sum(y[:-1] * x[1:] - y[1:] * x[:-1])
    # Centroid Calculation
    Cx = np.sum((x[:-1] + x[1:]) * (y[:-1] * x[1:] - y[1:] * x[:-1])) / (6*A)
    Cy = np.sum((y[:-1] + y[1:]) * (y[:-1] * x[1:] - y[1:] * x[:-1])) / (6*A)
    return A, (Cx, Cy)

def check_moments(N_Ed, M_top, M_bottom):
    e_i = 60/400 # eccentricity in mm
    
    if abs (M_top) > abs(M_bottom):
        if M_top >0:
            M_02 = M_top + N_Ed * e_i*1e-3
        else:
            M_02 = M_top - N_Ed * e_i*1e-3
        if M_bottom >0:
            M_01 = M_bottom + N_Ed*1e-3
        else:
            M_01 = M_bottom - N_Ed * e_i*1e-3 
    else: 
        if M_top > 0:
            M_01 = M_top + N_Ed * e_i*1e-3
        else:
            M_01 = M_top - N_Ed * e_i*1e-3
        if M_bottom > 0:
            M_02 = M_bottom + N_Ed * e_i*1e-3
        else:
            M_02 = M_bottom - N_Ed * e_i*1e-3
    return M_01, M_02

def check_slenderness (column, N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, phi_eff = None):
    M_01y, M_02y =  check_moments(N_Ed, M_y_top, M_y_bottom)
    M_01z, M_02z =  check_moments(N_Ed, M_z_top, M_z_bottom)
    
    # Check absolute slenderness
    slenderness_ratio_y = column.L_eff_y*1e3 / column.section.i_y
    slenderness_ratio_z = column.L_eff_z*1e3 / column.section.i_z

    if column.section.shape == "rectangular":
        e_0 = max(column.section.h/30, 20)*1e-3 # eccentricity dimension in m 
    elif column.section.shape == "circular":
        e_0 = max(column.section.diameter/30, 20)*1e-3 # eccentricity dimension in m
    elif column.section.shape == "arbitrary":
        e_0  = max((column.section.top_of_section - column.section.bottom_of_section)/30, 20)*1e-3 # eccentricity dimension in m
    # Check slenderness limits
    print(f"reinforcement area: {column.reinforcement.A_s}")
    mechanical_reinforcement_ratio = (column.reinforcement.A_s*column.reinforcement.f_yd) / (column.section.A*column.concrete_properties.f_cd)
    print(f'reinforcement ratio: {mechanical_reinforcement_ratio}')
    if phi_eff:
        A = 1 / (1+0.2*phi_eff)
    else:
        A = 0.7
    B = (1+2*mechanical_reinforcement_ratio)**0.5
    r_m_y = M_01y / M_02y
    r_m_z = M_01z / M_02z
    C_y = 1.7 - r_m_y
    C_z = 1.7 - r_m_z
    n = N_Ed / (column.section.A * column.concrete_properties.f_cd * 1e-3)  # axial force utilisation
    print(f'A: {A}')
    print(f'B: {B}')
    print(f'C_y: {C_y}')
    print(f'C_z: {C_z}')
    print(f'n: {n}')
    slenderness_limit_y = 20*A*B*C_y / (n**0.5)  
    print(f'slenderness ratio [y]: {slenderness_ratio_y}')
    print(f'slenderness limit [y]: {slenderness_limit_y}')
    slenderness_limit_z = 20*A*B*C_z / (n**0.5)
    print(f'slenderness ratio [z]: {slenderness_ratio_z}')
    print(f'slenderness limit [z]: {slenderness_limit_z}')
    # return major axis moment dependant on slenderness in major axis
    if slenderness_ratio_y > slenderness_limit_y:
        slenderness_y = True
        print("Column is slender is about the major axis")
        M_Edy = M_02y/abs(M_02y) * compute_major_axis_slender_moments(column, slenderness_ratio_y, N_Ed, n, M_01y, M_02y, phi_eff = 2.15) #MEdy is returned in kNm
    else:
        slenderness_y = False
        M_Edy = M_02y/abs(M_02y) * max(M_02y, N_Ed * e_0)
        print("Column is not slender about major axis")
    # return minor axis moment dependant on slenderness in minor axis
    if slenderness_ratio_z > slenderness_limit_z:
        slenderness_z = True
        print("Column is slender is about the minor axis")
        M_Edz = M_02z/abs(M_02z) *compute_minor_axis_slender_moments(column, slenderness_ratio_z, N_Ed, n, M_01z, M_02z, phi_eff = 2.15)
    else:
        slenderness_z = False
        print("Column is not slender about minor axis")
        M_Edz = M_02y/abs(M_02y) * max(M_02y, N_Ed * e_0)
    return M_01y, M_02y, M_01z, M_02z, M_Edy, M_Edz, slenderness_y, slenderness_z, slenderness_ratio_y, slenderness_ratio_z

def compute_major_axis_slender_moments(column, slenderness_ratio, N_Ed, n, M_01y, M_02y, phi_eff = 2.15):
    print(f"reinforcement area: {column.reinforcement.A_s}")
    mechanical_reinforcement_ratio = (column.reinforcement.A_s*column.reinforcement.f_yd) / (column.section.A*column.concrete_properties.f_cd) 
    n_u = 1 + mechanical_reinforcement_ratio
    n_bal = 0.4
    K_r = min((n_u - n)/(n_u - n_bal), 1)
    print(f'K_r: {K_r}')
    if column.section.shape == "rectangular":
        e_0 = max(column.section.h/30, 20)*1e-3 # eccentricity dimension in m 
    elif column.section.shape == "circular":
        e_0 = max(column.section.diameter/30, 20)*1e-3 # eccentricity dimension in m
    elif column.section.shape == "arbitrary":
        e_0  = max((column.section.top_of_section - column.section.bottom_of_section)/30, 20)*1e-3 # eccentricity dimension in m

    beta = 0.35 + column.concrete_properties.f_ck/200 - slenderness_ratio/150
    K_phi = max (1, (1+beta*phi_eff))
    e_2y = (0.1* (K_r * K_phi * column.reinforcement.f_yd) / (0.45 * column.dy * column.reinforcement.E_s*1e3) * (column.L_eff_y*1e3)**2) *1e-3 # eccentricity dimension in m
    M_0e_y = 0.6*M_02y + 0.4*M_01y
    M_2y = N_Ed * e_2y
    M_Ed_y = max(M_02y, (M_0e_y + M_2y), (M_01y + 0.2*M_2y), N_Ed*e_0)
    
    return M_Ed_y

def compute_minor_axis_slender_moments(column, slenderness_ratio, N_Ed, n, M_01z, M_02z, phi_eff = 2.15):
    mechanical_reinforcement_ratio = (column.reinforcement.A_s*column.reinforcement.f_yd) / (column.section.A*column.concrete_properties.f_cd)
    n_u = 1 + mechanical_reinforcement_ratio
    n_bal = 0.4
    if (n_u - n)/(n_u - n_bal) >=1:
        K_r = 1.0
    else:
        K_r = (n_u - n)/(n_u - n_bal)
    
    if column.section.shape == "rectangular":
        e_0 = max(column.section.h/30, 20)*1e-3 # eccentricity dimension in m 
    elif column.section.shape == "circular":
        e_0 = max(column.section.diameter/30, 20)*1e-3 # eccentricity dimension in m
    elif column.section.shape == "arbitrary":
        e_0  = max((column.section.top_of_section - column.section.bottom_of_section)/30, 20)*1e-3 # eccentricity dimension in m

    beta = 0.35 + column.concrete_properties.f_ck/200 - slenderness_ratio/150
    K_phi = max (1, (1+beta*phi_eff))
    e_2z = (0.1* (K_r*K_phi*column.reinforcement.f_yd/(0.45*column.dz*column.reinforcement.E_s*1e3))*(column.L_eff_z*1e3)**2) *1e-3 # eccentricity dimension in m
    M_0e_z = 0.6*M_02z + 0.4*M_01z
    M_2z = N_Ed * e_2z
    M_Ed_z = max(M_02z, (M_0e_z + M_2z), (M_01z + 0.2*M_2z), N_Ed*e_0)
    return M_Ed_z

def compute_creep_coefficients(column, RH, M_Ed_y, M_Ed_z, M_0Eqp_y, M_0Eqp_z, t_0, h_0):
    f_cm = column.concrete_properties.f_cm
    alpha_1 = (35/f_cm)**0.7
    alpha_2 = (35/f_cm)**0.2
    alpha_3 = (35/f_cm)**0.5
    if f_cm <=35:
        phi_RH = 1 + (1-(RH/100)) / (0.1*h_0**(1/3))
    else:
        phi_RH = (1 + (1-(RH/100)) / (0.1*h_0**(1/3))* alpha_1 )* alpha_2
    beta_fcm = 16.8 / f_cm**0.5
    beta_t_0 =  1 / (0.1 + t_0**0.2)
    h_0 = 2*column.section.A / column.section.polygon.length
    phi_0 = phi_RH * beta_fcm * beta_t_0
    phi_inf_t0 = phi_RH * beta_fcm * beta_t_0  
    phi_effy = phi_inf_t0 * M_0Eqp_y / M_Ed_y
    phi_effz = phi_inf_t0 * M_0Eqp_z / M_Ed_z

    return phi_effy, phi_effz
    
def determine_envelope_value_major_axis_positive(column, lambd, neutral_axis_y):       # Determines envelope value for major axis bending
    N_Rd = 0 # Axial capacity in kN
    M_Rdy = 0 # Moment capacity in kNm

    if column.section.shape == 'rectangular':
        section_centroid = [column.section.b/2, column.section.h/2]
        d_cy = min(lambd*neutral_axis_y, column.section.h)
        if neutral_axis_y <= column.section.h:
            in_section = True
            lever_arm_concrete_y = neutral_axis_y * (1 - lambd/2)
            N_Rd += d_cy * column.section.b * column.concrete_properties.f_cd *1e-3
            M_Rdy += d_cy * column.section.b * column.concrete_properties.f_cd * lever_arm_concrete_y *1e-6
        else:
            in_section = False
            lever_arm_concrete_y = (neutral_axis_y - d_cy/2)
            N_Rd += d_cy * column.section.b * column.concrete_properties.f_cd *1e-3
            M_Rdy += d_cy * column.section.b * column.concrete_properties.f_cd * lever_arm_concrete_y * 1e-6
    
    elif column.section.shape == 'circular':
        section_centroid = [column.section.diameter/2, column.section.diameter/2]
        if neutral_axis_y <= column.section.diameter:
            in_section = True
            segment_area, centroid_y, lever_arm_concrete_y = integrate_part_of_circle(neutral_axis_y, column.section.diameter, lambd) # area in mm^2, centroid distance in mm
            N_Rd += segment_area*column.concrete_properties.f_cd *1e-3
            M_Rdy += segment_area*column.concrete_properties.f_cd*lever_arm_concrete_y *1e-6
        else:
            in_section = False
            lever_arm_concrete_y = neutral_axis_y - section_centroid[1]
            N_Rd += column.section.A*column.concrete_properties.f_cd * 1e-3
            M_Rdy += column.section.A*column.concrete_properties.f_cd * lever_arm_concrete_y *1e-6
        
        #Define concrete strain according to neutral axis location
    
    elif column.section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.section.polygon)[1]
        if neutral_axis_y <= column.section.h:
            in_section = True
        else:
            in_section = False
        
        clipped_concrete = clip_polygon_at_y(column.section.polygon, column.section.bottom_of_section + lambd*neutral_axis_y)
        #cycle through the geometries created by clipping the arbitrary shape
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                # determine area and centroid of these shapes created
                area, centroid = compute_area_and_centroid(geom) # distance to centroid in mm, are in mm^2
                # determine whether the centroid is above the neutral axis (i.e. in compression)
                if centroid[1] < (column.section.bottom_of_section + neutral_axis_y):
                    N_Rd += area * column.concrete_properties.f_cd *1e-3
                    M_Rdy += (column.section.bottom_of_section + neutral_axis_y - centroid[1]) * area * column.concrete_properties.f_cd *1e-6
    
    steel_strains = []
    steel_stresses = []

    # determnine concrete strain based on whether neutral axis is within section or not
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2
    else:
        concrete_strain = column.concrete_properties.eps_c3*neutral_axis_y/(neutral_axis_y - section_centroid[1])
    
    # loop through steel reinforcement, determine stresses and add to moment capacity / axial capacity
    for rebar_coords in column.reinforcement.arrangement:
        if in_section == True:
            if column.section.shape == "rectangular" or column.section.shape == "circular":
                steel_strain = (neutral_axis_y - rebar_coords[1])/neutral_axis_y * concrete_strain
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * ((column.section.bottom_of_section + neutral_axis_y) - rebar_coords[1]) / (neutral_axis_y)   # Steel strain computed if neutral axis is within the section. Sim triangles 
        else:
            if column.section.shape == "rectangular":
                steel_strain = (neutral_axis_y - rebar_coords[1]) / (neutral_axis_y) * concrete_strain
            elif column.section.shape == "circular":
                steel_strain = concrete_strain * (neutral_axis_y - rebar_coords[1]) / (neutral_axis_y - column.section.diameter/2)
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * ((column.section.bottom_of_section + neutral_axis_y) - rebar_coords[1]) / (neutral_axis_y - section_centroid[1])

        steel_stress = max(min(steel_strain*column.reinforcement.E_s*1e3, column.reinforcement.f_yd), -column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain)
        steel_stresses.append(steel_stress) # steel stress in MPa
        if column.section.shape == "rectangular":
            lever_arm_y = neutral_axis_y - rebar_coords[1]
        elif column.section.shape == "circular":
            lever_arm_y = neutral_axis_y - rebar_coords[1]
        elif column.section.shape == "arbitrary":
            lever_arm_y = column.section.bottom_of_section + neutral_axis_y - rebar_coords[1]# lever arm from neutral axis in mm

        N_Rd += steel_stress * math.pi * (column.reinforcement.bar_diameter)**2 /4 * 1e-3 #kN
        M_Rdy += steel_stress * math.pi * (column.reinforcement.bar_diameter)**2 /4 * lever_arm_y *1e-6 #kNm

    # compute moment contribution from axial force acting on section
    if column.section.shape == "rectangular":
        M_Rdy += N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3 #kNm
    elif column.section.shape == "circular":
        M_Rdy += N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3 #kNm
    elif column.section.shape == "arbitrary":
        M_Rdy += N_Rd*(section_centroid[1] - (column.section.bottom_of_section + neutral_axis_y)) * 1e-3 #kNm

    return N_Rd, M_Rdy, steel_stresses, steel_strains

def determine_envelope_value_major_axis_negative(column, lambd, neutral_axis_y):       # Determines envelope value for major axis bending
    N_Rd = 0 # Axial capacity in kN
    M_Rdy = 0 # Moment capacity in kNm
    
    # DEFINE NEUTRAL AXIS FROM OTHER DIRECTION IN ORDER TO GET NEGATIVE MOMENTS
    if column.section.shape == 'rectangular':
        section_centroid = [column.section.b/2, column.section.h/2]
        d_cy = min(lambd * neutral_axis_y, column.section.h)
        lever_arm_concrete_y = d_cy/2 - neutral_axis_y
        N_Rd += d_cy * column.section.b * column.concrete_properties.f_cd *1e-3
        if (column.section.h - neutral_axis_y) > 0:
            in_section = True
            M_Rdy += d_cy * column.section.b * column.concrete_properties.f_cd * lever_arm_concrete_y *1e-6
        else:
            in_section = False
            M_Rdy += d_cy * column.section.b * column.concrete_properties.f_cd * lever_arm_concrete_y * 1e-6
    
    elif column.section.shape == 'circular':
        section_centroid = [column.section.diameter/2, column.section.diameter/2]
        if (column.section.diameter - neutral_axis_y) > 0:
            in_section = True
            segment_area, centroid_y, lever_arm_y = integrate_part_of_circle(neutral_axis_y, column.section.diameter, lambd) # area in mm^2, centroid distance in mm
            N_Rd += segment_area* column.concrete_properties.f_cd *1e-3
            M_Rdy += - lever_arm_y * segment_area * column.concrete_properties.f_cd *1e-6
        else:
            in_section = False
            N_Rd += column.section.A*column.concrete_properties.f_cd * 1e-3
            M_Rdy += column.section.A*column.concrete_properties.f_cd * (section_centroid[0] - neutral_axis_y) *1e-6
    
    elif column.section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.section.polygon)[1]
        if (column.section.top_of_section - neutral_axis_y) > column.section.bottom_of_section:
            in_section = True
        else:
            in_section = False
        
        clipped_concrete = clip_polygon_at_y(column.section.polygon, (column.section.top_of_section - lambd * neutral_axis_y))
        # cycle through the geometries created by clipping the arbitrary shape
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                # determine area and centroid of these shapes created
                area, centroid = compute_area_and_centroid(geom) # distance to centroid in mm, are in mm^2
                # determine whether the centroid is above the neutral axis (i.e. in compression)
                if centroid[1] > (column.section.top_of_section - neutral_axis_y):
                    N_Rd += area * column.concrete_properties.f_cd *1e-3
                    lever_arm_y = (column.section.top_of_section - neutral_axis_y - centroid[1])
                    M_Rdy += lever_arm_y * area * column.concrete_properties.f_cd *1e-6
    
    steel_strains = []
    steel_stresses = []

    # determnine concrete strain based on whether neutral axis is within section or not
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2
    else:
        concrete_strain = column.concrete_properties.eps_c3*neutral_axis_y/(neutral_axis_y - section_centroid[1])
    
    # loop through steel reinforcement, determine stresses and add to moment capacity / axial capacity
    for rebar_coords in column.reinforcement.arrangement:
        if in_section == True:
            if column.section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.section.h - neutral_axis_y)) / neutral_axis_y
            elif column.section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.section.diameter - neutral_axis_y)) / neutral_axis_y
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.section.top_of_section - neutral_axis_y)) / neutral_axis_y
        else:
            if column.section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.section.h - neutral_axis_y)) / (neutral_axis_y)
            elif column.section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.section.diameter - neutral_axis_y)) / (neutral_axis_y)
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.section.top_of_section - neutral_axis_y)) / (neutral_axis_y - section_centroid[1])
        
        steel_stress = max(min(steel_strain*column.reinforcement.E_s*1e3, column.reinforcement.f_yd), - column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain) 
        steel_stresses.append(steel_stress) # steel stress in MPa

        if column.section.shape == "rectangular":
            lever_arm_y = (column.section.h - neutral_axis_y) - rebar_coords[1] 
        elif column.section.shape == "circular":
            lever_arm_y = (column.section.diameter - neutral_axis_y) - rebar_coords[1]
        elif column.section.shape == "arbitrary":
            lever_arm_y = (column.section.top_of_section - neutral_axis_y) - rebar_coords[1] # lever arm from neutral axis in mm
        
        N_Rd += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * 1e-3 #kN
        M_Rdy += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * lever_arm_y *1e-6 #kNm

    # compute moment contribution from axial force acting on section
    if column.section.shape == "rectangular":
        M_Rdy += N_Rd*(section_centroid[1] - (column.section.h - neutral_axis_y)) * 1e-3
    elif column.section.shape == "circular":
        M_Rdy += N_Rd*(section_centroid[1] - (column.section.diameter - neutral_axis_y)) * 1e-3        
    elif column.section.shape == "arbitrary":
        M_Rdy += N_Rd*(section_centroid[1] - (column.section.top_of_section - neutral_axis_y)) * 1e-3
    return N_Rd, M_Rdy, steel_stresses, steel_strains

def determine_envelope_value_minor_axis_positive(column, lambd, neutral_axis_x):       # Determines envelope value for major axis bending
    N_Rd = 0 #Axial capacity in kN
    M_Rdz = 0 #Moment capacity in kNm

    if column.section.shape == 'rectangular':
        section_centroid = [column.section.b/2, column.section.h/2]
        d_cx = min(lambd*neutral_axis_x, column.section.b)
        lever_arm_concrete_x = neutral_axis_x  - d_cx/2                                                # lever arm in mm
        N_Rd += d_cx * column.section.h * column.concrete_properties.f_cd *1e-3
        if neutral_axis_x <= column.section.b:
            in_section = True                              # Add axial force from concrete in compression
            M_Rdz += d_cx * column.section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6      # Add moment from concrete in compression acting at a lever arm from the neutral axis
        else:
            in_section = False                                           # lever arm in mm
            M_Rdz += d_cx * column.section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
    
    elif column.section.shape == 'circular':
        section_centroid = [column.section.diameter/2, column.section.diameter/2]
        if neutral_axis_x <= column.section.diameter:
            in_section = True
            segment_area, centroid_x, lever_arm_concrete_x = integrate_part_of_circle(neutral_axis_x, column.section.diameter, lambd) # area in mm^2, centroid distance in mm
            N_Rd += segment_area*column.concrete_properties.f_cd *1e-3
            M_Rdz += segment_area*column.concrete_properties.f_cd*lever_arm_concrete_x *1e-6
        else:
            in_section = False
            lever_arm_concrete_x = neutral_axis_x - section_centroid[0]
            N_Rd += column.section.A*column.concrete_properties.f_cd * 1e-3
            M_Rdz += column.section.A*column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
        
    
    elif column.section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.section.polygon)[1]                                 # get centroid of polygon of the whole section      
        if column.section.left_of_section + neutral_axis_x < column.section.right_of_section:                   # Determine whether neutral axis is within the section or not
            in_section = True
        else:
            in_section = False
    
        clipped_concrete = clip_polygon_at_x(column.section.polygon, column.section.left_of_section + lambd*neutral_axis_x) #Clip the polygon into multiple polygons
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                area, centroid = compute_area_and_centroid(geom)                                                # area in mm^2, centroid in mm from extreme compression fibre
                if centroid[0] < column.section.left_of_section + neutral_axis_x:                               # Find all the geometries created that are on the side of the NA under compression
                    N_Rd += area * column.concrete_properties.f_cd *1e-3                                        # Add axial force for concrete in compression 
                    M_Rdz += (column.section.left_of_section + neutral_axis_x - centroid[0]) * area * column.concrete_properties.f_cd * 1e-6     # Moment resistance contribution from concrete in compression acting at a lever arm from the neutral axis 
    
    steel_strains = []
    steel_stresses = []
    
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2                                                    # Concrete strain if neutral axis within section
    else:
        concrete_strain = column.concrete_properties.eps_c3*neutral_axis_x/(neutral_axis_x - section_centroid[1])                                                     # Concrete strain if neutral axis outside of the section. Note that this is the concrete strain at the centre of the section, therefore computing steel strains realtive to this should use sim triangle from centre of section
    
    for rebar_coords in column.reinforcement.arrangement:                                                       
        if in_section ==True:
            if column.section.shape == "rectangular" or column.section.shape == "circular":
                steel_strain = (neutral_axis_x - rebar_coords[0])/neutral_axis_x * concrete_strain
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * ((column.section.left_of_section + neutral_axis_x) - rebar_coords[0]) / (neutral_axis_x)   # Steel strain computed if neutral axis is within the section. Sim triangles 
        else:
            if column.section.shape == "rectangular":
                steel_strain = concrete_strain * (neutral_axis_x - rebar_coords[0]) / (neutral_axis_x - column.section.b/2)
            elif column.section.shape == "circular":
                steel_strain = concrete_strain * (neutral_axis_x - rebar_coords[0]) / (neutral_axis_x - column.section.diameter/2)
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * ((column.section.left_of_section + neutral_axis_x) - rebar_coords[0]) / (neutral_axis_x - section_centroid[0])
        
        steel_stress = max(min(steel_strain*column.reinforcement.E_s*1e3, column.reinforcement.f_yd), - column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain)
        steel_stresses.append(steel_stress) # steel stress in MPa
        if column.section.shape == "rectangular":
            lever_arm_x = neutral_axis_x - rebar_coords[0]
        elif column.section.shape == "circular":
            lever_arm_x = neutral_axis_x - rebar_coords[0]
        elif column.section.shape == "arbitrary":
            lever_arm_x = column.section.left_of_section + neutral_axis_x - rebar_coords[0] # lever arm in mm
        
        N_Rd += steel_stress * math.pi*(column.reinforcement.bar_diameter**2)/4 *1e-3
        M_Rdz += steel_stress * math.pi*(column.reinforcement.bar_diameter**2)/4 * lever_arm_x * 1e-6
    
    # compute moment contribution from axial force acting on section
    if column.section.shape == "rectangular":
        M_Rdz += N_Rd*(section_centroid[0] - neutral_axis_x) * 1e-3
    elif column.section.shape == "circular":
        M_Rdz += N_Rd*(section_centroid[0] - neutral_axis_x) * 1e-3
    elif column.section.shape == "arbitrary":
        M_Rdz += N_Rd*(section_centroid[0] - (column.section.left_of_section + neutral_axis_x)) * 1e-3

    return N_Rd, M_Rdz, steel_stresses, steel_strains

def determine_envelope_value_minor_axis_negative(column, lambd, neutral_axis_x):        # Determines envelope value for minor axis bending for positive moments (left half of section in compression)
    N_Rd = 0                                                                            #Initialise axial capacity in kN, value accumulates as cycling through concrete and rebar
    M_Rdz = 0             

    if column.section.shape == 'rectangular':
        section_centroid = [column.section.b/2, column.section.h/2]
        d_cx = min(lambd*neutral_axis_x, column.section.b)
        lever_arm_concrete_x = d_cx/2 - neutral_axis_x #lever arm in mm
        N_Rd += d_cx * column.section.h * column.concrete_properties.f_cd *1e-3
        if (column.section.b - neutral_axis_x) > 0:
            in_section = True
            M_Rdz += d_cx * column.section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
        else:
            in_section = False
            M_Rdz += d_cx * column.section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
    
    elif column.section.shape == 'circular':
        section_centroid = [column.section.diameter/2, column.section.diameter/2]
        if (column.section.diameter - neutral_axis_x) > 0:
            in_section = True
            segment_area, centroid_x, lever_arm_x = integrate_part_of_circle(neutral_axis_x, column.section.diameter, lambd) # area in mm^2, centroid distance in mm
            N_Rd += segment_area* column.concrete_properties.f_cd *1e-3
            M_Rdz += - lever_arm_x * segment_area * column.concrete_properties.f_cd *1e-6
        else:
            in_section = False
            N_Rd += column.section.A*column.concrete_properties.f_cd * 1e-3
            M_Rdz += column.section.A*column.concrete_properties.f_cd * (section_centroid[0] - neutral_axis_x) *1e-6
    
    elif column.section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.section.polygon)[1]
        if (column.section.right_of_section - neutral_axis_x) > column.section.left_of_section:
            in_section = True
        else:
            in_section = False
        
        clipped_concrete = clip_polygon_at_x(column.section.polygon, (column.section.right_of_section - lambd*neutral_axis_x))
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                area, centroid = compute_area_and_centroid(geom) #area in mm^2, centroid in mm from extreme compression fibre
                if centroid[0] > (column.section.right_of_section - neutral_axis_x):
                    N_Rd += area * column.concrete_properties.f_cd *1e-3
                    lever_arm_x = (column.section.right_of_section - neutral_axis_x - centroid[0])
                    M_Rdz +=  lever_arm_x * area * column.concrete_properties.f_cd * 1e-6
    
    steel_strains = []
    steel_stresses = []
    
    # DETERMINE CONCRETE STRAINS BASED ON WHETHER NEUTRAL AXIS IS WITHIN SECTION OR NOT
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2
    else:
        concrete_strain = column.concrete_properties.eps_c3*neutral_axis_x/(neutral_axis_x - section_centroid[0])
    
    for rebar_coords in column.reinforcement.arrangement:
        if in_section == True:
            if column.section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.section.b - neutral_axis_x)) / neutral_axis_x
            elif column.section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.section.diameter - neutral_axis_x)) / neutral_axis_x
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.section.right_of_section - neutral_axis_x))/ neutral_axis_x
        else:
            if column.section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.section.b - neutral_axis_x)) / (neutral_axis_x - section_centroid[0])
            elif column.section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.section.diameter - neutral_axis_x)) / (neutral_axis_x - section_centroid[0])
            elif column.section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.section.right_of_section - neutral_axis_x)) / (neutral_axis_x - section_centroid[0])

        steel_stress = max(min(steel_strain * column.reinforcement.E_s*1e3, column.reinforcement.f_yd), -column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain)
        steel_stresses.append(steel_stress) # steel stress in MPa

        if column.section.shape == "rectangular":
            lever_arm_x = (column.section.b - neutral_axis_x) - rebar_coords[0]
        elif column.section.shape == "circular":
            lever_arm_x = (column.section.diameter - neutral_axis_x) - rebar_coords[0] # lever arm in mm
        elif column.section.shape == "arbitrary":
            lever_arm_x = (column.section.right_of_section - neutral_axis_x) - rebar_coords[0] # lever arm in mm

        N_Rd += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * 1e-3
        M_Rdz += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * lever_arm_x * 1e-6
        
    # compute moment contribution from axial force acting on section
    if column.section.shape == "rectangular":
        M_Rdz += N_Rd*(section_centroid[0] - (column.section.b - neutral_axis_x)) * 1e-3
    elif column.section.shape == "circular":
        M_Rdz += N_Rd*(section_centroid[0] - (column.section.diameter - neutral_axis_x)) * 1e-3
    elif column.section.shape == "arbitrary":
        M_Rdz += N_Rd*(section_centroid[0] - (column.section.right_of_section - neutral_axis_x)) * 1e-3
    
    return N_Rd, M_Rdz, steel_stresses, steel_strains

root = tk.Tk()
root.title('LABS - RC Column design Tool')
#root.iconbitmap(r'./assets/Logo.ico')
#root.tk.call("source", "azure.tcl")
#root.tk.call("set_theme", "dark")
app = InputApp(root)
root.mainloop()
root.destroy()