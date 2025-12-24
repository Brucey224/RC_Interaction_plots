import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import math
from .utils import get_concrete_properties
from .plotting import plot_rectangular_section, plot_circular_section, plot_arbitrary_section, plot_major_axis_failure_envelope, plot_minor_axis_failure_envelope
from .codeChecks import check_slenderness
from .classes import Column, Concrete_Section, Reinforcement, Concrete_Material

class InputApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RC Column Designer")
        
        self.slenderness_values = {
            "slenderness_x": None,
            "slenderness_y": None,
            "slenderness_ratio_x": None,
            "slenderness_ratio_y": None,
        }
        self.shape = None
        self.b = None
        self.h = None
        self.diameter = None
        self.vertices = None
        self.cover = None
        self.n_x = None
        self.n_y = None
        self.radial_num_bars = None
        self.concrete_grade = None
        self.f_yk = None
        self.E_s = None
        self.link_dia = None
        self.bar_dia = None
        self.L_effy = None
        self.L_effz = None
        self.N_Ed = None
        self.M_y_top = None
        self.M_y_bottom = None
        self.M_z_top = None
        self.M_z_bottom = None

        # Radio Button Variable
        self.shape_var = tk.StringVar(value="rectangular")
        
        self.coord_list = []
        self.bar_list = []

        # Create the main layout
        self.create_main_layout()

    def update_variables_and_plot(self, event):
        # Validate and update variables based on the selected shape
        self.shape = self.shape_var.get()       
        try:
            self.link_dia = int(self.link_dia_entry.get()) if hasattr(self, "link_dia_entry") else None
            self.bar_dia = int(self.bar_dia_entry.get()) if hasattr(self, "bar_dia_entry") else None
            self.f_yk = float(self.f_yk_entry.get()) if hasattr(self, "f_yk_entry") else None
            self.E_s = float(self.E_s_entry.get()) if hasattr(self, "E_s_entry") else None
        except ValueError:
            print("Bar parameters must be valid inputs.")
        try:
            self.L_effy = float(self.L_effy_entry.get()) if hasattr(self, "L_effy_entry") else None
            self.L_effz = float(self.L_effz_entry.get()) if hasattr(self, "L_effz_entry") else None
        except ValueError:
            print("Effective lengths must be numeric values.")
        try:
            self.N_Ed = float(self.N_Ed_entry.get()) if hasattr(self, "N_Ed_entry") else None
            self.M_y_top = float(self.M_y_top_entry.get()) if hasattr(self, "M_y_top_entry") else None
            self.M_y_bottom = float(self.M_y_bottom_entry.get()) if hasattr(self, "M_y_bottom_entry") else None
            self.M_z_top = float(self.M_z_top_entry.get()) if hasattr(self, "M_z_top_entry") else None
            self.M_z_bottom = float(self.M_z_bottom_entry.get()) if hasattr(self, "M_z_bottom_entry") else None
        except ValueError:
            print("make sure forces are numeric values.")
        
        if self.shape == "rectangular":
            try:
                self.b = int(self.b_entry.get()) if hasattr(self, "b_entry") else None
                self.h = int(self.h_entry.get()) if hasattr(self, "h_entry") else None
            except ValueError:
                print("Error: Section width (b) and height (h) must be integers.")
                return 
            # Validate that n_x and n_y are integers
            try:
                self.n_x = int(self.n_x_entry.get())
                self.n_y = int(self.n_y_entry.get())
            except ValueError:
                print("Error: Number of bars (n_x and n_y) must be integers.")
                return 
            try:
                self.cover = float(self.cover_entry.get())
            except ValueError:
                print("Error: Cover must be a numeric value.")
                return 
            
        elif self.shape == "circular":
            try:
                self.diameter = float(self.diameter_entry.get()) if hasattr(self, "diameter_entry") else None
                self.cover = float(self.cover_entry.get()) if hasattr(self, "cover_entry") else None
            except ValueError:
                print("Error: Ensure cover and dimaeters are numeric values.")
                return
            try:
                self.radial_num_bars = int(self.radial_num_bars_entry.get()) if hasattr(self, "radial_num_bars_entry") else None
                print(f"radial num bars: {self.radial_num_bars}")
            except ValueError:
                print("Error: Number of radial bars must be an integer")
                return      
               
        elif self.shape == "arbitrary":
            self.vertices = self.vertices_input.get("1.0", tk.END).strip()
            self.rebar_coords = self.rebar_coords_input.get("1.0", tk.END).strip()


        type_mapping = {
            "shape": str,
            "b": float,
            "h": float,
            "diameter": float,
            "cover": float,
            "n_x": int,
            "n_y": int,
            "radial_num_bars": int,
            "concrete_grade": str,
            "f_yk": float,
            "E_s": float,
            "link_dia": float,
            "bar_dia": float,
            "L_effy": float,
            "L_effz": float,
            "N_Ed": float,
            "M_y_top": float,
            "M_y_bottom": float,
            "M_z_top": float,
            "M_z_bottom": float
        }
        
        # Convert values to their expected types
        for attribute in type_mapping:
            value = getattr(self, attribute, None)
            if value is not None:
                try:
                    converted_value = type_mapping[attribute](value)
                    setattr(self, attribute, converted_value)
                    print(f"Updated {attribute} to {converted_value}")  # Debugging
                except (ValueError, TypeError):
                    print(f"Failed to convert {attribute} to {type_mapping[attribute].__name__}. Value: {value}")

        # Call the appropriate plotting function based on the shape
        if self.shape == "rectangular":
            plot_rectangular_section(self.ax0, self.canvas0,self.b,self.h,self.cover,self.link_dia,self.bar_dia,self.n_x,self.n_y)
        elif self.shape == "circular":
            plot_circular_section(self.ax0,self.canvas0,self.diameter,self.radial_num_bars,self.cover )
        elif self.shape == "arbitrary":
            self.update_section_plot()  # For arbitrary shapes, use the existing method

    def create_main_layout(self):
        # Options input for concrete grade
        options = ["C12/15","C16/20","C20/25","C25/30","C30/37","C35/45","C40/50","C45/55","C50/60","C55/67","C60/75","C70/85","C80/95","C90/105"]
        #Create row for concrete choice
        concrete_choice_label = ttk.Label(self.root, text="Concrete Grade:")
        concrete_choice_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        
        # Create the Combobox
        self.concrete_dropdown = ttk.Combobox(self.root, values=options)
        self.concrete_dropdown.current()  # Set the default option (index 0)
        self.concrete_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky="w")

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
        self.shape_label = ttk.Label(self.root, text="Section shape:")
        self.shape_label.grid(row=5, column=0, padx=10, pady=10, sticky="w")

        # --- Add Right-Hand Plots below the buttons ---
        # Create a frame for the two graphs below the buttons
        
        self.lengthloads_frame = ttk.Frame(self.root)
        self.lengthloads_frame.grid(row=0, column=6, rowspan=5, columnspan=5, padx=0, pady=0, sticky="nsew")

        self.L_effy_label = ttk.Label(self.lengthloads_frame, text="L_eff,y:")
        self.L_effy_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")
        self.L_effy_entry = ttk.Entry(self.lengthloads_frame)
        self.L_effy_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")
        self.L_effy_units = ttk.Label(self.lengthloads_frame, text="m")
        self.L_effy_units.grid(row=0, column=2, padx=10, pady=5, sticky="w")
        self.L_effy_entry.bind("<Return>", self.update_variables_and_plot)
        self.L_effy_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.L_effz_label = ttk.Label(self.lengthloads_frame, text="L_eff,z:")
        self.L_effz_label.grid(row=0, column=3, padx=10, pady=5, sticky="e")
        self.L_effz_entry = ttk.Entry(self.lengthloads_frame)
        self.L_effz_entry.grid(row=0, column=4, padx=10, pady=5, sticky="ew")
        self.L_effz_units = ttk.Label(self.lengthloads_frame, text="m")
        self.L_effz_units.grid(row=0, column=5, padx=10, pady=5, sticky="w")
        self.L_effz_entry.bind("<Return>", self.update_variables_and_plot)
        self.L_effz_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.N_Ed_label = ttk.Label(self.lengthloads_frame, text="N_Ed:")
        self.N_Ed_label.grid(row=2, column=0, padx=10, pady=5, sticky="e")
        self.N_Ed_entry = ttk.Entry(self.lengthloads_frame)
        self.N_Ed_entry.grid(row=2, column=1, padx=10, pady=5, sticky="ew")
        self.N_Ed_units = ttk.Label(self.lengthloads_frame, text="kN")
        self.N_Ed_units.grid(row=2, column=2, padx=10, pady=5, sticky="w")
        self.N_Ed_entry.bind("<Return>", self.update_variables_and_plot)
        self.N_Ed_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.First_order_moments_label = ttk.Label(self.lengthloads_frame, text="First Order Moments:")
        self.First_order_moments_label.grid(row=3, column=0, padx=10, pady=5, sticky="w")

        self.M_y_top_label = ttk.Label(self.lengthloads_frame, text="M_y [top]:")
        self.M_y_top_label.grid(row=4, column=0, padx=10, pady=5, sticky="e")
        self.M_y_top_entry = ttk.Entry(self.lengthloads_frame)
        self.M_y_top_entry.grid(row=4, column=1, padx=10, pady=5, sticky="ew")
        self.M_y_top_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_y_top_units.grid(row=4, column=2, padx=10, pady=5, sticky="w")
        self.M_y_top_entry.bind("<Return>", self.update_variables_and_plot)
        self.M_y_top_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.M_y_bottom_label = ttk.Label(self.lengthloads_frame, text="M_y [bottom]:")
        self.M_y_bottom_label.grid(row=5, column=0, padx=10, pady=5, sticky="e")
        self.M_y_bottom_entry = ttk.Entry(self.lengthloads_frame)
        self.M_y_bottom_entry.grid(row=5, column=1, padx=10, pady=5, sticky="ew")
        self.M_y_bottom_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_y_bottom_units.grid(row=5, column=2, padx=10, pady=5, sticky="w")
        self.M_y_bottom_entry.bind("<Return>", self.update_variables_and_plot)
        self.M_y_bottom_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.M_z_top_label = ttk.Label(self.lengthloads_frame, text="M_z [top]:")
        self.M_z_top_label.grid(row=4, column=3, padx=10, pady=5, sticky="e")
        self.M_z_top_entry = ttk.Entry(self.lengthloads_frame)
        self.M_z_top_entry.grid(row=4, column=4, padx=10, pady=5, sticky="ew")
        self.M_z_top_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_z_top_units.grid(row=4, column=5, padx=10, pady=5, sticky="w")
        self.M_z_top_entry.bind("<Return>", self.update_variables_and_plot)
        self.M_z_top_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.M_z_bottom_label = ttk.Label(self.lengthloads_frame, text="M_z [botom]:")
        self.M_z_bottom_label.grid(row=5, column=3, padx=10, pady=5, sticky="e")
        self.M_z_bottom_entry = ttk.Entry(self.lengthloads_frame)
        self.M_z_bottom_entry.grid(row=5, column=4, padx=10, pady=5, sticky="ew")
        self.M_z_bottom_units = ttk.Label(self.lengthloads_frame, text="kNm")
        self.M_z_bottom_units.grid(row=5, column=5, padx=10, pady=5, sticky="w")
        self.M_z_bottom_entry.bind("<Return>", self.update_variables_and_plot)
        self.M_z_bottom_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.plot_results_button = ttk.Button(self.root, text="Plot Failure Envelopes", command = self.plot_envelopes)
        self.plot_results_button.grid(row=5, column=6, padx=10, pady=5, sticky="ew")

        self.plot_frame = ttk.Frame(self.root)
        self.plot_frame.grid(row=6, column=6, columnspan=4, rowspan=10, padx=10, pady=10, sticky="nsew")

        # Add first plot frame
        self.plot1_frame = ttk.Frame(self.plot_frame)
        self.plot1_frame.grid(row=0, column=0, padx=5, pady=5)
        self.fig1, self.ax1 = plt.subplots(figsize=(6, 6))
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.plot1_frame)
        self.canvas1.get_tk_widget().pack()
        self.ax1.grid(True)
        self.ax1.set_title("Major Axis interaction Envelope")
        self.ax1.set_xlabel("M [kNm]")
        self.ax1.set_ylabel("N [kN]")

        # Add second plot frame next to the first
        self.plot2_frame = ttk.Frame(self.plot_frame)
        self.plot2_frame.grid(row=0, column=2, padx=5, pady=5)
        self.fig2, self.ax2 = plt.subplots(figsize=(6, 6))
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.plot2_frame)
        self.canvas2.get_tk_widget().pack()
        self.ax2.grid(True)
        self.ax2.set_title("Minor Axis Interaction Envelope")
        self.ax2.set_xlabel("M [kNm]")
        self.ax2.set_ylabel("N [kN]")

        # Radio buttons for shape selection
        self.radio_frame = ttk.Frame(self.root)
        self.radio_frame.grid(row=5, column=1, padx=10, pady=10)
        
        shapes = ["rectangular", "circular", "arbitrary"]
        for shape in shapes:
            ttk.Radiobutton(self.radio_frame, text=shape.capitalize(), value=shape, variable=self.shape_var, command=self.update_input_fields).pack(side="left", padx=5)

        # Bind the "f_yk" input field
        self.f_yk_entry.bind("<Return>", self.update_variables_and_plot)
        self.f_yk_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "E_s" input field
        self.E_s_entry.bind("<Return>", self.update_variables_and_plot)
        self.E_s_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "bar diameter" entry field
        self.bar_dia_entry.bind("<Return>", self.update_variables_and_plot)
        self.bar_dia_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "link diameter" entry field
        self.link_dia_entry.bind("<Return>", self.update_variables_and_plot)
        self.link_dia_entry.bind("<FocusOut>", self.update_variables_and_plot)
        
        # Create frame for dynamic input fields
        self.input_frame = ttk.Frame(self.root)
        self.input_frame.grid(row=6, column=0, columnspan = 3, padx=5,  pady=10)

        # Initialize input fields for the default selection (rectangular)
        self.update_input_fields()

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
        self.fig0, self.ax0 = plt.subplots(figsize=(3,3))
        self.canvas0 = FigureCanvasTkAgg(self.fig0, master=self.input_frame)
        self.canvas0.get_tk_widget().grid(row=5, column=0, columnspan=4, padx=5, pady=5)
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)
        #self.ax0.set_axis_on()
        self.ax0.grid(True)

        # Rectangular input fields
        ttk.Label(self.input_frame, text="Section width, b:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.b_entry = ttk.Entry(self.input_frame)
        self.b_entry.insert(0, "500")
        self.b_entry.grid(row=0, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="mm").grid(row=0, column=3, padx=5, pady=5, sticky="w")
        
        ttk.Label(self.input_frame, text="Section height, h:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.h_entry = ttk.Entry(self.input_frame)
        self.h_entry.insert(0, "500")
        self.h_entry.grid(row=1, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="mm").grid(row=1, column=3, padx=5, pady=5, sticky="w")
        
        ttk.Label(self.input_frame, text="Cover:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.cover_entry = ttk.Entry(self.input_frame)
        self.cover_entry.insert(0, "25")
        self.cover_entry.grid(row=2, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="mm").grid(row=2, column=3, padx=5, pady=5, sticky="w")

        ttk.Label(self.input_frame, text="Number of bars distributed along the width of the column\nn_x:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.n_x_entry = ttk.Entry(self.input_frame)
        self.n_x_entry.grid(row=3, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="#").grid(row=3, column=3, padx=5, pady=5, sticky="w")

        ttk.Label(self.input_frame, text="Number of bars distributed along the height of the column\n n_y:").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.n_y_entry = ttk.Entry(self.input_frame)
        self.n_y_entry.grid(row=4, column=2, padx=5, pady=5)                           
        
        ## BIND FUNCTIONS TO UPDATE VARIABLES AND PLOT WHEN VALUES ARE CHANGED
        # Bind the "b" entry field
        self.b_entry.bind("<Return>", self.update_variables_and_plot)
        self.b_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "h" entry field
        self.h_entry.bind("<Return>", self.update_variables_and_plot)
        self.h_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "cover" entry field
        self.cover_entry.bind("<Return>", self.update_variables_and_plot)
        self.cover_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "n_x" entry field
        self.n_x_entry.bind("<Return>", self.update_variables_and_plot)
        self.n_x_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "n_y" entry field
        self.n_y_entry.bind("<Return>", self.update_variables_and_plot)
        self.n_y_entry.bind("<FocusOut>", self.update_variables_and_plot)
        ttk.Label(self.input_frame, text="#").grid(row=4, column=3, padx=5, pady=5, sticky="w")

        self.canvas0.draw()

    def create_circular_inputs(self):
        # Circular input fields
        ttk.Label(self.input_frame, text="Diameter, d:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.diameter_entry = ttk.Entry(self.input_frame)
        self.diameter_entry.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(self.input_frame, text="Number of bars:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.radial_num_bars_entry = ttk.Entry(self.input_frame)
        self.radial_num_bars_entry.grid(row=1, column=1, padx=5, pady=5)
        
        ttk.Label(self.input_frame, text="Cover:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.cover_entry = ttk.Entry(self.input_frame)
        self.cover_entry.grid(row=2, column=1, padx=5, pady=5)
        
        # Bind the "diameter" entry field
        self.diameter_entry.bind("<Return>", self.update_variables_and_plot)
        self.diameter_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "radial_num_bars" entry field
        self.radial_num_bars_entry.bind("<Return>", self.update_variables_and_plot)
        self.radial_num_bars_entry.bind("<FocusOut>", self.update_variables_and_plot)
        # Bind the "cover" entry field
        self.cover_entry.bind("<Return>", self.update_variables_and_plot)
        self.cover_entry.bind("<FocusOut>", self.update_variables_and_plot)

        self.fig0, self.ax0 = plt.subplots(figsize=(3,3))
        self.canvas0 = FigureCanvasTkAgg(self.fig0, master=self.input_frame)
        self.canvas0.get_tk_widget().grid(row=3, column=0, columnspan=4, padx=5, pady=5)
        self.fig0.subplots_adjust(left=0.125, right=0.9, top=0.95, bottom=0.1)
        #self.ax0.set_axis_on()
        self.ax0.grid(True)
        self.canvas0.draw()

    def create_arbitrary_inputs(self):
        # Arbitrary input fields
        ttk.Label(self.input_frame, text="Input vertices for section boundary [mm]:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.vertex_x_input = ttk.Entry(self.input_frame, width=10)
        self.vertex_x_input.grid(row=0, column=1, padx=5, pady=5)
        self.vertex_y_input = ttk.Entry(self.input_frame, width=10)
        self.vertex_y_input.grid(row=0, column=2, padx=5, pady=5)
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
            x = float(self.vertex_x_input.get())
            y = float(self.vertex_y_input.get())
            self.coord_list.append([x, y])
            
            # Update the list of coordiantes
            self.update_coord_display()
            # Update the displayed geometry
            self.update_section_plot()
            
            # Clear the input fields
            self.vertex_x_input.delete(0, tk.END)
            self.vertex_y_input.delete(0, tk.END)
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
        print("rebar_coordinates:", self.bar_list)
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

    def plot_envelopes (self):
        concrete_properties = get_concrete_properties(self.concrete_dropdown.get())
        concrete_section = Concrete_Section(shape=self.shape, h=self.h, b=self.b, diameter = self.diameter, vertices = self.vertices)
        reinforcement = Reinforcement(self.f_yk, self.E_s, self.bar_dia, self.link_dia, shape=self.shape, cover = self.cover, b = self.b, h = self.h, diameter = self.diameter, num_of_rows_of_rebar = self.n_y, num_of_cols_of_rebar = self.n_x, radial_number = self.radial_num_bars, bar_list = self.bar_list, gamma_s=1.15)
        print(f'L_effy: {self.L_effy}, L_effz: {self.L_effz}')
        print(f'concrete_section: {concrete_section}')
        print(f'reinforcement: {reinforcement}')
        column = Column(concrete_section, reinforcement, concrete_properties, self.L_effy, self.L_effz)
        print(f'moments collected: M_y_top: {self.M_y_top}, M_y_bottom: {self.M_y_bottom}, M_z_top: {self.M_z_top}, M_z_bottom: {self.M_z_bottom}')
        My_01, My_02, Mz_01, Mz_02, M_Edy, M_Edz, slenderness_y, slenderness_z, slenderness_ratio_y, slenderness_ratio_z = check_slenderness(column, self.N_Ed, self.M_y_top, self.M_y_bottom, self.M_z_top, self.M_z_bottom)
        plot_major_axis_failure_envelope(self.ax1, self.canvas1, column, self.N_Ed, M_Edy, My_02)
        plot_minor_axis_failure_envelope(self.ax2, self.canvas2, column, self.N_Ed, M_Edz, Mz_02)
