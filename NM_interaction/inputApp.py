import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import math
from .utils import get_concrete_properties, collect_user_input
from .plotting import plot_rectangular_section, plot_circular_section, plot_arbitrary_section, plot_major_axis_failure_envelope, plot_minor_axis_failure_envelope
from .codeChecks import check_slenderness

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

        self.plot_results_button = ttk.Button(self.root, text="Plot Failure Envelopes", command = self.plot_envelopes)
        self.plot_results_button.grid(row=5, column=6, padx=10, pady=5, sticky="ew")

        self.plot_frame = ttk.Frame(self.root)
        self.plot_frame.grid(row=6, column=6, columnspan=4, rowspan=10, padx=10, pady=10, sticky="nsew")

        # Add first plot frame
        self.plot1_frame = ttk.Frame(self.plot_frame)
        self.plot1_frame.grid(row=0, column=0, padx=5, pady=5)
        self.fig1, self.ax1 = plt.subplots(figsize=(4, 4))
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.plot1_frame)
        self.canvas1.get_tk_widget().pack()
        self.ax1.set_title("Major Axis interaction Envelope")
        self.ax1.set_xlabel("M [kNm]")
        self.ax1.set_ylabel("N [kN]")

        # Add second plot frame next to the first
        self.plot2_frame = ttk.Frame(self.plot_frame)
        self.plot2_frame.grid(row=0, column=2, padx=5, pady=5)
        self.fig2, self.ax2 = plt.subplots(figsize=(4, 4))
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.plot2_frame)
        self.canvas2.get_tk_widget().pack()
        self.ax2.set_title("Minor Axis Interaction Envelope")
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
        self.b_input = ttk.Entry(self.input_frame)
        self.b_input.insert(0, "500")
        self.b_input.grid(row=0, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="mm").grid(row=0, column=3, padx=5, pady=5, sticky="w")
        
        ttk.Label(self.input_frame, text="Section height, h:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.h_input = ttk.Entry(self.input_frame)
        self.h_input.insert(0, "500")
        self.h_input.grid(row=1, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="mm").grid(row=1, column=3, padx=5, pady=5, sticky="w")
        
        ttk.Label(self.input_frame, text="Cover:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.cover_input = ttk.Entry(self.input_frame)
        self.cover_input.insert(0, "25")
        self.cover_input.grid(row=2, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="mm").grid(row=2, column=3, padx=5, pady=5, sticky="w")

        ttk.Label(self.input_frame, text="Number of bars distributed along the width of the column\nn_x:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.n_x_input = ttk.Entry(self.input_frame)
        self.n_x_input.grid(row=3, column=2, padx=5, pady=5)
        ttk.Label(self.input_frame, text="#").grid(row=3, column=3, padx=5, pady=5, sticky="w")

        ttk.Label(self.input_frame, text="Number of bars distributed along the height of the column\n n_y:").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.n_y_input = ttk.Entry(self.input_frame)
        self.n_y_input.grid(row=4, column=2, padx=5, pady=5)                           
        
        ## BIND FUNCTIONS TO UPDATE PLOT WHEN VALUES ARE CHANGED
        self.b_input.bind("<Return>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.b_input.bind("<FocusOut>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.h_input.bind("<Return>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.h_input.bind("<FocusOut>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.cover_input.bind("<Return>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.cover_input.bind("<FocusOut>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.n_x_input.bind("<Return>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.n_x_input.bind("<FocusOut>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.n_y_input.bind("<Return>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        self.n_y_input.bind("<FocusOut>", lambda event: plot_rectangular_section(self.ax0, self.canvas0, self.b_input, self.h_input, self.cover_input, self.link_dia_entry, self.bar_dia_entry, self.n_x_input, self.n_y_input))
        ttk.Label(self.input_frame, text="#").grid(row=4, column=3, padx=5, pady=5, sticky="w")

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

    def plot_envelopes (self):
        if self.shape_var.get() == "rectangular":
            N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column = collect_user_input(self.shape_var, 
                                                                                        N_Ed_entry = self.N_Ed_entry, M_y_top_entry = self.M_y_top_entry, M_y_bottom_entry = self.M_y_bottom_entry, M_z_top_entry = self.M_z_top_entry, M_z_bottom_entry = self.M_z_bottom_entry, 
                                                                                        concrete_grade_entry = self.concrete_dropdown, 
                                                                                        L_effy_entry = self.L_Effy_entry, L_effz_entry = self.L_effz_entry, 
                                                                                        f_yk_entry=self.f_yk_entry, E_s_entry=self.E_s_entry, 
                                                                                        bar_dia_entry = self.bar_dia_entry, link_dia_entry = self.link_dia_entry, cover_entry = self.cover_input, 
                                                                                        h_input = self.h_input, b_input = self.b_input, n_x_input = self.n_x_input, n_y_input = self.n_y_input)
        if self.shape_var.get() == "circular":
            N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column = collect_user_input(self.shape_var, 
                                                                                        N_Ed_entry = self.N_Ed_entry, M_y_top_entry = self.M_y_top_entry, M_y_bottom_entry = self.M_y_bottom_entry, M_z_top_entry = self.M_z_top_entry, M_z_bottom_entry = self.M_z_bottom_entry, 
                                                                                        concrete_grade_entry = self.concrete_dropdown, 
                                                                                        L_effy_entry = self.L_Effy_entry, L_effz_entry = self.L_effz_entry, 
                                                                                        f_yk_entry=self.f_yk_entry, E_s_entry=self.E_s_entry, 
                                                                                        bar_dia_entry = self.bar_dia_entry, link_dia_entry = self.link_dia_entry, cover_entry = self.cover_input, 
                                                                                        diameter_entry = self.diameter_entry, radial_num_bars = self.radial_num_bars_input)
        if self.shape_var.get() == "arbitrary":
            N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column = collect_user_input(self.shape_var, 
                                                                                        N_Ed_entry = self.N_Ed_entry, M_y_top_entry = self.M_y_top_entry, M_y_bottom_entry = self.M_y_bottom_entry, M_z_top_entry = self.M_z_top_entry, M_z_bottom_entry = self.M_z_bottom_entry, 
                                                                                        concrete_grade_entry = self.concrete_dropdown, 
                                                                                        L_effy_entry = self.L_Effy_entry, L_effz_entry = self.L_effz_entry, 
                                                                                        f_yk_entry=self.f_yk_entry, E_s_entry=self.E_s_entry, 
                                                                                        bar_dia_entry = self.bar_dia_entry, link_dia_entry = self.link_dia_entry, cover_entry = self.cover_input, 
                                                                                        coord_list = self.coord_list,bar_list =  self.bar_list)
        My_01, My_02, Mz_01, Mz_02, M_Edy, M_Edz, slenderness_y, slenderness_z, slenderness_ratio_y, slenderness_ratio_z = check_slenderness(column, N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom)
        plot_major_axis_failure_envelope(self.ax1, self.canvas1, column, N_Ed, M_Edy, My_02)
        plot_minor_axis_failure_envelope(self.ax2, self.canvas2, column, N_Ed, M_Edz, Mz_02)


