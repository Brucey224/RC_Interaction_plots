from .classes import Concrete_Material, Section, Column, Reinforcement
import pandas as pd
import importlib.resources as pkg_resources

def get_concrete_properties(concrete_grade, gamma_c = 1.5, alpha_cc=0.85):
    #extract material properties from database
    with pkg_resources.open_text("NM_interaction.data", "concrete_grade_DB.csv") as file:
        df = pd.read_csv(file)
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

def collect_user_input(shape_var, **kwargs):
    """
    Collect user inputs based on the selected shape (rectangular, circular, or arbitrary).
    Only relevant inputs for the selected shape are processed.
    """
    shape = shape_var.get()
    
    # Common inputs for all shapes
    N_Ed = float(kwargs["N_Ed_entry"].get())
    M_y_top = float(kwargs["M_y_top_entry"].get())
    M_y_bottom = float(kwargs["M_y_bottom_entry"].get())
    M_z_top = float(kwargs["M_z_top_entry"].get())
    M_z_bottom = float(kwargs["M_z_bottom_entry"].get())
    concrete_grade = kwargs["concrete_grade_entry"].get()
    concrete_properties = get_concrete_properties(concrete_grade)
    L_eff_y = float(kwargs["L_effy_entry"].get())
    L_eff_z = float(kwargs["L_effz_entry"].get())
    f_yk = float(kwargs["f_yk_entry"].get())
    E_s = float(kwargs["E_s_entry"].get())
    bar_dia = float(kwargs["bar_dia_entry"].get())
    link_dia = float(kwargs["link_dia_entry"].get())

    # Shape-specific inputs
    if shape == "rectangular":
        h = float(kwargs["h_input"].get())
        b = float(kwargs["b_input"].get())
        n_x = int(kwargs["n_x_input"].get())
        n_y = int(kwargs["n_y_input"].get())
        cover = float(kwargs["cover_entry"].get())
        section = Section(shape=shape, h=h, b=b)
        rebar = Reinforcement(f_yk, E_s, bar_dia, link_dia, cover=cover, shape=shape, b=b, h=h, num_of_rows_of_rebar=n_y, num_of_cols_of_rebar=n_x)

    elif shape == "circular":
        diameter = float(kwargs["diameter_entry"].get())
        radial_number_bars = int(kwargs["radial_num_bars_input"].get())
        cover = float(kwargs["cover_entry"].get())
        section = Section(shape=shape, diameter=diameter)
        rebar = Reinforcement(f_yk, E_s, bar_dia, link_dia, cover=cover, shape=shape, diameter=diameter, radial_number=radial_number_bars)

    elif shape == "arbitrary":
        vertices = kwargs["coord_list"]
        rebar_coords = kwargs["bar_list"]
        section = Section(shape=shape, vertices=vertices)
        rebar = Reinforcement(f_yk, E_s, bar_dia, link_dia, shape=shape, rebar_coords=rebar_coords)

    else:
        raise ValueError(f"Unknown shape: {shape}")

    # Create the column object
    column = Column(section, rebar, concrete_properties, L_eff_y, L_eff_z)

    # Return the collected inputs
    return N_Ed, M_y_top, M_y_bottom, M_z_top, M_z_bottom, column



