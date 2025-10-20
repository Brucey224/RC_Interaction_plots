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




