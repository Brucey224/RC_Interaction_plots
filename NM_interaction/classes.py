import math
from shapely.geometry import Polygon
from .geometryUtils import compute_second_moment_area

class Concrete_Section():
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
    def __init__(self, concrete_section, reinforcement, concrete_properties, L_eff_y, L_eff_z):
        self.concrete_section = concrete_section
        self.reinforcement = reinforcement
        self.concrete_properties = concrete_properties
        self.L_eff_y = float(L_eff_y)#*m
        self.L_eff_z = float(L_eff_z)#*m
        if concrete_section.shape == 'rectangular':
            self.dy = (concrete_section.h - reinforcement.cover - reinforcement.link_diameter - reinforcement.bar_diameter/2)#*mm
            self.d_2 = (reinforcement.cover + reinforcement.link_diameter + reinforcement.bar_diameter/2)#*mm
            self.dz = (concrete_section.b - reinforcement.cover - reinforcement.link_diameter - reinforcement.bar_diameter/2)#*mm
        elif concrete_section.shape == 'circular':
            self.dy = (concrete_section.diameter - reinforcement.cover - reinforcement.link_diameter - reinforcement.bar_diameter/2)#*mm
            self.dz = self.dy
            self.d_2 = (reinforcement.cover + reinforcement.link_diameter + reinforcement.bar_diameter/2)#*mm
        elif concrete_section.shape == 'arbitrary':
            self.dy = self.concrete_section.h - min(coord[1] for coord in self.reinforcement.arrangement)
            self.dy_2 = min(coord[1] for coord in self.reinforcement.arrangement)
            self.dz = self.concrete_section.b - min(coord[0] for coord in self.reinforcement.arrangement)
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
    def __init__(self, f_yk, E_s, bar_diameter, link_diameter, cover = None, shape=None, b =None, h = None, diameter = None, num_of_rows_of_rebar = None, num_of_cols_of_rebar = None, radial_number = None, bar_list = None, gamma_s=1.15):
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
            bar_list = []
            width_x = b - 2*cover - 2*link_diameter - bar_diameter
            height_y = h - 2*cover - 2*link_diameter - bar_diameter
            spacing_x = width_x / (num_of_cols_of_rebar-1)
            spacing_y = height_y / (num_of_rows_of_rebar - 1)
            #create rebar in top + bottom layer
            for i in range(num_of_cols_of_rebar):
                x1 = (cover + link_diameter + bar_diameter / 2 + i*spacing_x)#*mm
                y1 = (cover + link_diameter + bar_diameter/2)#*mm
                bar_list.append([x1, y1])
                x2 = x1
                y2 = h - cover - link_diameter - bar_diameter/2
                bar_list.append([x2,y2])
                A_s += 2* math.pi * (bar_diameter/2)**2
            # create side bars
            for j in range(num_of_rows_of_rebar-2):
                x1 = cover + link_diameter + bar_diameter/2
                y1 = cover+link_diameter + bar_diameter/2 + spacing_y*(j+1)
                bar_list.append([x1, y1])
                x2 = b - cover - link_diameter - bar_diameter/2
                y2 = y1
                bar_list.append([x2,y2])
                A_s += 2* math.pi * (bar_diameter/2)**2

        elif self.shape == 'circular':
            bar_list = []
            r = diameter / 2 - cover - link_diameter - bar_diameter/2
            radial_spacing = 2*math.pi / radial_number
            for i in range(radial_number):
                x = r*math.cos(radial_spacing*i) + diameter/2
                y = r*math.sin(radial_spacing*i) + diameter/2
                bar_list.append([x,y])
                A_s += math.pi * (bar_diameter/2)**2

        elif self.shape == 'arbitrary':
            A_s = len(bar_list)*math.pi * (bar_diameter/2)**2

        self.arrangement = bar_list
        self.A_s = A_s
