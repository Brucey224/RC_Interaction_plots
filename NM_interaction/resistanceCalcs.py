from shapely.geometry import Polygon
import math
from .geometryUtils import compute_area_and_centroid, clip_polygon_at_y, clip_polygon_at_x, integrate_part_of_circle

def compute_plastic_axial_capacity(column):
    Npl_Rd = (column.concrete_section.A * column.concrete_properties.f_cd + len(column.reinforcement.arrangement) * math.pi * column.reinforcement.bar_diameter**2 / 4 * column.concrete_properties.eps_c3 * column.reinforcement.E_s*1e3)*1e-3 #kN
    T_pl_Rd = - (len(column.reinforcement.arrangement) * math.pi * column.reinforcement.bar_diameter**2 / 4 * column.reinforcement.f_yd)*1e-3 #kN
    return T_pl_Rd, Npl_Rd

def determine_generic_major_axis_envelope_value(column, lambd, neutral_axis, moment_direction = 'positive'):
    N_Rd = 0 # Axial capacity in kN
    M_Rd = 0 # Moment capacity in kNm
    if moment_direction == 'positive':
        curvature = 1
    elif moment_direction == 'negative':
        curvature = -1 
    
    if column.concrete_section.shape == 'rectangular':
        section_centroid = [column.concrete_section.b/2, column.concrete_section.h/2]
        d_cy = min(lambd*neutral_axis, column.concrete_section.h)
        if neutral_axis <= column.concrete_section.h:
            in_section = True
            lever_arm_concrete_y = neutral_axis * (1 - lambd/2)
            F_cc = d_cy * column.concrete_section.b * column.concrete_properties.f_cd *1e-3

    
def determine_envelope_value_major_axis_positive(column, lambd, neutral_axis_y):       # Determines envelope value for major axis bending
    N_Rd = 0 # Axial capacity in kN
    M_Rdy = 0 # Moment capacity in kNm
    steel_contribution_M = 0
    steel_contribution_N = 0
    if column.concrete_section.shape == 'rectangular':
        section_centroid = [column.concrete_section.b/2, column.concrete_section.h/2]
        d_cy = min(lambd*neutral_axis_y, column.concrete_section.h)
        if neutral_axis_y <= column.concrete_section.h:
            in_section = True
            lever_arm_concrete_y = neutral_axis_y * (1 - lambd/2)
            F_cc = d_cy * column.concrete_section.b * column.concrete_properties.f_cd *1e-3
            N_Rd += F_cc
            M_Rdy += F_cc * lever_arm_concrete_y *1e-3
            concrete_contribution_N = F_cc
            concrete_contribution_M = F_cc * lever_arm_concrete_y *1e-3
        else:
            in_section = False
            lever_arm_concrete_y = (neutral_axis_y - d_cy/2)
            F_cc = d_cy * column.concrete_section.b * column.concrete_properties.f_cd *1e-3
            N_Rd += F_cc
            M_Rdy += F_cc * lever_arm_concrete_y * 1e-3
            concrete_contribution_N = F_cc
            concrete_contribution_M = F_cc * lever_arm_concrete_y *1e-3
    
    elif column.concrete_section.shape == 'circular':
        d_cy = min(lambd*neutral_axis_y, column.concrete_section.diameter)
        section_centroid = [column.concrete_section.diameter/2, column.concrete_section.diameter/2]
        if neutral_axis_y <= column.concrete_section.diameter:
            in_section = True
            segment_area, centroid_y, lever_arm_concrete_y = integrate_part_of_circle(neutral_axis_y, column.concrete_section.diameter, lambd) # area in mm^2, centroid distance in mm
            F_cc = segment_area * column.concrete_properties.f_cd * 1e-3
            N_Rd += F_cc
            concrete_contribution_N = F_cc
            concrete_contribution_M = F_cc*lever_arm_concrete_y *1e-3
            M_Rdy += F_cc*lever_arm_concrete_y *1e-3
        else:
            in_section = False
            segment_area, centroid_y, lever_arm_concrete_y = integrate_part_of_circle(neutral_axis_y, column.concrete_section.diameter, lambd) # area in mm^2, centroid distance in mm
            F_cc = segment_area*column.concrete_properties.f_cd * 1e-3
            N_Rd += F_cc
            concrete_contribution_N = F_cc
            concrete_contribution_M = F_cc*lever_arm_concrete_y *1e-3
            M_Rdy += F_cc*lever_arm_concrete_y *1e-3

        #Define concrete strain according to neutral axis location
    
    elif column.concrete_section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.concrete_section.polygon)[1]
        if neutral_axis_y <= column.concrete_section.h:
            in_section = True
        else:
            in_section = False
        
        clipping_plane = column.concrete_section.bottom_of_section + lambd*neutral_axis_y
        clipped_concrete = clip_polygon_at_y(column.concrete_section.polygon, clipping_plane)
        #cycle through the geometries created by clipping the arbitrary shape
        concrete_contribution_M = 0
        concrete_contribution_N = 0
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                # determine area and centroid of these shapes created
                area, centroid = compute_area_and_centroid(geom) # distance to centroid in mm, are in mm^2
                # determine whether the centroid is above the neutral axis (i.e. in compression)
                if centroid[1] < (clipping_plane):
                    concrete_contribution_N += area * column.concrete_properties.f_cd *1e-3
                    lever_arm_concrete_y = (neutral_axis_y - (centroid[1] - column.concrete_section.bottom_of_section))
                    concrete_contribution_M += concrete_contribution_N
                    N_Rd += area * column.concrete_properties.f_cd *1e-3
                    M_Rdy += lever_arm_concrete_y * area * column.concrete_properties.f_cd *1e-6
    
    steel_strains = []
    steel_stresses = []

    # determnine concrete strain based on whether neutral axis is within section or not
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2
    else:
        if column.concrete_section.shape == "rectangular" or column.concrete_section.shape == "arbitrary":
            concrete_strain = column.concrete_properties.eps_c3*neutral_axis_y/(neutral_axis_y - column.concrete_section.h/2)
        elif column.concrete_section.shape == "circular":
            concrete_strain = column.concrete_properties.eps_c3*neutral_axis_y/(neutral_axis_y - column.concrete_section.diameter/2)
    
    # loop through steel reinforcement, determine stresses and add to moment capacity / axial capacity
    for rebar_coords in column.reinforcement.arrangement:
        if in_section == True:
            if column.concrete_section.shape == "rectangular" or column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (neutral_axis_y - rebar_coords[1])/neutral_axis_y
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * ((neutral_axis_y) - (rebar_coords[1]- column.concrete_section.bottom_of_section)) / (neutral_axis_y)   # Steel strain computed if neutral axis is within the section. Sim triangles 
        else:
            if column.concrete_section.shape == "rectangular":
                steel_strain = concrete_strain * (neutral_axis_y - rebar_coords[1]) / (neutral_axis_y)
            elif column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (neutral_axis_y - rebar_coords[1]) / (neutral_axis_y)
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * (neutral_axis_y - (rebar_coords[1]- column.concrete_section.bottom_of_section)) / (neutral_axis_y)
                print(f'bottom of section: {column.concrete_section.bottom_of_section}, neutral axis: {neutral_axis_y}, rebar coords: {rebar_coords}, steel_strain: {steel_strain}')
        steel_stress = max(min(steel_strain*column.reinforcement.E_s*1e3, column.reinforcement.f_yd), -column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain)
        steel_stresses.append(steel_stress) # steel stress in MPa
        if column.concrete_section.shape == "rectangular":
            lever_arm_y = neutral_axis_y - rebar_coords[1]
        elif column.concrete_section.shape == "circular":
            lever_arm_y = neutral_axis_y - rebar_coords[1]
        elif column.concrete_section.shape == "arbitrary":
            lever_arm_y = neutral_axis_y - (rebar_coords[1] - column.concrete_section.bottom_of_section)# lever arm from neutral axis in mm
        steel_force = steel_stress * math.pi * (column.reinforcement.bar_diameter)**2 /4 # force in N by multiplying MPa by area in mm^2
        N_Rd += steel_force * 1e-3 #kN
        M_Rdy += steel_force * lever_arm_y *1e-6 #kNm
        steel_contribution_M += steel_stress * math.pi * (column.reinforcement.bar_diameter)**2 /4 * lever_arm_y *1e-6
        steel_contribution_N += steel_stress * math.pi * (column.reinforcement.bar_diameter)**2 /4 * 1e-3
    
    # compute moment contribution from axial force acting on section
    if column.concrete_section.shape == "rectangular":
        M_Rdy += N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3 #kNm
        axial_contribution_M = N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3
    elif column.concrete_section.shape == "circular":
        M_Rdy += N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3 #kNm
        axial_contribution_M = N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3
    elif column.concrete_section.shape == "arbitrary":
        M_Rdy += N_Rd*(section_centroid[1] - (neutral_axis_y - column.concrete_section.bottom_of_section)) * 1e-3 #kNm
        axial_contribution_M = N_Rd*(section_centroid[1] - neutral_axis_y) * 1e-3

    return N_Rd, M_Rdy, steel_stresses, steel_strains

def determine_envelope_value_major_axis_negative(column, lambd, neutral_axis_y):       # Determines envelope value for major axis bending
    N_Rd = 0 # Axial capacity in kN
    M_Rdy = 0 # Moment capacity in kNm
    
    # DEFINE NEUTRAL AXIS FROM OTHER DIRECTION IN ORDER TO GET NEGATIVE MOMENTS
    if column.concrete_section.shape == 'rectangular':
        section_centroid = [column.concrete_section.b/2, column.concrete_section.h/2]
        d_cy = min(lambd * neutral_axis_y, column.concrete_section.h)
        lever_arm_concrete_y = d_cy/2 - neutral_axis_y
        N_Rd += d_cy * column.concrete_section.b * column.concrete_properties.f_cd *1e-3
        if (column.concrete_section.h - neutral_axis_y) > 0:
            in_section = True
            M_Rdy += d_cy * column.concrete_section.b * column.concrete_properties.f_cd * lever_arm_concrete_y *1e-6
        else:
            in_section = False
            M_Rdy += d_cy * column.concrete_section.b * column.concrete_properties.f_cd * lever_arm_concrete_y * 1e-6
    
    elif column.concrete_section.shape == 'circular':
        section_centroid = [column.concrete_section.diameter/2, column.concrete_section.diameter/2]
        if (column.concrete_section.diameter - neutral_axis_y) > 0:
            in_section = True
            segment_area, centroid_y, lever_arm_y = integrate_part_of_circle(neutral_axis_y, column.concrete_section.diameter, lambd) # area in mm^2, centroid distance in mm
            F_cc = segment_area * column.concrete_properties.f_cd * 1e-3
            N_Rd += F_cc
            M_Rdy += - F_cc * lever_arm_y *1e-3
        else:
            in_section = False
            segment_area, centroid_y, lever_arm_concrete_y = integrate_part_of_circle(neutral_axis_y, column.concrete_section.diameter, lambd) # area in mm^2, centroid distance in mm
            F_cc = segment_area*column.concrete_properties.f_cd * 1e-3
            N_Rd += F_cc
            M_Rdy += - F_cc * lever_arm_concrete_y *1e-3
    
    elif column.concrete_section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.concrete_section.polygon)[1]
        if (column.concrete_section.top_of_section - neutral_axis_y) > column.concrete_section.bottom_of_section:
            in_section = True
        else:
            in_section = False
        
        clipping_plane = column.concrete_section.top_of_section - lambd*neutral_axis_y
        clipped_concrete = clip_polygon_at_y(column.concrete_section.polygon, clipping_plane)
        # cycle through the geometries created by clipping the arbitrary shape
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                # determine area and centroid of these shapes created
                area, centroid = compute_area_and_centroid(geom) # distance to centroid in mm, are in mm^2
                # determine whether the centroid is above the neutral axis (i.e. in compression)
                if centroid[1] > (clipping_plane):
                    N_Rd += area * column.concrete_properties.f_cd *1e-3
                    lever_arm_y = (column.concrete_section.top_of_section - neutral_axis_y - centroid[1])
                    M_Rdy += lever_arm_y * area * column.concrete_properties.f_cd *1e-6
    
    steel_strains = []
    steel_stresses = []

    # determnine concrete strain based on whether neutral axis is within section or not
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2
    else:
        if column.concrete_section.shape == "rectangular" or column.concrete_section.shape == "arbitrary":
            concrete_strain = column.concrete_properties.eps_c3*neutral_axis_y/(neutral_axis_y - column.concrete_section.h/2)
        elif column.concrete_section.shape == "circular":
            concrete_strain = column.concrete_properties.eps_c3*neutral_axis_y/(neutral_axis_y - column.concrete_section.diameter/2)
    
    # loop through steel reinforcement, determine stresses and add to moment capacity / axial capacity
    for rebar_coords in column.reinforcement.arrangement:
        if in_section == True:
            if column.concrete_section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.concrete_section.h - neutral_axis_y)) / neutral_axis_y
            elif column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.concrete_section.diameter - neutral_axis_y)) / neutral_axis_y
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.concrete_section.top_of_section - neutral_axis_y)) / neutral_axis_y
        else:
            if column.concrete_section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.concrete_section.h - neutral_axis_y)) / (neutral_axis_y)
            elif column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.concrete_section.diameter - neutral_axis_y)) / (neutral_axis_y)
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[1] - (column.concrete_section.top_of_section - neutral_axis_y)) / (neutral_axis_y)
        
        steel_stress = max(min(steel_strain*column.reinforcement.E_s*1e3, column.reinforcement.f_yd), - column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain) 
        steel_stresses.append(steel_stress) # steel stress in MPa

        if column.concrete_section.shape == "rectangular":
            lever_arm_y = (column.concrete_section.h - neutral_axis_y) - rebar_coords[1] 
        elif column.concrete_section.shape == "circular":
            lever_arm_y = (column.concrete_section.diameter - neutral_axis_y) - rebar_coords[1]
        elif column.concrete_section.shape == "arbitrary":
            lever_arm_y = (column.concrete_section.top_of_section - neutral_axis_y) - rebar_coords[1] # lever arm from neutral axis in mm
        
        N_Rd += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * 1e-3 #kN
        M_Rdy += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * lever_arm_y *1e-6 #kNm

    # compute moment contribution from axial force acting on section
    if column.concrete_section.shape == "rectangular":
        M_Rdy += N_Rd*(section_centroid[1] - (column.concrete_section.h - neutral_axis_y)) * 1e-3
    elif column.concrete_section.shape == "circular":
        M_Rdy += N_Rd*(section_centroid[1] - (column.concrete_section.diameter - neutral_axis_y)) * 1e-3        
    elif column.concrete_section.shape == "arbitrary":
        M_Rdy += N_Rd*(section_centroid[1] - (column.concrete_section.top_of_section - neutral_axis_y)) * 1e-3
    
    return N_Rd, M_Rdy, steel_stresses, steel_strains

def determine_envelope_value_minor_axis_positive(column, lambd, neutral_axis_x):       # Determines envelope value for major axis bending
    N_Rd = 0 #Axial capacity in kN
    M_Rdz = 0 #Moment capacity in kNm

    if column.concrete_section.shape == 'rectangular':
        section_centroid = [column.concrete_section.b/2, column.concrete_section.h/2]
        d_cx = min(lambd*neutral_axis_x, column.concrete_section.b)
        lever_arm_concrete_x = neutral_axis_x  - d_cx/2                                                # lever arm in mm
        N_Rd += d_cx * column.concrete_section.h * column.concrete_properties.f_cd *1e-3
        if neutral_axis_x <= column.concrete_section.b:
            in_section = True                              # Add axial force from concrete in compression
            M_Rdz += d_cx * column.concrete_section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6      # Add moment from concrete in compression acting at a lever arm from the neutral axis
        else:
            in_section = False                                           # lever arm in mm
            M_Rdz += d_cx * column.concrete_section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
    
    elif column.concrete_section.shape == 'circular':
        section_centroid = [column.concrete_section.diameter/2, column.concrete_section.diameter/2]
        if neutral_axis_x <= column.concrete_section.diameter:
            in_section = True
            segment_area, centroid_x, lever_arm_concrete_x = integrate_part_of_circle(neutral_axis_x, column.concrete_section.diameter, lambd) # area in mm^2, centroid distance in mm
            N_Rd += segment_area*column.concrete_properties.f_cd *1e-3
            M_Rdz += segment_area*column.concrete_properties.f_cd*lever_arm_concrete_x *1e-6
        else:
            in_section = False
            lever_arm_concrete_x = neutral_axis_x - section_centroid[0]
            N_Rd += column.concrete_section.A*column.concrete_properties.f_cd * 1e-3
            M_Rdz += column.concrete_section.A*column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
        
    
    elif column.concrete_section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.concrete_section.polygon)[1]                                 # get centroid of polygon of the whole section      
        if column.concrete_section.left_of_section + neutral_axis_x < column.concrete_section.right_of_section:                   # Determine whether neutral axis is within the section or not
            in_section = True
        else:
            in_section = False

        clipping_plane = column.concrete_section.left_of_section + lambd*neutral_axis_x
        clipped_concrete = clip_polygon_at_x(column.concrete_section.polygon, clipping_plane) #Clip the polygon into multiple polygons
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                area, centroid = compute_area_and_centroid(geom)                                                # area in mm^2, centroid in mm from extreme compression fibre
                if centroid[0] < clipping_plane:                               # Find all the geometries created that are on the side of the NA under compression
                    N_Rd += area * column.concrete_properties.f_cd *1e-3                                        # Add axial force for concrete in compression 
                    M_Rdz += (column.concrete_section.left_of_section + neutral_axis_x - centroid[0]) * area * column.concrete_properties.f_cd * 1e-6     # Moment resistance contribution from concrete in compression acting at a lever arm from the neutral axis 
    
    steel_strains = []
    steel_stresses = []
    
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2                                                    # Concrete strain if neutral axis within section
    else:
        concrete_strain = column.concrete_properties.eps_c3*neutral_axis_x/(neutral_axis_x - section_centroid[1])                                                     # Concrete strain if neutral axis outside of the section. Note that this is the concrete strain at the centre of the section, therefore computing steel strains realtive to this should use sim triangle from centre of section
    
    for rebar_coords in column.reinforcement.arrangement:                                                       
        if in_section ==True:
            if column.concrete_section.shape == "rectangular" or column.concrete_section.shape == "circular":
                steel_strain = (neutral_axis_x - rebar_coords[0])/neutral_axis_x * concrete_strain
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * ((column.concrete_section.left_of_section + neutral_axis_x) - rebar_coords[0]) / (neutral_axis_x)   # Steel strain computed if neutral axis is within the section. Sim triangles 
        else:
            if column.concrete_section.shape == "rectangular":
                steel_strain = concrete_strain * (neutral_axis_x - rebar_coords[0]) / neutral_axis_x
            elif column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (neutral_axis_x - rebar_coords[0]) / neutral_axis_x
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * ((column.concrete_section.left_of_section + neutral_axis_x) - rebar_coords[0]) / (neutral_axis_x)
        
        steel_stress = max(min(steel_strain*column.reinforcement.E_s*1e3, column.reinforcement.f_yd), - column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain)
        steel_stresses.append(steel_stress) # steel stress in MPa
        if column.concrete_section.shape == "rectangular":
            lever_arm_x = neutral_axis_x - rebar_coords[0]
        elif column.concrete_section.shape == "circular":
            lever_arm_x = neutral_axis_x - rebar_coords[0]
        elif column.concrete_section.shape == "arbitrary":
            lever_arm_x = column.concrete_section.left_of_section + neutral_axis_x - rebar_coords[0] # lever arm in mm
        
        N_Rd += steel_stress * math.pi*(column.reinforcement.bar_diameter**2)/4 *1e-3
        M_Rdz += steel_stress * math.pi*(column.reinforcement.bar_diameter**2)/4 * lever_arm_x * 1e-6
    
    # compute moment contribution from axial force acting on section
    if column.concrete_section.shape == "rectangular":
        M_Rdz += N_Rd*(section_centroid[0] - neutral_axis_x) * 1e-3
    elif column.concrete_section.shape == "circular":
        M_Rdz += N_Rd*(section_centroid[0] - neutral_axis_x) * 1e-3
    elif column.concrete_section.shape == "arbitrary":
        M_Rdz += N_Rd*(section_centroid[0] - (column.concrete_section.left_of_section + neutral_axis_x)) * 1e-3

    return N_Rd, M_Rdz, steel_stresses, steel_strains

def determine_envelope_value_minor_axis_negative(column, lambd, neutral_axis_x):        # Determines envelope value for minor axis bending for positive moments (left half of section in compression)
    N_Rd = 0                                                                            #Initialise axial capacity in kN, value accumulates as cycling through concrete and rebar
    M_Rdz = 0             

    if column.concrete_section.shape == 'rectangular':
        section_centroid = [column.concrete_section.b/2, column.concrete_section.h/2]
        d_cx = min(lambd*neutral_axis_x, column.concrete_section.b)
        lever_arm_concrete_x = d_cx/2 - neutral_axis_x #lever arm in mm
        N_Rd += d_cx * column.concrete_section.h * column.concrete_properties.f_cd *1e-3
        if (column.concrete_section.b - neutral_axis_x) > 0:
            in_section = True
            M_Rdz += d_cx * column.concrete_section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
        else:
            in_section = False
            M_Rdz += d_cx * column.concrete_section.h * column.concrete_properties.f_cd * lever_arm_concrete_x *1e-6
    
    elif column.concrete_section.shape == 'circular':
        section_centroid = [column.concrete_section.diameter/2, column.concrete_section.diameter/2]
        if (column.concrete_section.diameter - neutral_axis_x) > 0:
            in_section = True
            segment_area, centroid_x, lever_arm_x = integrate_part_of_circle(neutral_axis_x, column.concrete_section.diameter, lambd) # area in mm^2, centroid distance in mm
            N_Rd += segment_area* column.concrete_properties.f_cd *1e-3
            M_Rdz += - lever_arm_x * segment_area * column.concrete_properties.f_cd *1e-6
        else:
            in_section = False
            N_Rd += column.concrete_section.A*column.concrete_properties.f_cd * 1e-3
            M_Rdz += column.concrete_section.A*column.concrete_properties.f_cd * (section_centroid[0] - neutral_axis_x) *1e-6
    
    elif column.concrete_section.shape == 'arbitrary':
        section_centroid = compute_area_and_centroid(column.concrete_section.polygon)[1]
        if (column.concrete_section.right_of_section - neutral_axis_x) > column.concrete_section.left_of_section:
            in_section = True
        else:
            in_section = False

        clipping_plane = column.concrete_section.right_of_section - lambd*neutral_axis_x
        clipped_concrete = clip_polygon_at_x(column.concrete_section.polygon, clipping_plane) 
        for geom in clipped_concrete:
            if isinstance(geom, Polygon):
                area, centroid = compute_area_and_centroid(geom) #area in mm^2, centroid in mm from extreme compression fibre
                if centroid[0] > clipping_plane:
                    N_Rd += area * column.concrete_properties.f_cd *1e-3
                    lever_arm_x = (column.concrete_section.right_of_section - neutral_axis_x - centroid[0])
                    M_Rdz +=  lever_arm_x * area * column.concrete_properties.f_cd * 1e-6
    
    steel_strains = []
    steel_stresses = []
    
    # DETERMINE CONCRETE STRAINS BASED ON WHETHER NEUTRAL AXIS IS WITHIN SECTION OR NOT
    if in_section == True:
        concrete_strain = column.concrete_properties.eps_cu2
    else:
        if column.concrete_section.shape == "rectangular" or column.concrete_section.shape == "arbitrary":
            concrete_strain = column.concrete_properties.eps_c3*neutral_axis_x/(neutral_axis_x - column.concrete_section.b/2)
        elif column.concrete_section.shape == "circular":
            concrete_strain = column.concrete_properties.eps_c3*neutral_axis_x/(neutral_axis_x - column.concrete_section.diameter/2)
    
    for rebar_coords in column.reinforcement.arrangement:
        if in_section == True:
            if column.concrete_section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.concrete_section.b - neutral_axis_x)) / neutral_axis_x
            elif column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.concrete_section.diameter - neutral_axis_x)) / neutral_axis_x
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.concrete_section.right_of_section - neutral_axis_x))/ neutral_axis_x
        else:
            if column.concrete_section.shape == "rectangular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.concrete_section.b - neutral_axis_x)) / (neutral_axis_x)
            elif column.concrete_section.shape == "circular":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.concrete_section.diameter - neutral_axis_x)) / (neutral_axis_x)
            elif column.concrete_section.shape == "arbitrary":
                steel_strain = concrete_strain * (rebar_coords[0] - (column.concrete_section.right_of_section - neutral_axis_x)) / (neutral_axis_x)

        steel_stress = max(min(steel_strain * column.reinforcement.E_s*1e3, column.reinforcement.f_yd), -column.reinforcement.f_yd) # Steel stress in MPa, limited to design yield strength of steel bars
        steel_strains.append(steel_strain)
        steel_stresses.append(steel_stress) # steel stress in MPa

        if column.concrete_section.shape == "rectangular":
            lever_arm_x = (column.concrete_section.b - neutral_axis_x) - rebar_coords[0]
        elif column.concrete_section.shape == "circular":
            lever_arm_x = (column.concrete_section.diameter - neutral_axis_x) - rebar_coords[0] # lever arm in mm
        elif column.concrete_section.shape == "arbitrary":
            lever_arm_x = (column.concrete_section.right_of_section - neutral_axis_x) - rebar_coords[0] # lever arm in mm

        N_Rd += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * 1e-3
        M_Rdz += steel_stress * math.pi/4*(column.reinforcement.bar_diameter**2) * lever_arm_x * 1e-6
        
    # compute moment contribution from axial force acting on section
    if column.concrete_section.shape == "rectangular":
        M_Rdz += N_Rd*(section_centroid[0] - (column.concrete_section.b - neutral_axis_x)) * 1e-3
    elif column.concrete_section.shape == "circular":
        M_Rdz += N_Rd*(section_centroid[0] - (column.concrete_section.diameter - neutral_axis_x)) * 1e-3
    elif column.concrete_section.shape == "arbitrary":
        M_Rdz += N_Rd*(section_centroid[0] - (column.concrete_section.right_of_section - neutral_axis_x)) * 1e-3
    
    return N_Rd, M_Rdz, steel_stresses, steel_strains

def compute_UTR(column, N_Ed, M_Edy, M_Edz):
    pass

