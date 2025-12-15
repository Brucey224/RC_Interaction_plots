

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
    slenderness_ratio_y = column.L_eff_y*1e3 / column.concrete_section.i_y
    slenderness_ratio_z = column.L_eff_z*1e3 / column.concrete_section.i_z

    if column.concrete_section.shape == "rectangular":
        e_0 = max(column.concrete_section.h/30, 20)*1e-3 # eccentricity dimension in m 
    elif column.concrete_section.shape == "circular":
        e_0 = max(column.concrete_section.diameter/30, 20)*1e-3 # eccentricity dimension in m
    elif column.concrete_section.shape == "arbitrary":
        e_0  = max((column.concrete_section.top_of_section - column.concrete_section.bottom_of_section)/30, 20)*1e-3 # eccentricity dimension in m
    # Check slenderness limits
    print(f"reinforcement area: {column.reinforcement.A_s}")
    mechanical_reinforcement_ratio = (column.reinforcement.A_s*column.reinforcement.f_yd) / (column.concrete_section.A*column.concrete_properties.f_cd)
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
    n = N_Ed / (column.concrete_section.A * column.concrete_properties.f_cd * 1e-3)  # axial force utilisation
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
    mechanical_reinforcement_ratio = (column.reinforcement.A_s*column.reinforcement.f_yd) / (column.concrete_section.A*column.concrete_properties.f_cd) 
    n_u = 1 + mechanical_reinforcement_ratio
    n_bal = 0.4
    K_r = min((n_u - n)/(n_u - n_bal), 1)
    print(f'K_r: {K_r}')
    if column.concrete_section.shape == "rectangular":
        e_0 = max(column.concrete_section.h/30, 20)*1e-3 # eccentricity dimension in m 
    elif column.concrete_section.shape == "circular":
        e_0 = max(column.concrete_section.diameter/30, 20)*1e-3 # eccentricity dimension in m
    elif column.concrete_section.shape == "arbitrary":
        e_0  = max((column.concrete_section.top_of_section - column.concrete_section.bottom_of_section)/30, 20)*1e-3 # eccentricity dimension in m

    beta = 0.35 + column.concrete_properties.f_ck/200 - slenderness_ratio/150
    K_phi = max (1, (1+beta*phi_eff))
    e_2y = (0.1* (K_r * K_phi * column.reinforcement.f_yd) / (0.45 * column.dy * column.reinforcement.E_s*1e3) * (column.L_eff_y*1e3)**2) *1e-3 # eccentricity dimension in m
    M_0e_y = 0.6*M_02y + 0.4*M_01y
    M_2y = N_Ed * e_2y
    M_Ed_y = max(M_02y, (M_0e_y + M_2y), (M_01y + 0.2*M_2y), N_Ed*e_0)
    
    return M_Ed_y

def compute_minor_axis_slender_moments(column, slenderness_ratio, N_Ed, n, M_01z, M_02z, phi_eff = 2.15):
    
    mechanical_reinforcement_ratio = (column.reinforcement.A_s*column.reinforcement.f_yd) / (column.concrete_section.A*column.concrete_properties.f_cd)
    n_u = 1 + mechanical_reinforcement_ratio
    n_bal = 0.4
    if (n_u - n)/(n_u - n_bal) >=1:
        K_r = 1.0
    else:
        K_r = (n_u - n)/(n_u - n_bal)
    
    if column.concrete_section.shape == "rectangular":
        e_0 = max(column.concrete_section.h/30, 20)*1e-3 # eccentricity dimension in m 
    elif column.concrete_section.shape == "circular":
        e_0 = max(column.concrete_section.diameter/30, 20)*1e-3 # eccentricity dimension in m
    elif column.concrete_section.shape == "arbitrary":
        e_0  = max((column.concrete_section.top_of_section - column.concrete_section.bottom_of_section)/30, 20)*1e-3 # eccentricity dimension in m

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
    h_0 = 2*column.concrete_section.A / column.concrete_section.polygon.length
    phi_0 = phi_RH * beta_fcm * beta_t_0
    phi_inf_t0 = phi_RH * beta_fcm * beta_t_0  
    phi_effy = phi_inf_t0 * M_0Eqp_y / M_Ed_y
    phi_effz = phi_inf_t0 * M_0Eqp_z / M_Ed_z

    return phi_effy, phi_effz