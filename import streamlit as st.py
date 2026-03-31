import streamlit as st
import numpy as np
import math

# ==============================================================================
# WEB PAGE CONFIGURATION & HEADER
# ==============================================================================
# Added a water droplet icon to the browser tab
st.set_page_config(page_title="Air-Stripping | Osama Nimri", layout="wide", page_icon="💧")

# Create two columns: a wide one for the text (size 4), a narrow one for the logo (size 1)
col_text, col_logo = st.columns([4, 1])

with col_text:
    # 1. Professional Custom Header using HTML/CSS
    st.markdown("""
    <div style="padding-bottom: 10px;">
        <h1 style="font-size: 3.2rem; margin-bottom: 0px; padding-bottom: 0px;">Air-Stripping Calculator</h1>
        <h3 style="font-size: 1.4rem; margin-top: 0px; color: #666666; font-weight: 400;">
            Developed by <b>Osama Nimri</b>
        </h3>
        <hr style="border: none; height: 3px; background-color: #005eb8; width: 120px; margin-left: 0; margin-top: 15px;">
    </div>
    """, unsafe_allow_html=True)

    st.markdown("#### Multicomponent Packed Tower Aeration Design")
    
    # 2. Clean UI: Hiding the technical explanation inside a dropdown menu
    with st.expander("📖 Read Technical Methodology", expanded=False):
        st.markdown("""
        This reactive application computes strict hydrodynamic sizing, phase equilibria, and 
        mass transfer kinetics using the **Eckert Generalized Pressure Drop Correlation** and the 
        **Onda Mass Transfer Model**. It implements the multicomponent $Z_{max}$ bounding logic 
        to ensure regulatory compliance across all target profiles.
        """)

with col_logo:
    # Keeps your university logo in the top right
    st.image("R.png", width=150)

# ==============================================================================
# SIDEBAR: SYSTEM INPUT PARAMETERS
# ==============================================================================
st.sidebar.header("System Operational Parameters")



# Operating Conditions
T_C = st.sidebar.slider("Design Temp (°C)", 
                        min_value=5.0, max_value=30.0, value=15.0, step=0.5)
Q_water_m3day = st.sidebar.number_input("Well Yield Flow Rate (m³/day)", value=1000.0)
R_s = st.sidebar.number_input("Design Stripping Factor (Rs)", 
                              min_value=1.1, max_value=10.0, value=3.0, step=0.1)

st.sidebar.header("Contaminant Profiling")
# TCE Influent targets
st.sidebar.markdown("**Trichloroethylene (TCE)**")
C_in_TCE = st.sidebar.number_input("TCE Influent (µg/L)", value=200.0)
C_out_TCE = st.sidebar.number_input("TCE Effluent Target (µg/L)", value=5.0)

# PCE Influent targets
st.sidebar.markdown("**Tetrachloroethylene (PCE)**")
C_in_PCE = st.sidebar.number_input("PCE Influent (µg/L)", value=275.0)
C_out_PCE = st.sidebar.number_input("PCE Effluent Target (µg/L)", value=1.0)

st.sidebar.header("Packing Media Selection")

# FIX 1: Added choices to the selectbox
packing_type = st.sidebar.selectbox("Ceramic Topology", 
                                    ("1.5-inch Berl Saddles", "1-inch Berl Saddles"))

# Assign Empirical Packing Constants 
# (a_t: m2/m3, C_f: dimensionless, d_p: m, sigma_c: N/m)
if packing_type == "1.5-inch Berl Saddles":
    a_t = 110.0
    C_f = 60.0
    d_p = 0.038
    sigma_c = 0.061
else:
    a_t = 190.0
    C_f = 137.0
    d_p = 0.025
    sigma_c = 0.061

# Global System Constraints
dP_limit = 100.0  # Allowable pressure drop in N/m2/m 
P_atm = 1.0       # Ambient pressure in atm

# ==============================================================================
# GLOBAL THERMODYNAMIC & RHEOLOGICAL FUNCTIONS
# ==============================================================================
T_K = T_C + 273.15
R_u_atm = 0.08206   # Universal gas constant (L*atm/mol*K)
R_u_kcal = 1.987    # Universal gas constant (kcal/kmol*K)
g = 9.81            # Gravitational acceleration (m/s^2)
MW_air = 28.966     # Molar mass of air (kg/kmol)
MW_water = 18.015   # Molar mass of water (kg/kmol)

# 1. Fluid Rheology (Polynomial Temperature Approximations)
# Aqueous phase density and dynamic viscosity
rho_L = 1000 * (1 - ((T_C + 288.9414) / (508929.2 * (T_C + 68.12963))) * (T_C - 3.9863)**2)
mu_L = 0.001 * np.exp(1.003 - 1.011 * T_C / 20.0) 
sigma_L = (75.6 - 0.14 * T_C) / 1000.0            
C_water = rho_L / MW_water                        

# Vapor phase density and dynamic viscosity (Ideal Gas Assumption)
rho_G = (P_atm * MW_air) / (0.08206 * T_K)
mu_G = 1.716e-5 * ((T_K / 273.15)**1.5) * (383.55 / (T_K + 110.4))

def solve_kinetics_for_compound(MW_voc, dH_kcal, van_t_hoff_C, V_B, sum_V_voc, C_in, C_out):
    """
    Computes Henry's Equilibrium, Molecular Diffusivities, Eckert Hydrodynamics, 
    and Onda Mass Transfer coefficients for a specific contaminant profile.
    """
    # A. Phase Equilibrium (van't Hoff Relationship)
    log10_H = (-dH_kcal / (R_u_kcal * T_K)) + van_t_hoff_C
    H_atm = 10 ** log10_H
    H_D = H_atm / (R_u_atm * T_K * C_water)
    
    # B. Molecular Diffusivity 
    # Liquid Phase (Hayduk-Laudie converted from cm2/s to m2/s)
    D_L_cm2_s = (13.26e-5) / ((mu_L * 1000)**1.14 * V_B**0.589)
    D_L = D_L_cm2_s * 1e-4
    
    # Gas Phase (Fuller-Schettler-Giddings converted from cm2/s to m2/s)
    sum_V_air = 20.1
    M_AB = (MW_air + MW_voc) / (MW_air * MW_voc)
    D_G_cm2_s = (1e-3 * T_K**1.75 * math.sqrt(M_AB)) / (P_atm * (sum_V_air**(1/3) + sum_V_voc**(1/3))**2)
    D_G = D_G_cm2_s * 1e-4

    # C. Hydrodynamics (Eckert GPDC at 100 N/m2/m constraint)
    G_to_L_molar = R_s * P_atm / H_atm
    L_to_G_molar = 1.0 / G_to_L_molar
    L_to_G_mass = L_to_G_molar * (MW_water / MW_air)
    
    # Flow Parameter (Abscissa X)
    X = L_to_G_mass * math.sqrt(rho_G / (rho_L - rho_G))
    
    # Capacity Parameter (Ordinate Y) via curve-fit regression
    ln_X = math.log(X)
    ln_Y = -3.386 - 1.081 * ln_X - 0.127 * (ln_X**2)
    Y = math.exp(ln_Y)
    
    # Solved Gas and Liquid Superficial Mass Fluxes
    G_prime = math.sqrt((Y * rho_G * (rho_L - rho_G)) / (C_f * (mu_L * 1000)**0.1))
    L_prime = G_prime * L_to_G_mass
    
    # D. Onda Mass Transfer Correlation
    # 1. Effective Wetted Area
    term1 = (sigma_c / sigma_L)**0.75
    term2 = (L_prime / (a_t * mu_L))**0.1
    term3 = ((L_prime**2 * a_t) / (rho_L**2 * g))**(-0.05)
    term4 = ((L_prime**2) / (rho_L * sigma_L * a_t))**0.2
    a_w = a_t * (1 - math.exp(-1.45 * term1 * term2 * term3 * term4))
    
    # 2. Local Transfer Coefficients
    kL_term1 = 0.005 * (L_prime / (a_w * mu_L))**(2/3)
    kL_term2 = (mu_L / (rho_L * D_L))**(-0.5)
    kL_term3 = (a_t * d_p)**0.4
    k_L = (kL_term1 * kL_term2 * kL_term3) / ((rho_L / (mu_L * g))**(1/3))
    
    kG_term1 = 5.23 * (G_prime / (a_t * mu_G))**0.7
    kG_term2 = (mu_G / (rho_G * D_G))**(1/3)
    k_G = (kG_term1 * kG_term2 * ((a_t * d_p)**(-2))) * (a_t * D_G)
    
    # 3. Overall Volumetric Coefficient
    K_L = 1 / ((1 / k_L) + (1 / (H_D * k_G)))
    K_L_a = K_L * a_w
    
    # E. Unit Sizing
    HTU = L_prime / (K_L_a * rho_L)
    NTU = (R_s / (R_s - 1)) * math.log(((C_in / C_out) * (R_s - 1) + 1) / R_s)
    Z_height = HTU * NTU
    
    return {
        'H_atm': H_atm, 'H_D': H_D, 'D_L': D_L, 'D_G': D_G, 
        'X': X, 'L_prime': L_prime, 'G_prime': G_prime, 
        'a_w': a_w, 'k_L': k_L, 'k_G': k_G, 'K_L_a': K_L_a, 
        'HTU': HTU, 'NTU': NTU, 'Z': Z_height
    }

# ==============================================================================
# MULTICOMPONENT OPTIMIZATION LOGIC (Z_max Approach)
# ==============================================================================

# Evaluate TCE (Data source: Kavanaugh and Trussell properties)
res_TCE = solve_kinetics_for_compound(
    MW_voc=131.389, dH_kcal=3410.0, van_t_hoff_C=8.59, 
    V_B=107.1, sum_V_voc=93.48, C_in=C_in_TCE, C_out=C_out_TCE
)

# Evaluate PCE (Data source: Approximated structural properties)
res_PCE = solve_kinetics_for_compound(
    MW_voc=165.83, dH_kcal=4100.0, van_t_hoff_C=9.10, 
    V_B=131.0, sum_V_voc=115.0, C_in=C_in_PCE, C_out=C_out_PCE
)

# 1. Determine Limiting Compound (Lowest Henry's Constant) dictates Area
limiting_res = res_TCE if res_TCE['H_atm'] < res_PCE['H_atm'] else res_PCE

Q_water_m3_s = Q_water_m3day / (24 * 3600)
mass_flow_water = Q_water_m3_s * rho_L
Area_required = mass_flow_water / limiting_res['L_prime']
Diameter_final = math.sqrt((4 * Area_required) / math.pi)

# 2. Extract Max Height (Z_max) ensuring compliance for both
Z_max = max(res_TCE['Z'], res_PCE['Z'])
Total_Volume = Area_required * Z_max

# ==============================================================================
# REACTIVE DASHBOARD RENDER
# ==============================================================================
st.markdown("---")
st.header("Final Multicomponent Sizing Specifications")

col1, col2, col3, col4 = st.columns(4)
col1.metric("Required Column Diameter", f"{Diameter_final:.3f} m")
col2.metric("Total Packed Height (Z_max)", f"{Z_max:.2f} m")
col3.metric("Total Packing Volume", f"{Total_Volume:.2f} m³")
col4.metric("Geometric Cross-Section", f"{Area_required:.3f} m²")

st.markdown("### Granular Kinetic Resolution (TCE vs PCE)")

# FIX 2: Added exact dictionary keys for the formatting engine to parse
table_md = f"""

| System Parameter | Trichloroethylene (TCE) | Tetrachloroethylene (PCE) |
| :--- | :--- | :--- |
| **Dimensionless Henry (H_D)** | {res_TCE['H_D']:.4f} | {res_PCE['H_D']:.4f} |
| **Liquid Diffusivity (D_L) [m²/s]** | {res_TCE['D_L']:.3e} | {res_PCE['D_L']:.3e} |
| **Gas Diffusivity (D_G) [m²/s]** | {res_TCE['D_G']:.3e} | {res_PCE['D_G']:.3e} |
| **Superficial Liquid Flux (L')** | {res_TCE['L_prime']:.3f} | {res_PCE['L_prime']:.3f} |
| **Effective Wetted Area (a_w)** | {res_TCE['a_w']:.2f} | {res_PCE['a_w']:.2f} |
| **Overall Volumetric (K_L a) [1/s]**| {res_TCE['K_L_a']:.4f} | {res_PCE['K_L_a']:.4f} |
| **Height of Transfer Unit (HTU)** | {res_TCE['HTU']:.3f} m | {res_PCE['HTU']:.3f} m |
| **Number of Transfer Units (NTU)** | {res_TCE['NTU']:.2f} | {res_PCE['NTU']:.2f} |
| **Required Depth (Z_i)** | **{res_TCE['Z']:.2f} m** | **{res_PCE['Z']:.2f} m** |
"""
st.markdown(table_md)
