import numpy as np
from CoolProp.CoolProp import PropsSI

# -----------------------------------------------------------------------------------
# CONSTANTES COMUNES
# -----------------------------------------------------------------------------------
g = 9.81                # m/s²
rho_H2O = 1000          # kg/m³
rho_diesel = 848        # kg/m³
PCI_diesel = 42590      # kJ/kg
C_p_air = 1.005         # kJ/kgK
C_e_comb = 1.922        # kJ/kgK
eta_ge = 0.80
d = 0.015               # m
Do = 0.0446             # m
gamma = 1.4
rho_air = 1.225         # kg/m³
mu_air = 1.8e-5         # Pa.s
P_atm_mmHg = 762.8
T_amb_C = 16.0
P_atm_Pa = P_atm_mmHg * 133.322
T_amb_K = T_amb_C + 273.15
T_b_sec_C = 15.2
T_b_hum_C = 11.2
P1 = P_atm_Pa
R_universal = 8.314462618  # J/molK

# Datos de emisiones y costos (punto h e i)
factor_emision_CO2 = 2.69  # kgCO2 por litro de diesel quemado
precio_litro_diesel = 958.17  # pesos/litro
precio_kWh_SIC = 244.12      # pesos/kWh

# -----------------------------------------------------------------------------------
# Datos de los ensayos
# -----------------------------------------------------------------------------------
ensayos = [
    {
        "V_diesel_ml": 2,
        "tau_diesel_s": 13.078,
        "N_imp": 0,
        "tau_imp_s": 0,
        "delta_h_cmH2O": 18.5,
        "T_gases_esc_C": 172
    },
    {
        "V_diesel_ml": 2,
        "tau_diesel_s": 7.646,
        "N_imp": 62,
        "tau_imp_s": 120,
        "delta_h_cmH2O": 17.8,
        "T_gases_esc_C": 252
    },
    {
        "V_diesel_ml": 2,
        "tau_diesel_s": 5.366,
        "N_imp": 114,
        "tau_imp_s": 113,
        "delta_h_cmH2O": 17.1,
        "T_gases_esc_C": 380
    }
]

# -----------------------------------------------------------------------------------
# Funciones auxiliares
# -----------------------------------------------------------------------------------
def presion_vapor_agua(T_K):
    T_C = T_K - 273.15
    return 610.78 * np.exp(17.2694 * T_C / (T_C + 237.3))

def mass_flow_iso_orifice(rho1, mu1, d, Do, gamma, P1, deltaP):
    beta = d / Do
    P2 = P1 + deltaP
    L1 = 0.0254 / Do
    Re_Do = 1e-6
    M_dot = 1e-3
    tol = 1e-8
    max_iter = 100

    for _ in range(max_iter):
        Re_Do = 4 * M_dot / (np.pi * mu1 * Do)
        A_term = ((19000 * beta) / Re_Do) ** 0.8
        M_p2 = (2 * L1) / (1 - beta)
        term1 = 0.5961 + 0.0261 * beta ** 2 - 0.216 * beta ** 8
        term2 = 0.000521 * (1e6 * beta / Re_Do) ** 0.7
        term3 = (0.0188 + 0.0063 * A_term) * beta ** 3.5 * (1e6 / Re_Do) ** 0.3
        exp1 = np.exp(-10 * L1)
        exp2 = np.exp(-7 * L1)
        term4 = (0.043 + 0.080 * exp1 - 0.123 * exp2) * (1 - 0.11 * A_term) * (beta ** 4 / (1 - beta ** 4))
        term5 = -0.031 * (M_p2 - 0.8 * M_p2 ** 1.1) * beta ** 1.3
        D_mm = Do * 1000
        correction = 0.0
        if D_mm < 71.12:
            correction = 0.011 * (0.75 - beta) * (2.8 - D_mm / 25.4)
        C = term1 + term2 + term3 + term4 + term5 + correction
        tau = P2 / P1
        epsilon = 1 - (0.351 + 0.256 * beta ** 4 + 0.93 * beta ** 8) * (1 - tau ** (1 / gamma))
        M_dot_new = (C / np.sqrt(1 - beta ** 4)) * epsilon * (np.pi / 4) * d ** 2 * np.sqrt(2 * deltaP * rho1)
        if abs(M_dot_new - M_dot) < tol:
            return M_dot_new, C, Re_Do, epsilon
        M_dot = M_dot_new
    return M_dot, C, Re_Do, epsilon

def calcular_Tadiabatica(M_dot_air, M_dot_diesel, y_O2, y_CO2, y_H2O, y_N2, 
                         h_air_Tr, h_comb_Tr, PCI, T_amb_K, P_atm_Pa, tol=0.1):
    Q_dot_r = M_dot_diesel * PCI
    Q_dot_a = M_dot_air * h_air_Tr
    Q_dot_c = M_dot_diesel * h_comb_Tr
    Q_dot_g = Q_dot_r - Q_dot_a - Q_dot_c
    M_dot_g = M_dot_air + M_dot_diesel
    T_low, T_high = 1500, 3000
    while (T_high - T_low) > tol:
        T_mid = (T_low + T_high) / 2
        h_O2 = PropsSI('H', 'T', T_mid, 'P', P_atm_Pa, 'O2') / 1000
        h_CO2 = PropsSI('H', 'T', T_mid, 'P', P_atm_Pa, 'CO2') / 1000
        h_H2O = PropsSI('H', 'T', T_mid, 'P', P_atm_Pa, 'Water') / 1000
        h_N2 = PropsSI('H', 'T', T_mid, 'P', P_atm_Pa, 'N2') / 1000
        h_gases = (y_O2 * h_O2) + (y_CO2 * h_CO2) + (y_H2O * h_H2O) + (y_N2 * h_N2)
        h_O2_amb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'O2') / 1000
        h_CO2_amb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'CO2') / 1000
        h_H2O_amb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'Water') / 1000
        h_N2_amb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'N2') / 1000
        h_gases_amb = (y_O2 * h_O2_amb) + (y_CO2 * h_CO2_amb) + (y_H2O * h_H2O_amb) + (y_N2 * h_N2_amb)
        Q_dot_gases = M_dot_g * (h_gases - h_gases_amb)
        if Q_dot_gases < Q_dot_g:
            T_low = T_mid
        else:
            T_high = T_mid
    return (T_low + T_high) / 2

# -----------------------------------------------------------------------------------
# CÁLCULOS PARA CADA ENSAYO
# -----------------------------------------------------------------------------------
for i, ensayo in enumerate(ensayos, start=1):

    # Flujo volumétrico y másico de diésel
    V_dot_diesel_m3_s = (ensayo["V_diesel_ml"] * 1e-6) / ensayo["tau_diesel_s"]
    M_dot_diesel = rho_diesel * V_dot_diesel_m3_s

    # Potencia eléctrica generada
    W_elec = ensayo["N_imp"] * 3600 / ensayo["tau_imp_s"] / 1000 if ensayo["N_imp"] > 0 else 0.0

    # Diferencia de presión convertida a Pa
    deltaP = ensayo["delta_h_cmH2O"] / 100 * rho_H2O * g

    # Calcular flujo másico de aire
    M_dot_air, C, Re_Do, epsilon = mass_flow_iso_orifice(rho_air, mu_air, d, Do, gamma, P1, deltaP)

    # Masas molares (g/mol)
    M_C = 12.0107
    M_H = 1.00784
    M_aire = 28.97

    # Composición del combustible (C16H34)
    m_C = 16
    n_H = 34
    A_est = m_C + n_H / 4

    # Moles estequiométricos
    x_est = m_C
    y_est = n_H / 2
    z_est = A_est * 3.76

    # Calcular RAC estequiométrico
    masa_combustible_g = m_C * M_C + n_H * M_H  # g/mol
    RAC_est = (A_est * (1 + 3.76) * M_aire) / masa_combustible_g  # [g aire / g combustible]

    # Calcular RAC experimental
    RAC_exp = M_dot_air / M_dot_diesel

    # Calcular exceso de aire
    exceso_aire = (RAC_exp / RAC_est - 1) * 100
    
    # Cálculo de humedad relativa y n_w (vapor de agua ambiente)
    T_b_sec_K = T_b_sec_C + 273.15
    T_b_hum_K = T_b_hum_C + 273.15
    P_ws = presion_vapor_agua(T_b_sec_K)
    P_w = P_ws - (P_atm_Pa * (T_b_sec_K - T_b_hum_K) * 0.00066 * (1 + 0.00115 * T_b_hum_K))
    humedad_relativa = (P_w / P_ws) * 100

    # n_w en los productos
    n_aire_total = (1 + (RAC_exp / RAC_est - 1)) * A_est * (1 + 3.76)
    n_w = (P_w * n_aire_total) / (P_atm_Pa - P_w)

    # Moles reales de productos (con exceso de aire y humedad)
    e = (RAC_exp / RAC_est) - 1
    n_CO2 = x_est
    n_H2O = y_est + n_w
    n_O2 = e * A_est
    n_N2 = (1 + e) * A_est * 3.76

    # Porcentajes en volumen medidos por un analizador (base seca, sin H2O)
    n_total_seco = n_CO2 + n_O2 + n_N2
    vol_O2 = 100 * n_O2 / n_total_seco
    vol_CO2 = 100 * n_CO2 / n_total_seco
    vol_N2 = 100 * n_N2 / n_total_seco

    # Masas molares en g/mol
    MM_O2 = 31.9988
    MM_CO2 = 44.0095
    MM_H2O = 18.01528
    MM_N2 = 28.0134

    # Calcula la masa total de productos en gramos
    M_total_products = (n_O2 * MM_O2 +
                n_CO2 * MM_CO2 +
                n_H2O * MM_H2O +
                n_N2 * MM_N2)
    
    M_total_products_kg = M_total_products / 1000  # kg

    # Fracciones másicas de gases
    y_O2   = (n_O2 * MM_O2 / 1000) / M_total_products_kg
    y_CO2  = (n_CO2 * MM_CO2 / 1000) / M_total_products_kg
    y_H2O  = (n_H2O * MM_H2O / 1000) / M_total_products_kg
    y_N2   = (n_N2 * MM_N2 / 1000) / M_total_products_kg

    # Calor para llevar reactantes a 25°C
    h_air_Tr = C_p_air * (25)
    h_comb_Tr = C_e_comb * (25)

    # Cálculos de energías (en kW)
    Q_dot_r = M_dot_diesel * PCI_diesel
    Q_dot_a = M_dot_air * h_air_Tr
    Q_dot_c = M_dot_diesel * h_comb_Tr
    Q_dot_g = Q_dot_r - Q_dot_a - Q_dot_c

    # Temperatura adiabática
    T_ad_K = calcular_Tadiabatica(M_dot_air, M_dot_diesel, n_O2, n_CO2, n_H2O, n_N2,
                                  h_air_Tr, h_comb_Tr, PCI_diesel, T_amb_K, P_atm_Pa)

    # Flujo de energía de gases de escape usando las diferencias de entalpía
    M_dot_g = M_dot_air + M_dot_diesel

    # Entalpías específicas en la temperatura de salida (kJ/kg)
    h_O2_Tad = PropsSI('H', 'T', T_ad_K, 'P', P_atm_Pa, 'O2') / 1000
    h_CO2_Tad = PropsSI('H', 'T', T_ad_K, 'P', P_atm_Pa, 'CO2') / 1000
    h_H2O_Tad = PropsSI('H', 'T', T_ad_K, 'P', P_atm_Pa, 'Water') / 1000
    h_N2_Tad = PropsSI('H', 'T', T_ad_K, 'P', P_atm_Pa, 'N2') / 1000

    # Entalpías específicas a temperatura ambiente (kJ/kg)
    h_O2_Tamb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'O2') / 1000
    h_CO2_Tamb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'CO2') / 1000
    h_H2O_Tamb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'Water') / 1000
    h_N2_Tamb = PropsSI('H', 'T', T_amb_K, 'P', P_atm_Pa, 'N2') / 1000

    # Diferencias de entalpía (kJ/kg)
    delta_h_O2 = h_O2_Tad - h_O2_Tamb
    delta_h_CO2 = h_CO2_Tad - h_CO2_Tamb
    delta_h_H2O = h_H2O_Tad - h_H2O_Tamb
    delta_h_N2 = h_N2_Tad - h_N2_Tamb

    # Energía de los gases (kW) usando fracciones másicas
    h_gases_diff = (y_O2 * delta_h_O2 +
                    y_CO2 * delta_h_CO2 +
                    y_H2O * delta_h_H2O +
                    y_N2 * delta_h_N2)

    Q_dot_gases = M_dot_g * h_gases_diff  # kW

    # Potencia al eje del motor
    W_eje = W_elec / eta_ge

    # Rendimiento global
    eta_global = (W_elec / (M_dot_diesel * PCI_diesel)) * 100

    # Resultados
    print(f"\n=== ENSAYO {i} ===")
    print(f"Flujo másico aire: {M_dot_air:.7f} kg/s")
    print(f"Flujo másico diésel: {M_dot_diesel:.7f} kg/s")
    print(f"Razón aire-combustible (RAC): {RAC_exp:.7f}")
    print(f"RAC estequiométrico: {RAC_est:.7f}")
    print(f"Exceso de aire: {exceso_aire:.7f} %")
    print(f"Humedad relativa del aire: {humedad_relativa:.2f} %")

    print("Moles estequiométricos de productos:")
    print(f"  CO2: {x_est:.2f} moles")
    print(f"  H2O: {y_est:.2f} moles")
    print(f"  N2: {z_est:.2f} moles")
    print(f"  O2: 0.00 moles (todo reacciona en teoría)")

    print("Moles reales de productos (con exceso de aire y humedad):")
    print(f"  CO2: {n_CO2:.2f} moles")
    print(f"  H2O: {n_H2O:.2f} moles")
    print(f"  N2: {n_N2:.2f} moles")
    print(f"  O2: {n_O2:.2f} moles (exceso)")

    print("Composición que mediría el analizador de gases (base seca):")
    print(f"  O2: {vol_O2:.2f} %")
    print(f"  CO2: {vol_CO2:.2f} %")
    print(f"  N2: {vol_N2:.2f} %")

    print("Balance energético (en kW):")
    print(f"  Q̇_r (calor liberado por la combustión): {Q_dot_r:.2f} kW")
    print(f"  Q̇_a (calor para precalentar el aire): {Q_dot_a:.2f} kW")
    print(f"  Q̇_c (calor para precalentar el combustible): {Q_dot_c:.2f} kW")
    print(f"  Q̇_g (balance global, método 1) = {Q_dot_g:.2f} kW")
    print(f"  Q̇_gases (diferencias entalpías, método 2) = {Q_dot_gases:.2f} kW")

    print(f"T adiabática de combustión: {T_ad_K - 273.15:.7f} °C")

    print(f"Potencia eje motor: {W_eje:.2f} kW")
    print(f"Potencia eléctrica generada: {W_elec:.2f} kW")
    print(f"Rendimiento global: {eta_global:.2f} %")

# -----------------------------------------------------------------------------------
# PUNTO (h) Consumo mensual, emisiones y costo económico
# -----------------------------------------------------------------------------------
energia_mensual_kWh = 3000  # kWh/mes
energia_mensual_kJ = energia_mensual_kWh * 3600  # kJ
ef_gp = 0.2698
rho_diesel_L = 0.848  # kg/L

# Consumo específico de diésel (litros/kWh)
energia_diesel_kJ_por_litro = PCI_diesel * rho_diesel_L  # kJ/litro
litros_diesel_mes = energia_mensual_kJ / ( energia_diesel_kJ_por_litro * ef_gp )

# Emisiones de CO2
emisiones_CO2_mes_kg = litros_diesel_mes * factor_emision_CO2

# Costo económico mensual
costo_diesel_mes = litros_diesel_mes * precio_litro_diesel

print("\n=== PUNTO (h): Consumo mensual grupo electrógeno ===")
print(f"Litros diésel requeridos: {litros_diesel_mes:.3f} litros")
print(f"Emisiones de CO2: {emisiones_CO2_mes_kg:.3f} kgCO2")
print(f"Costo mensual diésel: ${costo_diesel_mes:.3f} pesos")

# -----------------------------------------------------------------------------------
# PUNTO (i) Consumo mensual y costo con red SEN
# -----------------------------------------------------------------------------------
costo_SEN_mes = energia_mensual_kWh * precio_kWh_SIC
emisiones_CO2_SEN = energia_mensual_kWh * 0.2421  # 0.2421 kgCO2/kWh

print("\n=== PUNTO (i): Consumo mensual con SEN ===")
print(f"Costo mensual con SEN: ${costo_SEN_mes:.3f} pesos")
print(f"Emisiones de CO2 con SEN: {emisiones_CO2_SEN:.3f} kgCO2")

