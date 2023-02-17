#########################################################################################################################################################################
################################################### This is the code to apply SEBS to the images of Italy ###############################################################
#########################################################################################################################################################################

#*********************************************** Calculation of NDVI, Albedo and Emissivity from Planet images ********************************************************** 
import numpy as np
import openpyxl as xl
import rasterio as rio
from rasterio.plot import show
from metpy.units import units
import pandas as pd
from matplotlib import pyplot
import scipy
from osgeo import gdal
import math
# from valensat.github.io.heat_fluxes_italo import KB_1, z0h, Cw, PSIh_y, u_pbl, FRUstar

#Formulas (esto hay que resolverlo)
import numpy as np
from sklearn import linear_model
from rasterio.plot import show

#SEBAL and METRIC functions

def sensible_heat(RN,LST_cold,LST_hot,LST, rah, ustar, air_dens, coldX, coldY, hotX, hotY,LET_cold):
    # Monin-Obukhov length (m):
    k_vk = 0.41
    ro = 1.23  # [kg/m3] dry air density
    Cp = 1004  # [JK/kg] air specific heat
    RNHot = float(RN[[hotX], :][:, [hotY]])
    rahHot = float(rah[[hotX], :][:, [hotY]])
    RN_cold = float(RN[[coldX], :][:, [coldY]])
    #G0_cold = float(G0[[coldX], :][:, [coldY]])
    # Hhot_h = RNHot-GHot #hourly
    Hhot = RNHot  # daily
    if LET_cold == 0:
        Hcold = 0
    else:
        Hcold = RN_cold - LET_cold
    rahCold = float(rah[[coldX], :][:, [coldY]])
    if str(rahHot == 'nan'):
        rahHot = np.nanmean(rah)
        if str(rahHot == 'nan'):
            rahHot = 110  # aerodynamic resistance for dry soil Qiu et al 1998
    dThot = (Hhot * rahHot) / (ro * Cp)  # dT daily first iteration
    dTcold = (Hcold * rahCold) / (ro * Cp)  # dT daily first iteration
    regr = linear_model.LinearRegression()
    XLST = np.zeros((2, 1))
    XLST[0] = LST_cold
    XLST[1] = LST_hot
    YdT_d = np.zeros((2, 1))
    YdT_d[0] = dTcold
    YdT_d[1] = dThot
    regr.fit(XLST, YdT_d)
    aD = regr.coef_
    bD = regr.intercept_
    dT = aD * LST + bD  # dT daily

    h = air_dens * 1004 * dT / rah
    L_MO = ((-1.0 * air_dens * 1004 * np.power(ustar, 3) * LST) /
            (k_vk * 9.81 * h))
    L_MO[L_MO < -1000] = -1000

    # Stability correction for momentum, stable conditions (L_MO >= 0):
    psi_200_stable = -0.05 * 200 / L_MO
    # Stability correction for momentum and heat transport, unstable
    # conditions (L_MO < 0):
    x2 = np.power((1.0 - 16.0 * (2.0 / L_MO)), 0.25)  # x at 2m
    x2[np.isnan(x2)] = 1
    x200 = np.power(1.0 - 16.0 * (200 / L_MO), 0.25)  # x at 200m
    x200[np.isnan(x200)] = 1
    psi_h = 2 * np.log((1 + np.power(x2, 2)) / 2)
    psi_m200 = (2 * np.log((1 + x200) / 2) + np.log((1 + np.power(x200, 2)) /
                                                    2) - 2 * np.arctan(x200) + 0.5 * np.pi)
    print('Sensible Heat ', np.nanpercentile(h, 50))
    print('dT', np.nanpercentile(dT, 50))
    return L_MO, psi_200_stable, psi_h, psi_m200, h, dT

def Iterate_Friction_Velocity(RN,COLD,HOT,LST,rah,uF,ro,coldX,coldY,hotX,hotY,u_200, Surf_roughness, L, psi, psi_m200,
                              psi_m200_stable,LET_cold):
    """
    Function to correct the windspeed and aerodynamic resistance for the iterative process the output can be used as the new input for this model
    """
    k_vk = 0.41
    air_dens = 1.23
    # Sensible heat 2 (Step 6)
    # Corrected value for the friction velocity, unstable
    ustar_corr_unstable = (k_vk * u_200 / (np.log(200.0 / Surf_roughness) -
                                           psi_m200))
    # Corrected value for the friction velocity, stable
    ustar_corr_stable = (k_vk * u_200 / (np.log(200.0 / Surf_roughness) -
                                         psi_m200_stable))
    ustar_corr = np.where(L > 0.0, ustar_corr_stable, ustar_corr_unstable)
    ustar_corr[ustar_corr < 0.02] = 0.02

    rah_corr_unstable = (np.log(2.0 / 0.01) - psi) / (k_vk * ustar_corr)  # unstable
    rah_corr_stable = (np.log(2.0 / 0.01) - 0.0) / (k_vk * ustar_corr)  # stable

    L[np.isnan(L)] = -1000
    rah_corr = np.where(L > 0.0, rah_corr_stable, rah_corr_unstable)
    rah_corr[rah_corr > 477] = 477 #Jia references values
    rah_corr[rah_corr < 110] = 110 #Qiu reference values
    rah_corr[np.isnan(rah_corr)] = np.nanmean(rah_corr)

    print('Aerodynamic resistance 10:', np.nanpercentile(rah_corr, 10))
    print('Aerodynamic resistance 50:', np.nanpercentile(rah_corr, 50))
    print('Aerodynamic resistance 90:', np.nanpercentile(rah_corr, 90))

    L_corr, psi_m200_corr_stable, psi_corr, psi_m200_corr, h, dT = sensible_heat(RN, COLD, HOT, LST, rah_corr, ustar_corr, ro, coldX, coldY, hotX,
                                                                  hotY,LET_cold)
    print('Sensible Heat percentile 10:', np.nanpercentile(h, 10))
    print('Sensible Heat percentile 50:', np.nanpercentile(h, 50))
    print('Sensible Heat percentile 90:', np.nanpercentile(h, 90))
    # return L_MO, psi_200_stable, psi_h, psi_m200, h, dT
    return (L_corr, psi_corr, psi_m200_corr, psi_m200_corr_stable, h, ustar_corr, rah_corr, dT)


#SEBS
# SEBS EVAPOTRANSPIRATION
# MOS STABILITY CORRECTION FUNCTIONS
def PSIma(f, g):
    a = 0.33
    b = 0.41
    pi = 3.141592654
    pre_tangens = (np.arctan((2.0 * g - 1.0) / (3.0) ** (1. / 2.))) * pi / 180
    # tangens = ifthenelse(tangens > pi/2.0, tangens - 2.0 * pi, tangens)
    # Initial conditional mask
    Mask_tangens1 = np.copy(pre_tangens)
    Mask_tangens1[Mask_tangens1 > np.pi / 2.0] = -9999.0
    Mask_tangens1[Mask_tangens1 != -9999.0] = 0
    Mask_tangens1[Mask_tangens1 == -9999.0] = 1
    tangens1 = np.copy(pre_tangens)
    tangens1[tangens1 == -np.inf] = 0
    tangens1[tangens1 == np.inf] = 0
    tangens1 = tangens1 * Mask_tangens1
    tangens1 = np.nan_to_num(tangens1)
    # Second conditional mask
    Mask_tangens2 = np.copy(pre_tangens)
    Mask_tangens2[Mask_tangens2 > np.pi / 2.0] = -9999.0
    Mask_tangens2[Mask_tangens2 != -9999.0] = 1
    Mask_tangens2[Mask_tangens2 == -9999.0] = 0
    tangens2 = np.copy(pre_tangens)
    tangens2[tangens2 == -np.inf] = 0
    tangens2[tangens2 == np.inf] = 0
    tangens2 = tangens2 * Mask_tangens2
    tangens2 = np.nan_to_num(tangens2)
    # SUM CONDITIONAL
    tangens = tangens1 + tangens2
    PSIma = np.log(a + f) - 3.0 * b * f ** (1.0 / 3.0) + b * a ** (1.0 / 3.0) / 2.0 * np.log(
        (1 + g) ** 2.0 / (1.0 - g + (g) ** 2)) + (3.0) ** (1. / 2.) * b * a ** (1.0 / 3.0) * tangens
    return PSIma


def PSIm_y(Y):
    # Integrated stability correction function for momentum
    # Inputs
    # Y = -z/L, where z is the height, L the Obukhov length
    # test values

    # Constants (Brutsaert, 1999)
    a = 0.33
    b = 0.41
    # m = 1.0
    pi = 3.141592654

    # Calculation
    # //HK 040902
    Y = abs(Y)  # abs(Y)
    x = (Y / a) ** (1.0 / 3.0)
    PSI0 = np.log(a) + (3.0) ** (1. / 2.) * b * a ** (1.0 / 3.0) * pi / 6.0
    b_3 = b ** -3.0
    # PSIm_y = ifthenelse(Y <= b_3, PSIma(Y, x) + PSI0, PSIma(b_3, ((b_3/a)**(1.0/3.0))) + PSI0)
    # PSIm_y = ifthenelse(Y <= b_3, PSIma(Y, x) + PSI0, (1.0 / (PSIma(b_3, ((b_3/a)**(1.0/3.0))))) + PSI0)
    # Initial conditional mask
    Mask_PSIm_y1 = np.copy(Y)
    Mask_PSIm_y1[Mask_PSIm_y1 <= b_3] = -9999.0
    Mask_PSIm_y1[Mask_PSIm_y1 != -9999.0] = 0
    Mask_PSIm_y1[Mask_PSIm_y1 == -9999.0] = 1
    PSIm_y1 = PSIma(Y, x) + PSI0
    PSIm_y1[PSIm_y1 == -np.inf] = 0
    PSIm_y1[PSIm_y1 == np.inf] = 0
    PSIm_y1 = PSIm_y1 * Mask_PSIm_y1
    PSIm_y1 = np.nan_to_num(PSIm_y1)
    # Second conditional mask
    Mask_PSIm_y2 = np.copy(Y)
    Mask_PSIm_y2[Mask_PSIm_y2 <= b_3] = -9999.0
    Mask_PSIm_y2[Mask_PSIm_y2 != -9999.0] = 1
    Mask_PSIm_y2[Mask_PSIm_y2 == -9999.0] = 0
    PSIm_y2 = PSIma(b_3, ((b_3 / a) ** (1.0 / 3.0))) + PSI0
    PSIm_y2 = PSIm_y2 * Mask_PSIm_y2
    PSIm_y2 = np.nan_to_num(PSIm_y2)
    # SUM CONDITIONAL
    PSIm_y = PSIm_y1 + PSIm_y2

    return PSIm_y


def PSIh_y(Y):
    # Integrated stability correction function for heat
    # Inputs
    # Y = -z/L, z is the height, L the Obukhov length
    # constants (Brutsaert, 1999)
    c = 0.33
    d = 0.057
    n = 0.78
    # Calculation
    Y = abs(Y)
    PSIh_y = (1.0 - d) / n * np.log((c + Y ** n) / c)
    return PSIh_y


# BAS STABILITY CORRECTION FUNCTIONS
def Bw(hi, L, z0):
    # constants (Brutsaert, 1999)
    alfa = 0.12
    beta = 125.0

    # calculations
    B0 = (alfa / beta) * hi
    B1 = -1.0 * z0 / L
    B11 = -alfa * hi / L
    B21 = hi / (beta * z0)
    B22 = -beta * z0 / L
    tempB11 = PSIm_y(B11)
    tempB1 = PSIm_y(B1)
    # B = ifthenelse(z0 < B0, -1.0 * np.log(alfa) + PSIm_y(B11) - PSIm_y(B1), np.log(B21) + PSIm_y(B22) - PSIm_y(B1))
    # Initial conditional mask
    Mask_B1 = np.copy(z0)
    Mask_B1[Mask_B1 <= B0] = -9999.0
    Mask_B1[Mask_B1 != -9999.0] = 0
    Mask_B1[Mask_B1 == -9999.0] = 1
    B_1 = -1.0 * np.log(alfa) + PSIm_y(B11) - PSIm_y(B1)
    B_1[B_1 == -np.inf] = 0
    B_1[B_1 == np.inf] = 0
    B_1 = B_1 * Mask_B1
    B_1 = np.nan_to_num(B_1)
    # Second conditional mask
    Mask_B2 = np.copy(z0)
    Mask_B2[Mask_B2 <= B0] = -9999.0
    Mask_B2[Mask_B2 != -9999.0] = 1
    Mask_B2[Mask_B2 == -9999.0] = 0
    B_2 = np.log(B21) + PSIm_y(B22) - PSIm_y(B1)
    B_2[B_2 == -np.inf] = 0
    B_2[B_2 == np.inf] = 0
    B_2 = B_2 * Mask_B2
    B_2 = np.nan_to_num(B_2)
    # SUM CONDITIONAL
    B = B_1 + B_2
    # Bw = ifthenelse(B < 0.0, 0.0, B) # This results from unfortunate parameter combination!
    Bw = np.copy(B)
    Bw[Bw < 0.0] = 0  # This results from unfortunate parameter combination!

    return Bw


def Cw(hi, L, z0, z0h):
    alfa = 0.12
    beta = 125.0
    C0 = (alfa / beta) * hi
    # C1 = pcrumin(z0h) / L
    z0h[z0h == -np.inf] = 9999
    z0h[z0h == np.inf] = 9999
    C1 = np.nanmin(z0h) / L
    C11 = -alfa * hi / L
    C21 = hi / (beta * z0)
    C22 = -beta * z0 / L
    # C = ifthenelse(z0 < C0, pcrumin(ln(alfa)) + PSIh_y(C11) - PSIh_y(C1), ln(C21) + PSIh_y(C22) - PSIh_y(C1))
    # C = ifthenelse(z0 < C0, np.log(alfa) + PSIh_y(C11) - PSIh_y(C1), np.log(C21) + PSIh_y(C22) - PSIh_y(C1))
    # Initial conditional mask
    Mask_C1 = np.copy(z0)
    Mask_C1[Mask_C1 <= C0] = -9999.0
    Mask_C1[Mask_C1 != -9999.0] = 0
    Mask_C1[Mask_C1 == -9999.0] = 1
    C_1 = np.log(alfa) + PSIh_y(C11) - PSIh_y(C1)
    C_1[C_1 == -np.inf] = 0
    C_1[C_1 == np.inf] = 0
    C_1 = C_1 * Mask_C1
    C_1 = np.nan_to_num(C_1)
    # Second conditional mask
    Mask_C2 = np.copy(z0)
    Mask_C2[Mask_C2 <= C0] = -9999.0
    Mask_C2[Mask_C2 != -9999.0] = 1
    Mask_C2[Mask_C2 == -9999.0] = 0
    C_2 = np.log(C21) + PSIh_y(C22) - PSIh_y(C1)
    C_2[C_2 == -np.inf] = 0
    C_2[C_2 == np.inf] = 0
    C_2 = C_2 * Mask_C2
    C_2 = np.nan_to_num(C_2)
    # SUM CONDITIONAL
    C = C_1 + C_2
    # Cw = ifthenelse(C < 0.0, 0.0, C) # This results from unfortunate parameter combination!
    Cw = np.copy(C)
    Cw[Cw < 0.0] = 0  # This results from unfortunate parameter combination!

    return Cw


def u_pbl(u_s, z_ms, z0m, z_pbl):  # wind speed, height of wind speed measurement, z0m, z_pbl
    """Calculates Planetary Boundary Layer wind speed [m s-1] from Fcover"""
    h = z0m / 0.136  # total height of vegetation (m)
    d = 2.0 / 3.0 * h  # zero plane displacement (m)
    u_c = np.log((z_pbl - d) / z0m) / np.log((z_ms - d) / z0m)
    u_pbl = u_s * u_c
    return u_pbl, d, h


def z0h(KB_1, z0m):
    """Calculates the scalar roughness height for heat transfer (z0h)
    KB_1 Input KB_1 values
    z0m Input scalar roughness height for momentum"""
    z0h = z0m / np.exp(KB_1)
    return z0h


def KB_1(u, uF, Zom, Fc, LAI, TaK, pa, h):
    """Initial determination of roughness length for heat transfer
    #FRICTION VELOCITY IN CONDITION OF STABILITY"""
    k = 0.41
    Ct = 0.00625 * 2  # it asumes a medium value of 0.00625, and two because the sides of leaf, (Su, 2002)
    Pr = 0.713  # Prandt Number
    Cd = 0.2  # Drag coefficient of the foliage elements
    hs = 5e-03
    # uFi= (k*u)/np.log(2/Zom)
    z = 2  # Height of wind measurement
    # AN EXTENDED MODEL FOR DETERMINATION OF THE ROUGHNESS LENGTH FOR HEAT TRANSFER
    z0 = 0.136 * h  # Brutsaert (1982)
    u_h = u * np.log(2.446) / np.log((z - 0.667 * h) / z0)  # wind speed at canopy height
    ust2u_h = 0.32 - 0.264 / np.exp(15.1 * Cd * LAI)  # wind speed friction at canopy height
    nec = (Cd * LAI) / ((2 * (uF ** 2)) / u_h ** 2)  # Within-canopy wind profile extinction coefficient
    Fs = 1 - Fc  # Fraction of Soil Cover

    v = 1.327 * (10 ** -5) * (101.3 / pa) * ((TaK / 273.15) ** 1.81)  # kinematic viscosity of air, Su, 2002 based on Massman,1999

    Re = (hs * uF) / v  # Su, 2002
    Cst = (Pr ** (-2. / 3)) * (Re ** (-1. / 2))  # Su, 2002
    kBs = 2.46 * ((Re) ** (1. / 4)) - np.log(7.4)  # Su, 2002
    kB = ((k * Cd) / (4 * Ct * ust2u_h * (1 - np.exp((-nec / 2)))) * (Fc ** 2)) + (2 * Fc * Fs * ((k * ust2u_h * (Zom / h)) / Cst)) + (kBs * (Fs ** 2))  # B = Inverse Stanton number, a dimensionless heat transfer coefficient
    show(kB)
    return kB


def FRUstar(z_pbl, u_pbl, hst, z0m, z0h, DEM, HR, Ta, Ts):
    """Iteration to calculate RUstar OR FRICTION VELOCITY
    z_pbl Input PBL depth [m]
    u_pbl Input Wind speed PBL [m/s]
    hst Input height of the ASL [m]
    h Input height of vegetation
    Air Temperature
    Relative Humidity
    Land Surface Temperature"""
    print("Starting iterations to derive stability parameters...")
    g=9.8
    k=0.41
    p_s = 101325  # Pa
    z_pbl_A = z_pbl
    alt_ms = 2  # height of measurement
    d = z0m * 4.9
    u_pbl_A = u_pbl
    # h = z0m / 0.136
    Rv = 461.05  # specific gas constant water vapour (J kg-1 K-1)
    Rd = 287.04  # specific gas constant dry air (J kg-1 K-1)
    Cp = 1005.0  # specific heat (J kg-1 K-1)
    p_pbl_A = p_s * ((44331.0 - (DEM + z_pbl_A)) / (44331.0 - alt_ms)) ** (1.0 / 0.1903)

    zd0 = z_pbl - d
    ku = 0.41 * u_pbl_A
    zdm = np.log(zd0 / z0m)
    zdh = np.log(zd0 / z0h)
    zdh[zdh == np.inf] = 375579840.0  # high resistance
    zdh[np.isnan(zdh)] = 375579840
    RUstar = ku / zdm

    helpvar1 = DEM / 44331.0
    helpvar2 = 1.0 - helpvar1
    T0 = Ts / helpvar2 ** 1.5029
    t_s = Ta
    hr_s = HR
    hr_pbl = hr_s / 100
    Tcn = Ta - 273.15
    esat = 611.0 * np.exp(17.502 * Tcn / (Tcn + 240.97))  # Pa
    eact = hr_pbl * esat  # actual vapour pressure
    ps = p_pbl_A
    t_c = np.log((z_pbl - d) / z0h) / np.log((alt_ms - d) / z0h)
    t_pbl_A = Ts * (1.0 - t_c) + t_s * t_c
    t_pbl_A = t_pbl_A / (1.0 - DEM / 44331.0) ** 1.5029
    t_pbl = t_pbl_A
    # T_0pbl = 0.5 * (T0 + t_pbl_A)   # mean potential temperature
    q_pbl_A = 5.0 / 8.0 * eact / p_pbl_A

    p_pbl = p_pbl_A
    q_pbl = q_pbl_A
    eact = p_pbl * q_pbl * (Rv / Rd)  # actual vapour pressure

    Theta_s = T0
    Theta_v = Ts * (1.0 + 0.61 * q_pbl)  # surface virtual temperature
    Theta_a = t_pbl * (101325 / p_pbl) ** 0.286
    rhoa = ps / (Rd * Theta_v)  # surface air density (kg m-3)
    rhoam = (ps / (Rd * Ts)) * (1.0 - 0.378 * eact / ps)  # moist air density (kg m-3)
    rhoacp = rhoa * Cp  # specific air heat capacity (J K-1 mÂ­3)
    T0ta = Theta_s - Theta_a
    CH = T0ta * k * rhoacp

    RH = CH * RUstar / zdh  # Initial sensible heat flux
    RH0 = RH
    Reps = 10.0
    Isteps = 0
    RHA = RH
    RHB = RH
    RH0A = RH0
    RH0B = RH0
    RUstarA = RUstar
    RUstarB = RUstar
    IstepsA = Isteps
    IstepsB = Isteps
    RepsA = Reps
    RepsB = Reps
    itNr = 100.0
    itThreshold = 0.01
    CL = np.nanmin(rhoam) * Cp * Theta_v / (k * g)  # Monin Obukhov Length without H and uF

    while RepsA > itThreshold and IstepsA < itNr:
        RLA = CL * RUstarA ** 3.0 / RHA
        tempBw = Bw(z_pbl, RLA, z0m)
        RUstarA = ku / (zdm - tempBw)
        tempCw = Cw(z_pbl, RLA, z0m, z0h)
        RHA = CH * RUstarA / (zdh - tempCw)
        RepsA = np.nanmax(abs(RH0A - RHA))
        # difa = abs(RH0A - RHA)
        # Min = np.nanmin(difa)
        # meandif = np.nanmean(difa)
        RH0A = RHA
        IstepsA = IstepsA + 1
        percentage = (IstepsA / itNr) * 100
        print("Iteration A:", int(percentage), "% completed\r")
        print("it Threshold: " + str(RepsA))

    while RepsB > itThreshold and IstepsB < itNr:
        RLB = CL * RUstarB ** 3.0 / RHB
        tempPSIm_y1 = PSIm_y(zd0 / RLB)
        tempPSIm_y2 = PSIm_y(z0m / RLB)
        RUstarB = ku / (zdm - tempPSIm_y1 + tempPSIm_y2)
        tempPSIh_y1 = PSIh_y(zd0 / RLB)
        tempPSIh_y2 = PSIh_y(z0h / RLB)
        RHB = CH * RUstarB / (zdh - tempPSIh_y1 + tempPSIh_y2)
        RepsB = np.nanmax(abs(RH0B - RHB))
        # difb = abs(RH0B - RHB)
        # meandif = np.nanmean(difb)
        # Min = np.nanmin(difb)
        RH0B = RHB
        IstepsB = IstepsB + 1
        percentage = (IstepsB / itNr) * 100
        print("Iteration B:", int(percentage), "% completed\r")
        print("it Threshold: " + str(RepsB))
    print
    MASKYRU1 = np.copy(z_pbl)
    MASKYRU2 = np.copy(z_pbl)
    MASKYRU1[MASKYRU1 >= hst] = 1
    MASKYRU1[MASKYRU1 < hst] = 0
    MASKYRU2[MASKYRU2 >= hst] = 0
    MASKYRU2[MASKYRU2 < hst] = 1
    RUstarA[RUstarA == -np.inf] = 0
    RUstarA[RUstarA == np.inf] = 0
    RUstarA1 = MASKYRU1 * RUstarA
    RUstar[RUstar == -np.inf] = 0
    RUstar[RUstar == np.inf] = 0
    RUstarB1 = MASKYRU2 * RUstar
    RUstar = RUstarA1 + RUstarB1
    # RUstar = ifthenelse(z_pbl >= hst, RUstarA, RUstarB)
    RLA[RLA == -np.inf] = 0
    RLA[RLA == np.inf] = 0
    RL1 = MASKYRU1 + RLA
    RLB[RLB == -np.inf] = 0
    RLB[RLB == np.inf] = 0
    RL2 = MASKYRU2 + RLB
    RL = RL1 + RL2
    # RL = ifthenelse(z_pbl >= hst, RLA, RLB)

    # dif = ifthenelse(z_pbl >= hst, difa, difb)
    return RUstar, RL





#Open the Planet image (mosaic created with QGIS) (geotiff) --- Modena sud (is the station that contain diamante images and the EC station)
#all the images are at 3 meter resolution and with CRS usgs 32632 - UTM 32N (the transformation has been done with QGIS)

image =rio.open(r"C:\Users\1\Desktop\ET_SEBS\20220713_SR_Planet_mosaic_ModenaSud.tif")

#naming the bands and converting to float64 values
band1 = image.read(1)/(10000)
band2 = image.read(2)/(10000)
band3 = image.read(3)/(10000)
band4 = image.read(4)/(10000)
band5 = image.read(5)/(10000)
band6 = image.read(6)/(10000)
band7 = image.read(7)/(10000)
band8 = image.read(8)/(10000)
#print the shape
print(band1.shape)

#open the mask image
# the masks are orribles, I assume there are no clouds (the time pass of planet and of the airplane is the same more or less)
# mask =rio.open(r"C:\Users\1\Desktop\ET_SEBS\OneDrive_2023-01-30\SMARTIES field campaign Burana - 2022\satellite data\planet\files\20220713_095824_68_2416_3B_udm2_clip.tif")
# mask1 = mask.read(1) #Band 1: clear mask (a value of “1” indicates the pixel is clear, a value of “0” indicates that the pixel is not clear and is snow, cloud, shadow...)


#create a function to apply the mask to the bands and replace it
def apply_mask_and_replace(band, mask):
    band_masked = band * mask
    #band_masked = np.nan_to_num(band_masked) #replaces any Not-a-Number (NaN) values in the array with a zero value
    return band_masked

# #Apply it
# band1_msk = apply_mask_and_replace(band1, mask1)
# band2_msk = apply_mask_and_replace(band2, mask1)
# band3_msk = apply_mask_and_replace(band3, mask1)
# band4_msk = apply_mask_and_replace(band4, mask1)
# band5_msk = apply_mask_and_replace(band5, mask1)
# band6_msk = apply_mask_and_replace(band6, mask1)
# band7_msk = apply_mask_and_replace(band7, mask1)
# band8_msk = apply_mask_and_replace(band8, mask1)


#NDVI calculation = (NIR-RED)/(NIR + RED); with NIR = band8 and RED = band6 for Planet images
NDVI = (band8 - band6)/(band8 + band6)
#NDVI = np.nan_to_num(NDVI)

# print("Mean NDVI:", np.nanmean(NDVI))
# print("Max NDVI:", np.nanmax(NDVI))
# print("Min NDVI:", np.nanmin(NDVI))

# Mean NDVI: 0.5795904
# Max NDVI: 0.9259408
# Min NDVI: -0.16762327

# pyplot.imshow(NDVI, 'gray', aspect = 'auto')
# pyplot.show()

# Create a binary water mask with NDVI (NDVI>0)
water_mask = NDVI > 0
water_mask = water_mask.astype(int) # Convert the mask to an integer array (water = 0, good values = 1)

#Apply the mask to the bands
band1_msk = apply_mask_and_replace(band1, water_mask)
band2_msk = apply_mask_and_replace(band2, water_mask)
band3_msk = apply_mask_and_replace(band3, water_mask)
band4_msk = apply_mask_and_replace(band4, water_mask)
band5_msk = apply_mask_and_replace(band5, water_mask)
band6_msk = apply_mask_and_replace(band6, water_mask)
band7_msk = apply_mask_and_replace(band7, water_mask)
band8_msk = apply_mask_and_replace(band8, water_mask)

#Now I calculate albedo and emissivity
#ALBEDO

#Formula de Drazen
#Este es el resultado del test: bias, standard deviation y RMSE   -0.000433150   0.00310082   0.00313092
#b1=442.5, b2=490, b3=560, b4=665, b8a=865
#albedoS2 = (b1*0.239855) + (b2*0.00138392) + (b3*0.0544263) + (b4*0.359608) + (b8a*0.313689) +0.00673963
albedo = (band1_msk*0.239855) + (band2_msk*0.00138392) + (band4_msk*0.0544263) + (band6_msk*0.359608) + (band8_msk*0.313689) +0.00673963
# Find all the elements equal to 0.00673963 (the old 0), and replace it with Nan value
albedo[albedo == 0.00673963] = [np.nan]
albedo = albedo

# print('albedo :', albedo.dtype)
# print("Mean albedo:", np.nanmean(albedo))
# print("Max albedo:", np.nanmax(albedo))
# print("Min albedo:", np.nanmin(albedo))

# Mean albedo: 0.16202799440626237
# Max albedo: 1.4656104005273753
# Min albedo: 0.04642329758883543

# pyplot.imshow(albedo, 'gray')
# pyplot.show()

#EMISSIVITY (Sobrino, 2008)
# NDVI_soil = 0.2
# NDVI_veg = 0.8
fvc = ((NDVI - 0.2) / (0.8 - 0.2)) ** 2  # Fraction of Vegetation Cover

fvc[fvc > 1] = [1]
fvc = fvc

fvc[fvc < 0] = [0]
fvc = fvc

# soil_emis = 0.971
# veg_emis = 0.982
emiss = (0.971 *(1-fvc)) + (0.982*fvc) #It have to be included between 0 (white body) and 1(black body)
#all the part around is = 0.978

# print('emiss :', emiss.dtype)

# mean_emiss = np.nanmean(emiss)
# maximum_emiss = np.nanmax(emiss)
# minimum_emiss = np.nanmin(emiss)
# print("Mean emiss:", mean_emiss)
# print("Max emiss:", maximum_emiss)
# print("Min emiss:", minimum_emiss)

# Mean emiss: 0.9764822
# Max emiss: 0.9871024
# Min emiss: 0.97099996

# pyplot.imshow(emiss, 'gray')
# pyplot.show()


#SAVE THE IMAGES

# Profile information for the GeoTIFF image
size1, size2 = image.shape
profile = {
    "driver": "GTiff",
    "height": image.shape[0],
    "width": image.shape[1],
    "count": 1,
    "dtype": np.float64,
    "transform": image.transform,
    "crs": 'EPSG:32632',
}



# Save the arrays as GeoTIFF images
with rio.open("C:/Users/1/Desktop/ET_git/ET_NDVI_220713_0958.tif", "w", **profile) as dst:
    dst.write(NDVI, 1)

with rio.open("C:/Users/1/Desktop/ET_git/ET_ALBEDO_220713_0958.tif", "w", **profile) as dst:
    dst.write(albedo, 1)

with rio.open("C:/Users/1/Desktop/ET_git/ET_EMISS_220713_0958.tif", "w", **profile) as dst:
    dst.write(emiss, 1)

#********************************************************************* Opening LST images ********************************************************************
#Open the LST image
LST = rio.open(r"C:\Users\1\Desktop\ET_SEBS\Modena_Sud_4m_8wn_L1-3_20220713_0825_3m.LST.tif")
LST = LST.read(1)

LST[LST == 0] = [np.nan]
LST = LST

#hay algun pixel que falla --- pongo una mascara
LST[LST > 340] = [np.nan]
LST = LST

LST[LST < 290] = [np.nan]
LST = LST

# Mean LST: 308.13271670014984
# Max LST: 339.69772951460845
# Min LST: 290.04785392586734


#********************************************************************* Opening ERA5 images ********************************************************************
#Open the ERA5 image and calculate other parameters

#Tref = t2m (2 metre temperature) ERA5 product --- 7 bands, one for every hour, from 8 to 17 --- units=K
t2m =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\t2m_3m_clip_ModenaSud.tif")

t2m_08 = t2m.read(2) #in this case i want to utilize the image at 08.00
t2m_08[t2m_08 == -32767] = [np.nan] #missing_value=-32767
t2m_08 = t2m_08


# print('t2m_08 :', t2m_08.dtype)
# print("Mean t2m_08:", np.nanmean(t2m_08))
# print("Max t2m_08:", np.nanmax(t2m_08))
# print("Min t2m_08:", np.nanmin(t2m_08))

# Mean t2m_08: 302.0274
# Max t2m_08: 303.77518
# Min t2m_08: 300.37268

# pyplot.imshow(t2m_08, 'gray')
# pyplot.show()


lai = np.sqrt(NDVI * (1.0 + NDVI) / (1.0 - NDVI + 1.0E-6))

# print('lai :', lai.dtype)
# print("Mean lai:", np.nanmean(lai))
# print("Max lai:", np.nanmax(lai))
# print("Min lai:", np.nanmin(lai))

# Mean lai: 1.6845813
# Max lai: 4.907052
# Min lai: 0.009609516

# pyplot.imshow(lai, 'gray')
# pyplot.show()


z0m = np.exp(-5.5 + 5.8*NDVI) #Tomelloso method, Bolle and Streckbach 1993
#z0m = 0.018*lai #brutsaert et al 1982 ---- sale bastante mal

# print('z0m :', z0m.dtype)
# print("Mean z0m:", np.nanmean(z0m))
# print("Max z0m:", np.nanmax(z0m))
# print("Min z0m:", np.nanmin(z0m))

# Mean z0m: 0.19184501
# Max z0m: 0.8784965
# Min z0m: 0.0015457977

# pyplot.imshow(z0m, 'gray')
# pyplot.show()


# hc = 0.15 * lai #SEBAL manual #muy mal
hc = z0m/0.136 #Brutsaert, 1982 #parece estar bien

# print('hc :', hc.dtype)
# print("Mean hc:", np.nanmean(hc))
# print("Max hc:", np.nanmax(hc))
# print("Min hc:", np.nanmin(hc))

# Mean hc: 1.4106233
# Max hc: 6.4595327
# Min hc: 0.011366159

# pyplot.imshow(hc, 'viridis')
# pyplot.show()


# WIND VELOCITY
# The u component represents the wind velocity in the east-west direction, and the v component represents the wind velocity in the north-south direction. To find the total wind velocity, we treat the u and v components as vectors and add them together using the Pythagorean theorem.
u10 =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\u10_3m_clip_ModenaSud.tif")

u10_08 = u10.read(2) #in this case i want to utilize the image at 08.00
u10_08[u10_08 == -32767] = [np.nan] #missing_value=-32767
u10_08 = u10_08

v10 =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\v10_3m_clip_ModenaSud.tif")

v10_08 = v10.read(2) #in this case i want to utilize the image at 08.00
v10_08[v10_08 == -32767] = [np.nan] #missing_value=-32767
v10_08 = v10_08

w_ref = np.sqrt(u10_08**2 + v10_08**2) # (m/s)

# print('w_ref :', hc.dtype)
# print("Mean w_ref:", np.nanmean(w_ref))
# print("Max w_ref:", np.nanmax(w_ref))
# print("Min w_ref:", np.nanmin(w_ref))

# Mean w_ref: 0.5987057
# Max w_ref: 0.9717973
# Min w_ref: 0.53785473

# pyplot.imshow(w_ref, 'viridis')
# pyplot.show()

#PRESSURE
press_ref = np.full_like(t2m_08, 101325) #Sea level standard atmospheric pressure (Pa)

#SURFACE PRESSURE
press_surf =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\sp_3m_clip_ModenaSud.tif") #Pa
press_surf_08 = press_surf.read(2) #in this case i want to utilize the image at 08.00
press_surf_08[press_surf_08 == -32767] = [np.nan] #missing_value=-32767
press_surf_08 = press_surf_08

# Mean press_surf: 98190.03
# Max press_surf: 100653.305
# Min press_surf: 96121.29


#DEWPOINT T
d2m =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\d2m_3m_clip_ModenaSud.tif") #Kelvin
d2m_08 = d2m.read(2) #in this case i want to utilize the image at 08.00
d2m_08[d2m_08 == -32767] = [np.nan] #missing_value=-32767
d2m_08 = d2m_08

# Mean d2m: 285.49936
# Max d2m: 286.72897
# Min d2m: 285.35596

#RELATIVE HUMIDITY
#where Td is the dew point temperature in K and T is the ambient temperature in K
#the equation assumes a standard pressure of 1013.25 hPa
#from WMO (World Meteorological Organization) Guide to Meteorological Instruments and Methods of Observation

def relative_humidity(T, Td):
    T = T - 273.15 #transformation to Celsius Degree
    Td = Td - 273.15 #transformation to Celsius Degree
    rh = 100 * (np.exp((17.625 * Td) / (243.04 + Td)) / np.exp((17.625 * T) / (243.04 + T)))
    return rh

q_ref = relative_humidity(t2m_08, d2m_08)  #[%]

# print('q_ref :', q_ref.dtype)
# print("Mean q_ref:", np.nanmean(q_ref))
# print("Max q_ref:", np.nanmax(q_ref))
# print("Min q_ref:", np.nanmin(q_ref))

# Mean q_ref: 36.110554
# Max q_ref: 39.784348
# Min q_ref: 34.35187

# pyplot.imshow(q_ref, 'viridis')
# pyplot.show()


# SWnet and LWnet
# remember: shortwave radiation is solar radiation, while longwave radiation is thermal radiation
# The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds
# to obtain SWin at 8: (SWin_08 - SWin_07)/3600).

SWin_08 =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\ssrd_3m_clip_ModenaSud.tif").read(2) #[J m-2]
SWin_08[SWin_08 == -32767] = [np.nan] #missing_value=-32767
SWin_08 = SWin_08

SWin_07 =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\ssrd_3m_clip_ModenaSud.tif").read(1) #[J m-2]
SWin_07[SWin_07 == -32767] = [np.nan] #missing_value=-32767
SWin_07 = SWin_07

SWin_final = (SWin_08 - SWin_07)/3600

LWin_08 =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\strd_3m_clip_ModenaSud.tif").read(2) #[J m-2]
LWin_08[LWin_08 == -32767] = [np.nan] #missing_value=-32767
LWin_08 = LWin_08

LWin_07 =rio.open(r"C:\Users\1\Desktop\ET_SEBS\ERA5\ERA5_20220713\strd_3m_clip_ModenaSud.tif").read(1) #[J m-2]
LWin_07[LWin_07 == -32767] = [np.nan] #missing_value=-32767
LWin_07 = LWin_07

LWin_final = (LWin_08 - LWin_07)/3600

sigma = 5.6704 * 10 ** -8  # boltzmann constant
SWnet = (1-albedo)*SWin_final
LWnet = (albedo*LWin_final) - (albedo*sigma*(LST**4))

# Mean SWin: 697.0475
# Max SWin: 702.196
# Min SWin: 689.7624

# Mean LWin: 374.59207
# Max LWin: 378.61472
# Min LWin: 374.28168

# MeanSWnet: 584.1046718435377
# MaxSWnet: 666.8312874127324
# MinSWnet: -325.0200448695023

# MeanLWnet: -22.27949224492215
# MaxLWnet: -3.7444079359536637
# MinLWnet: -172.75507708830335


# z_ref = 2m
z_ref = np.full_like(t2m_08, 2)

Rn = SWnet + LWnet #Net radiation [W/m2]
G0 = Rn * (0.05 + (1 - fvc) * (0.315 - 0.05)) #soil heat flux [W/m2]

# print ('Rn :', Rn.dtype)
# print("MeanRn:", np.nanmean(Rn))
# print("MaxRn:", np.nanmax(Rn))
# print("MinRn:", np.nanmin(Rn))

# # MeanRn: 564.2960743046193
# # MaxRn: 658.3841773306526
# # MinRn: -17.458414984320015

# pyplot.imshow(Rn, 'viridis')
# pyplot.show()

# print ('G0 :', G0.dtype)
# print("MeanG0:", np.nanmean(G0))
# print("MaxG0:", np.nanmax(G0))
# print("MinG0:", np.nanmin(G0))

# # MeanG0: 104.68610452158913
# # MaxG0: 206.82898139438467
# # MinG0: -5.483887374525834

# pyplot.imshow(G0, 'viridis')
# pyplot.show()


#****************************************************************** Rinomino le variabili di input **************************************************************************************

RN = Rn  # insert variable - Net Radiation [W/m2]
Ta = (t2m_08 - 273.15) #transformation to Celsius Degree  # insert variable - Air Temperature [Celsius]
LST = LST  # insert variable Land Surface Temperature [K]
LAI = lai  # insert variable LAI [m2/m2]
RH = q_ref  # Relative Humidity [%]
u = w_ref  # Wind speed at 2m [m/s]
TaK = t2m_08  # Air Temperature [K]
NDVI = NDVI  # insert - variable NDVI [-]
z0m = z0m
Fc = fvc
press_surf = press_surf_08


alt_ms = z_ref  # [m] Height of measurement
z_pbl = np.full_like(t2m_08, 1000)  # creo un array con valore 1000, [m] Height of Planetary Boundary Layer
hst = 0.12 * z_pbl  # Height of Atmospheric Boundary Layer Similary Stull.1988
z_ref = z_ref  # heightness reference
p_s = press_surf # Pa
ro = 1.23  # [kg/m3] dry air density
Cp = 1004  # [JK/kg] air specific heat
# DEFINITION OF ASL HEIGHT
alfa = 0.12  # Medium value of alpha defined by Brutsaert (1999)
beta = 125.0  # Medium value of beta defined by Brutsaert (1999)
hbpl = alfa * (beta / alfa)  # Height of Atmospheric Boundary Layer (ABL)
Cd = 0.2  # Drag coefficient of the foliage elements
Rv = 461.05  # specific gas constant water vapour (J kg-1 K-1)
Rd = 287.04  # specific gas constant dry air (J kg-1 K-1)
k = 0.41  # Von Karmann Constant
g = 9.81  # Gravity force [m/s2]
Alt = 150  # Altitude [m] - from google earth (manaro sul padano)
DOY = 194  # DAY OF THE YEAR
Cdi = -7 * (1e-06) * ((DOY)) ** 2 + 0.0027 * (DOY) + 0.124

es = 6.122 * np.exp((17.67 * (TaK - 273.16)) / ((TaK - 273.16) + 243.5))  # [HPa]
ea = (RH / 100) * es  # [HPa]
Ea = 1.24 * (ea / TaK) ** (1. / 7.)  # air emissivity
h = hc #1  # canopy height [1m]


#************************************************************************* Applico il modello ********************************************************************************************

d0 = 0.67 * h
z0 = 0.136 * h  # Brutsaert (1982)
pa = 101.3 * ((293 - 0.0065 * Alt) / 293) ** 5.26  # [Kpa]
es_kpa = es / 10  # [KPa]
ea_kpa = ea / 10  # [KPa]
ea_pa = ea * 100
ee = 0.622  # Water Vapor/Dry Air mix ratio
rv = (ee * es_kpa) / (pa - es_kpa)  # Mixing Ratio
Tv = TaK * (1 + (rv / ee)) / (1 + rv)  # Wallace, John M.; Hobbs, Peter V. (2006). Atmospheric Science. ISBN 0-12-732951-X.
ps = p_s * ((44331.0 - (Alt + z_pbl)) / (44331.0 - alt_ms)) ** (1.0 / 0.1903)
rhoam = (ps / (Rd * LST)) * (1.0 - 0.378 * ea_pa / ps)  # moist air density (kg m-3)
d = 2.0 / 3.0 * h  # zero plane displacement (m)

# FRICTION VELOCITY IN CONDITION OF STABILITY
uFi = (k * u) / np.log(2 / z0m)
nec = (Cd * LAI) / ((2 * (uFi ** 2)) / u ** 2)  # Within-canopy wind profile extinction coefficient
# Fs = 1-Fc  #Fracction of Soil Cover
kB = KB_1(u, uFi, z0m, Fc, LAI, TaK, pa, h)
show(kB)
Z0h = z0h(kB, z0m)
show(Z0h)
Z0h[Z0h == np.inf] = 28814991000.0  # 375579840.0 #high resistance
Z0h[Z0h > 3.093136e+10] = 28814991000.0
Z0h[np.isnan(Z0h)] = 1e-38
show(Z0h)
Upbl = u_pbl(u, 2, z0m, z_pbl)
upbl = Upbl[0]
show(upbl)
print("Calculating Dry Limit...")
H_d = np.copy(RN)
FRUstar_i = FRUstar(z_pbl, upbl, hst, z0m, Z0h, Alt, RH, TaK, LST)
RUstar = FRUstar_i[0]
show(RUstar)
RL = FRUstar_i[1]
show(RL)
kB = KB_1(u, RUstar, z0m, Fc, LAI, TaK, pa, h)
show(kB)
Z0h = z0h(kB, z0m)
show(Z0h)
print("Calculating Wet Limit...")
# For completely wet areas
# Wet-limit stability length
L_e = 2.430E+06  # Latent heat of vapourization (J kg-1) (Brutsaert, 1982)
L_w = (np.nanmin(RUstar) ** 3.0) * rhoam / (0.61 * k * g * RN / L_e)
# C_wet = ifthenelse(z_pbl >= hst, Cw(z_pbl, L_w, z0m, z0h), PSIh_y(pcrumin(z_pbl/L_w)))
# Initial conditional mask
Mask_Cwet1 = np.copy(z_pbl)
Mask_Cwet1[Mask_Cwet1 >= hst] = -9999.0
Mask_Cwet1[Mask_Cwet1 != 9999.0] = 0
Mask_Cwet1[Mask_Cwet1 == 9999.0] = 1
C_wet1 = Cw(z_pbl, L_w, z0m, Z0h)
C_wet1 = C_wet1 * Mask_Cwet1
C_wet1 = np.nan_to_num(C_wet1)
# Second conditional mask
Mask_Cwet2 = np.copy(z0)
Mask_Cwet2[Mask_Cwet2 >= hst] = -9999.0
Mask_Cwet2[Mask_Cwet2 != 9999.0] = 1
Mask_Cwet2[Mask_Cwet2 == 9999.0] = 0
C_wet2 = np.nanmin(PSIh_y(z_pbl / L_w))
C_wet2 = C_wet2 * Mask_Cwet2
C_wet2 = np.nan_to_num(C_wet2)
# SUM CONDITIONAL
C_wet = C_wet1 + C_wet2
zd0 = z_pbl - d
zdm = np.log(zd0 / z0m)
zdh = np.log(zd0 / Z0h)
zdh[zdh == np.inf] = 375579840.0  # high resistance
zdh[np.isnan(zdh)] = 375579840
# Wet-limit external resistance
re_w = (zdh - C_wet) / (k * RUstar)
# re_w = ifthenelse(re_w <= 0.0, zdh / (k * RUstar), re_w)
# Initial conditional mask
Mask_re_w1 = np.copy(re_w)
Mask_re_w1[Mask_re_w1 <= 0.0] = -9999.0
Mask_re_w1[Mask_re_w1 != -9999.0] = 0
Mask_re_w1[Mask_re_w1 == -9999.0] = 1
re_w1 = zdh / (k * RUstar)
re_w1[re_w1 == -np.inf] = 0
re_w1[re_w1 == np.inf] = 0
re_w1 = re_w1 * Mask_re_w1
re_w1 = np.nan_to_num(re_w1)
# Second conditional mask
Mask_re_w2 = np.copy(re_w)
Mask_re_w2[Mask_re_w2 <= 0.0] = -9999.0
Mask_re_w2[Mask_re_w2 != -9999.0] = 1
Mask_re_w2[Mask_re_w2 == -9999.0] = 0
re_w2 = re_w
re_w2[re_w2 == -np.inf] = 0
re_w2[re_w2 == np.inf] = 0
re_w2 = re_w2 * Mask_re_w2
re_w2 = np.nan_to_num(re_w2)
# SUM CONDITIONAL
rew = re_w1 + re_w2
rew[rew == 0] = np.nanmax(rew)
# Wet-limit heat flux
#slopef = 17.502 * 240.97 * ea / (TaK + 240.97) ** 2.0 #SEBS
slopef = (4098 * (0.6018 * np.exp((17.27 * (TaK - 273.15)) / ((TaK - 273.15) + 237.7))) / ((TaK - 273.15) + 237.7) ** 2) * 1000  # Slope of the saturation of vapor pressure [Pa/Â°C]: divided by 1000 get in Pascal
#slopef = (4098*ea/((TaK-273.15)+237.7)**2)
gamma = (((1.013 * 10 ** -3) * pa) / (ee * (L_e / 10 ** 6))) * 1000  # Psychrometric Constant [Pa/Â°C]:divided by 1000 to get in Pascal #FAO56
#gamma = 67.0   # psychrometric constant (Pa K-1)
rhoa = 1.23  # surface air density
rhoacp = rhoa * Cp
H_w = (RN - (rhoacp / rew) * ((ea - es) / gamma)) / (1.0 + slopef / gamma)
H_w[H_w < 0] = 0
LEwet = RN - H_w
# Sensible Heat flux
print("Calculating sensible heat flux...")
# C_i = ifthenelse(z_pbl >= hst, Cw(z_pbl, RL, z0m, z0h), PSIh_y(pcrumin(z_pbl)/RL))
# Initial conditional mask
RUstar[np.isnan(RUstar)] = np.nanmin(RUstar)
RL[np.isnan(RL)] = np.nanmin(RL)
# BAS MODELATION
Mask_C_i1 = np.copy(z_pbl)
Mask_C_i1[Mask_C_i1 >= hst] = -9999.0
Mask_C_i1[Mask_C_i1 != 9999.0] = 0
Mask_C_i1[Mask_C_i1 == 9999.0] = 1
C_i_1 = Cw(z_pbl, RL, z0m, Z0h)
C_i_1 = C_i_1 * Mask_C_i1
C_i_1 = np.nan_to_num(C_i_1)
# MOS MODELATION
# Second conditional mask
Mask_C_i2 = np.copy(z_pbl)
Mask_C_i2[Mask_C_i2 >= hst] = -9999.0
Mask_C_i2[Mask_C_i2 != 9999.0] = 1
Mask_C_i2[Mask_C_i2 == 9999.0] = 0
C_i_2 = PSIh_y(z_pbl / RL)
C_i_2 = C_i_2 * Mask_C_i2
C_i_2 = np.nan_to_num(C_i_2)
# SUM CONDITIONAL OF PIXELS OF MOS & BAS
C_i = C_i_1 + C_i_2
# external resistance
re_i = (zdh - C_i) / (k * RUstar)
"""T0Ta calculation"""
t_c = np.log((z_pbl - d) / Z0h) / np.log((alt_ms - d) / Z0h)
t_pbl_A = LST * (1.0 - t_c) + TaK * t_c
t_pbl_A = t_pbl_A / (1.0 - Alt / 44331.0) ** 1.5029
t_pbl = t_pbl_A
helpvar1 = Alt / 44331.0
helpvar2 = 1.0 - helpvar1
T0 = LST / helpvar2 ** 1.5029
z_pbl_A = 1000.0  # PBL height [m]
p_pbl_A = p_s * ((44331.0 - (Alt + z_pbl_A)) / (44331.0 - alt_ms)) ** (1.0 / 0.1903)
#Theta_a = t_pbl * (101325 / p_pbl_A) ** 0.286  # provo ad utilizzare questo Theta_a, ma mi vengono valori troppo elevati!
Theta_a = (TaK - 0.6 * 20) * (101325 / p_s) ** 0.286 #questo è quello di Italo ma non so perchè è così
show(Theta_a)
Theta_s = T0
T0ta = Theta_s - Theta_a
#T0ta = lst_b1 - TaK #questa è un'approssimazione
"""End calculation"""
H_i = rhoacp * (T0ta) / re_i
# H_i = ifthenelse(H_i > H_d, H_d, H_i)
print('Percentile 10 H: ' + str(np.nanpercentile(H_i, 10)))
print('Percentile 50 H: ' + str(np.nanpercentile(H_i, 50)))
print('Percentile 90 H: ' + str(np.nanpercentile(H_i, 90)))
# H_i = ifthenelse(H_i < H_w, H_w, H_i)
# report(H_i, hmap)
# Calculate evaporation variables
print("Calculating relative evaporation and evaporative fraction...")
# Calculate relative evaporation
# Ef =(EfR*LETwet)/(RN)
ev_r = 1.0 - (H_i - H_w) / (H_d - H_w)  # set water and wet surfaces to 1.0
# report(ev_r, evaprmap)
# Calculate evaporative fraction
ETF = ev_r * (1.0 - H_w / H_d)
# report(Evapfr, evapfrmap)
Lambda = 2501000 - (2361 * (TaK - 273.15))  # Latent heat of vaporization [J]
# Calculate latent energy flux
print("Calculating Latent Energy Flux...")
labdaE = ETF * (RN)  # Latent Heat Flux [W/m2] (Su, 2002)
labdaE2 = RN - H_i  # Latent Heat Flux [W/m2] (Su, 2002)
ET_Instant = RN - H_i - G0 #[W/m2]
# Evapfr = labdaE2/RN
ETR = (labdaE2 / Lambda) * Cdi * (24 * 3600)  # Actual Evapotranspiration [mm/day]
ETR[ETR > 12] = 12
ETR[ETR < 0] = 0
ETR = ETR
ETF[ETF > 1] = 1
ETF[ETF < 0] = 0
print('Actual Evapotranspiration average:' + str(np.nanmean(ETR)), 'mm/day')
print('Actual Evapotranspiration variability:' + str(np.nanstd(ETR)), 'mm/day')
show(ETR)


