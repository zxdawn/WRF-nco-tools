import sys
import itertools
import numpy as np
from numpy import exp
from numba import njit
from netCDF4 import Dataset
import multiprocessing as mp
from functools import partial
from wrf import getvar, to_np


def get_variables(filename):
    global f
    f = Dataset(filename,'r+')
    phys   = f.MP_PHYSICS
    P      = to_np(getvar(f, 'P'))
    PB     = to_np(getvar(f, 'PB'))
    t_phy  = to_np(getvar(f, 'temp'))
    p_phy  = to_np(getvar(f, 'p'))
    QV     = to_np(getvar(f, 'QVAPOR')) # water vapor mixing ratio (kg/kg)
    QC     = to_np(getvar(f, 'QCLOUD')) # cloud water mixing ratio (kg/kg)
    Qr     = to_np(getvar(f, 'QRAIN'))  # rain water mixing ratio (kg/kg)
    QI     = to_np(getvar(f, 'QICE'))   # cloud ice mixing ratio (kg/kg)
    QS     = to_np(getvar(f, 'QSNOW'))  # snow mixing ratio (kg/kg)
    return phys, P, PB, t_phy, p_phy, QV, QC, Qr, QI, QS


def create():
    # Create variable: CloudFraction (0~1)
    cldfra             = f.createVariable('cldfra','f4',('Time', 'south_north', 'west_east'))
    cldfra.description = 'CloudFraction generated from wrf'
    cldfra.FieldType   = '104'
    cldfra.MemoryOrder = 'XYZ'
    cldfra.units       = ''
    cldfra.stagger     = ''
    cldfra.coordinates = 'XLONG XLAT'

    return cldfra


def cal_cldfra(phys, P, PB, t_phy, p_phy, QV, QC, Qr, QI, QS):
    # Calculate cldfra_max between 350 hPa and 400 hPa
    pressure = (P + PB)/100 # hPa

    CLDFRA = np.zeros_like(t_phy)
    CLDFRA = cal_cldfra1(phys, t_phy, p_phy, QV, QC, Qr, QI, QS, CLDFRA)

    temp_array = np.zeros((1,pressure.shape[1],pressure.shape[2]))
    CLDFRA = cldfra_max(pressure, CLDFRA, temp_array)
    f.variables['cldfra'][:] = CLDFRA


@njit
def cal_cldfra1(phys, t_phy, p_phy, QV, QC, Qr, QI, QS, CLDFRA):
    # cal_cldfra1 - Compute cloud fraction
    # Code adapted from that in module_ra_gfdleta.F in WRF_v2.0.3 by James Done
    ##
    ##---  Cloud fraction parameterization follows Xu and Randall (JAS), 1996
    ##     (see Hong et al., 1998)
    ##     (modified by Ferrier, Feb '02)
    ##     (converted to py by Xin, Jun 2018) only tested for MP options 2

    # DESCRIPTION:
    # Compute cloud fraction from input ice and cloud water fields
    # if provided.
    #
    # Whether QI or QC is active or not is determined from the indices of
    # the fields into the 4D scalar arrays in WRF. These indices are 
    # P_QI and P_QC, respectively, and they are passed in to the routine
    # to enable testing to see if QI and QC represent active fields in
    # the moisture 4D scalar array carried by WRF.
    # 
    # If a field is active its index will have a value greater than or
    # equal to PARAM_FIRST_SCALAR, which is also an input argument to 
    # this routine.


    #-----------------------------------------------------------------------
    #---  COMPUTE GRID-SCALE CLOUD COVER FOR RADIATION
    #     (modified by Ferrier, Feb '02)
    #
    #---  Cloud fraction parameterization follows Randall, 1994
    #     (see Hong et al., 1998)
    #-----------------------------------------------------------------------
    # Note: ep_2=287./461.6 Rd/Rv
    # Note: R_D=287.

    # Alternative calculation for critical RH for grid saturation
    #     RHGRID=0.90+.08*((100.-DX)/95.)**.5

    # Calculate saturation mixing ratio weighted according to the fractions of
    # water and ice.
    # Following:
    # Murray, F.W. 1966. ``On the computation of Saturation Vapor Pressure''  J. Appl. Meteor.  6 p.204
    #    es (in mb) = 6.1078 . exp[ a . (T-273.16)/ (T-b) ]
    #
    #       over ice        over water
    # a =   21.8745584      17.2693882
    # b =   7.66            35.86

    for i in np.arange(t_phy.shape[0]):
        for j in np.arange(t_phy.shape[1]):
            for k in np.arange(t_phy.shape[2]):
                tc   = t_phy[i,j,k] - SVPT0
                esw  = 1000.0 * SVP1 * exp( SVP2  * tc / ( t_phy[i,j,k] - SVP3  ) )
                esi  = 1000.0 * SVP1 * exp( SVPI2 * tc / ( t_phy[i,j,k] - SVPI3 ) )
                QVSW = ep_2 * esw / ( p_phy[i,j,k] - esw )
                QVSI = ep_2 * esi / ( p_phy[i,j,k] - esi )

                # mji - For MP options 2, 4, 6, 7, 8, etc. (qc = liquid, qi = ice, qs = snow)            
                if (phys in (2, 4, 6, 7, 8)):
                    QCLD = QI[i,j,k]+QC[i,j,k]+QS[i,j,k]
                    if (QCLD < QCLDMIN):
                        weight = 0.
                    else:
                        weight = (QI[i,j,k]+QS[i,j,k]) / QCLD

                else:
                    CLDFRA[i,j,k] = 0.

                QVS_WEIGHT = (1-weight)*QVSW + weight*QVSI
                RHUM = QV[i,j,k]/QVS_WEIGHT   #--- Relative humidity
                #
                #--- Determine cloud fraction (modified from original algorithm)
                #
                if (QCLD < QCLDMIN):
                #
                #--- Assume zero cloud fraction if there is no cloud mixing ratio
                #
                    CLDFRA[i,j,k] = 0.
                elif (RHUM > RHGRID):
                #
                #--- Assume cloud fraction of unity if near saturation and the cloud
                #    mixing ratio is at or above the minimum threshold
                #
                    CLDFRA[i,j,k] = 1.
                else:
                #
                #--- Adaptation of original algorithm (Randall, 1994; Zhao, 1995)
                #    modified based on assumed grid-scale saturation at RH=RHgrid.
                #
                    SUBSAT = max(1.e-10,RHGRID*QVS_WEIGHT-QV[i,j,k])
                    DENOM  = (SUBSAT)**GAMMA
                    ARG    = max(-6.9, -ALPHA0*QCLD/DENOM)    # <-- exp(-6.9)=.001
                # prevent negative values  (new)
                    RHUM   = max(1.e-10, RHUM)
                    CLDFRA[i,j,k] = (RHUM/RHGRID)**PEXP*(1.-exp(ARG))
                ##              ARG=-1000*QCLD/(RHUM-RHGRID)
                ##              ARG=max(ARG, ARGMIN)
                ##              CLDFRA[i,j,k]=(RHUM/RHGRID)*(1.-exp(ARG))
                    if (CLDFRA[i,j,k] < .01):
                        CLDFRA[i,j,k] = 0.

    return CLDFRA


@njit
def cldfra_max(pressure, CLDFRA, temp_array):
    for j in np.arange(pressure.shape[1]):
        for k in np.arange(pressure.shape[2]):
            temp_array[:,j,k] = np.max(CLDFRA[np.where( (350 < pressure[:,j,k]) & (pressure[:,j,k] < 400))[0],j,k])

    return temp_array


def main(filename):
    phys, P, PB, t_phy, p_phy, QV, QC, Qr, QI, QS = get_variables(filename)
    cldfra = create()
    cal_cldfra(phys, P, PB, t_phy, p_phy, QV, QC, Qr, QI, QS)
    f.close()


if __name__ == '__main__':
    #---------------------------------------------------------------------
    # Constants
    ALPHA0 = 100.
    GAMMA = 0.49
    QCLDMIN = 1.e-12
    PEXP = 0.25
    RHGRID = 1.0
    SVP1 = 0.61078
    SVP2 = 17.2693882
    SVPI2 = 21.8745584
    SVP3 = 35.86
    SVPI3 = 7.66
    SVPT0 = 273.15
    r_d = 287.
    r_v = 461.6
    ep_2 = r_d/r_v

    #---------------------------------------------------------------------
    with open('read_wrf.conf') as name:
        namelist = name.read().replace('\n', '').split()
        namelist = [sys.argv[1]+s for s in namelist]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(main, namelist)