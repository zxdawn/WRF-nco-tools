import argparse
import numpy as np
from numpy import exp
from netCDF4 import Dataset
from wrf import getvar, to_np

def parse_args():
    parser = argparse.ArgumentParser(description="program to calculate CloudFraction")
    parser.add_argument("filename", help="name of wrfout file")
    args = parser.parse_args()
    return vars(args)

def cal_cldfra1(filename):
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
    ep_2=r_d/r_v

    #---------------------------------------------------------------------

    ncfile = Dataset(filename)

    t_phy  = to_np(getvar(ncfile, 'temp'))
    p_phy  = to_np(getvar(ncfile, 'p'))
    QV     = to_np(getvar(ncfile, 'QVAPOR')) # water vapor mixing ratio (kg/kg)
    QC     = to_np(getvar(ncfile, 'QCLOUD')) # cloud water mixing ratio (kg/kg)
    Qr     = to_np(getvar(ncfile, 'QRAIN'))  # rain water mixing ratio (kg/kg)
    QI     = to_np(getvar(ncfile, 'QICE'))   # cloud ice mixing ratio (kg/kg)
    QS     = to_np(getvar(ncfile, 'QSNOW'))  # snow mixing ratio (kg/kg)

    CLDFRA = np.zeros_like(t_phy)

    for i in np.arange(t_phy.shape[0]):
        for k in np.arange(t_phy.shape[1]):
            for j in np.arange(t_phy.shape[2]):
                tc   = t_phy[i,k,j] - SVPT0
                esw  = 1000.0 * SVP1 * exp( SVP2  * tc / ( t_phy[i,k,j] - SVP3  ) )
                esi  = 1000.0 * SVP1 * exp( SVPI2 * tc / ( t_phy[i,k,j] - SVPI3 ) )
                QVSW = ep_2 * esw / ( p_phy[i,k,j] - esw )
                QVSI = ep_2 * esi / ( p_phy[i,k,j] - esi )

                # mji - For MP options 2, 4, 6, 7, 8, etc. (qc = liquid, qi = ice, qs = snow)            
                if (ncfile.MP_PHYSICS in (2, 4, 6, 7, 8)):
                    QCLD = QI[i,k,j]+QC[i,k,j]+QS[i,k,j]
                    if (QCLD < QCLDMIN):
                        weight = 0.
                    else:
                        weight = (QI[i,k,j]+QS[i,k,j]) / QCLD

                # mji - For MP options 1 and 3, (qc only)
                #  For MP=1, qc = liquid, for MP=3, qc = liquid or ice depending on temperature
                elif (ncfile.MP_PHYSICS in (1, 3)):
                    QCLD = QC[i,k,j]
                    if (QCLD < QCLDMIN):
                       weight = 0.
                    else:
                       if (t_phy[i,k,j] > 273.15): weight = 0.
                       if (t_phy[i,k,j] <= 273.15): weight = 1.        

                # mji - For MP option 5; (qc = liquid, qs = ice)
                # Mixing ratios of cloud water & total ice (cloud ice + snow).
                # Mixing ratios of rain are not considered in this scheme.
                # F_ICE is fraction of ice
                # F_RAIN is fraction of rain
                elif (ncfile.MP_PHYSICS in (5)):
                    QIMID = QS[i,k,j]
                    QWMID = QC[i,k,j]
                    # old method
                    #           QIMID = QC[i,k,j]*F_ICE_PHY[i,k,j]
                    #           QWMID = (QC[i,k,j]-QIMID)*(1.-F_RAIN_PHY[i,k,j])
                    #
                    #--- Total "cloud" mixing ratio, QCLD.  Rain is not part of cloud,
                    #    only cloud water + cloud ice + snow
                    #
                    QCLD = QWMID + QIMID
                    if (QCLD < QCLDMIN):
                        weight = 0.
                    else:
                        weight = F_ICE_PHY[i,k,j]

                else:
                    CLDFRA[i,k,j] = 0.

                QVS_WEIGHT = (1-weight)*QVSW + weight*QVSI
                RHUM = QV[i,k,j]/QVS_WEIGHT   #--- Relative humidity
                #
                #--- Determine cloud fraction (modified from original algorithm)
                #
                if (QCLD < QCLDMIN):
                #
                #--- Assume zero cloud fraction if there is no cloud mixing ratio
                #
                    CLDFRA[i,k,j] = 0.
                elif (RHUM > RHGRID):
                #
                #--- Assume cloud fraction of unity if near saturation and the cloud
                #    mixing ratio is at or above the minimum threshold
                #
                    CLDFRA[i,k,j] = 1.
                else:
                #
                #--- Adaptation of original algorithm (Randall, 1994; Zhao, 1995)
                #    modified based on assumed grid-scale saturation at RH=RHgrid.
                #
                    SUBSAT = max(1.e-10,RHGRID*QVS_WEIGHT-QV[i,k,j])
                    DENOM  = (SUBSAT)**GAMMA
                    ARG    = max(-6.9, -ALPHA0*QCLD/DENOM)    # <-- exp(-6.9)=.001
                # prevent negative values  (new)
                    RHUM   = max(1.e-10, RHUM)
                    CLDFRA[i,k,j] = (RHUM/RHGRID)**PEXP*(1.-exp(ARG))
                ##              ARG=-1000*QCLD/(RHUM-RHGRID)
                ##              ARG=max(ARG, ARGMIN)
                ##              CLDFRA[i,k,j]=(RHUM/RHGRID)*(1.-exp(ARG))
                    if (CLDFRA[i,k,j] < .01):
                        CLDFRA[i,k,j] = 0.
    return CLDFRA

def create(filename):
    ncfile = Dataset(filename,'r+')

    # Create variable: CloudFraction (0~1)
    cldfra_l           = ncfile.createDimension('cldfra_l', 1)
    cldfra             = ncfile.createVariable('cldfra','f4',('Time', 'cldfra_l', 'south_north', 'west_east'))
    cldfra.description = 'CloudFraction generated from wrf'
    cldfra.FieldType   = '104'
    cldfra.MemoryOrder = 'XYZ'
    cldfra.units       = ''
    cldfra.stagger     = ''
    cldfra.coordinates = 'XLONG XLAT'
    return (cldfra)

def cldfra_max(filename,cldfra):
    # Calculate cldfra_max between 350 hPa and 400 hPa
    ncfile   = Dataset(filename,'r+')
    P        = to_np(getvar(ncfile, 'P'))
    PB       = to_np(getvar(ncfile, 'PB'))
    pressure = (P + PB)/100 # hPa

    CLDFRA = cal_cldfra1(filename)
    for j in np.arange(pressure.shape[1]):
        for k in np.arange(pressure.shape[2]):
            ncfile.variables['cldfra'][:,:,j,k] = np.max(CLDFRA[np.where( (350 < pressure[:,j,k]) & (pressure[:,j,k] < 400))[0],j,k])

def main(filename):
    cldfra = create(filename)
    cldfra_max(filename,cldfra)

if __name__ == '__main__':
    args = parse_args()
    main(**args)
