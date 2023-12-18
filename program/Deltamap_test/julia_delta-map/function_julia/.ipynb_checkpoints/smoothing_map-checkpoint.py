import healpy as hp
import math

def smootinhg_map(input_map, nside):
    
    alm = hp.sphtfunc.map2alm(input_map, lmax = 2*nside)

    input_fwhm = 2200 * (4/nside) ** 2
    
    smoothed_map = hp.sphtfunc.alm2map(alm, nside, lmax=2*nside, pixwin=True, verbose=False, fwhm=input_fwhm*math.pi/10800.)

    return smoothed_map

    
def smootinhg_map_fwhm(input_map, input_fwhm, nside):
    
    alm = hp.sphtfunc.map2alm(input_map, lmax = 2*nside)
    
    smoothed_map = hp.sphtfunc.alm2map(alm, nside, lmax=2*nside, pixwin=True, verbose=False, fwhm=input_fwhm*math.pi/10800.)

    return smoothed_map
