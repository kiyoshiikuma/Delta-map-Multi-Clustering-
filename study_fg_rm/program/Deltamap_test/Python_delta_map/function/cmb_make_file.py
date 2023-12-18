import healpy as hp
import numpy as np
import sys, platform, os
import math
import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))
# make sure the version and path is what you expect

def cmb_make_file(nside, r_input, random_seed_cmb, seed_syn):

    #Set up a new set of parameters for CAMB
    pars = camb.read_ini('/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/Python_delta_map/params_camb_for_PTEP.ini')
    #Set WantTensors to True
    pars.WantTensors = True
    #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
    pars.InitPower.set_params(As=2e-9, ns=0.965, r=1)
    #pars.set_for_lmax(3*nside - 1, lens_potential_accuracy=0)
    pars.set_for_lmax(2*nside, lens_potential_accuracy=0)
    pars.RandomSeed = random_seed_cmb;

    #calculate results for these parameters
    results = camb.get_results(pars)

    #get dictionary of CAMB power spectra
    #powers =results.get_cmb_power_spectra(pars, lmax = 3*nside - 1, CMB_unit='muK', raw_cl=True)
    powers =results.get_cmb_power_spectra(pars, lmax = 2*nside, CMB_unit='muK', raw_cl=True)
    #for name in powers: print(name)

    cl_scal = powers['unlensed_scalar']
    cl_tens = powers['tensor']
    cl_lens = powers['lensed_scalar']
    cl_pot = powers['lens_potential']
    cl_total = powers['total']

    cl = cl_scal + r_input * cl_tens + cl_lens

    #without lense
    #cl = cl_scal + r * cl_tens

    np.random.seed(seed_syn)
    cmb_map = hp.synfast(cl.T, nside, new = True)

    return cmb_map


def cmb_cell_make(nside, random_seed_cmb, seed_syn, r_input = 1):

    #Set up a new set of parameters for CAMB
    pars = camb.read_ini('/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/Python_delta_map/params_camb_for_PTEP.ini')
    #Set WantTensors to True
    pars.WantTensors = True
    #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
    pars.InitPower.set_params(As=2e-9, ns=0.965, r=1)
    pars.set_for_lmax(2*nside, lens_potential_accuracy=0)
    pars.RandomSeed = random_seed_cmb;

    #calculate results for these parameters
    results = camb.get_results(pars)

    #get dictionary of CAMB power spectra
    powers =results.get_cmb_power_spectra(pars, lmax = 2*nside, CMB_unit='muK', raw_cl=True)
    #for name in powers: print(name)

    cl_scal = powers['unlensed_scalar']
    cl_tens = powers['tensor']
    cl_lens = powers['lensed_scalar']
    cl_pot = powers['lens_potential']
    cl_total = powers['total']

    #without lense
    #cl = cl_scal + r * cl_tens

    return cl_scal, cl_lens, cl_tens
    
    