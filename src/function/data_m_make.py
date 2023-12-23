import numpy as np
import healpy as hp
from numpy.linalg import solve as bslash
import math

# ストークスパラメータ Q, U のデータ m　配列を作成

def make_data_m(cmb, freq_band, which_model, nside):

    band = int(len(freq_band))

    N_pix = hp.nside2npix(nside)

    data_m = np.zeros(int(2*N_pix*band))

    for ii, freq in enumerate(freq_band):

        # file name
        dir = "../../make_map/fg_map_make/fg_map/"
        nside_name = "nside_"
        nside_n = str(nside)
        npy =".npy"
    
        fg_synch = "Synch_"
        GHz = "_GHz_"

        # Synchrotron name
        Synch_name = f"{dir}{fg_synch}{freq}{GHz}{nside_name}{nside_n}{npy}"

        fg_dust = "Dust_"
        GHz = "_GHz_"

        # Dust name
        Dust_name = f"{dir}{fg_dust}{freq}{GHz}{nside_name}{nside_n}{npy}"

        synch_data = np.load(Synch_name)
        
        dust_data = np.load(Dust_name)
        
        if which_model == "s1":
            
            m = synch_data + cmb
            
            data_m[int((2*ii)*N_pix) : int((2*ii+1)*N_pix)] = m[1]
            data_m[int((2*ii+1)*N_pix) : int(((2*ii+2))*N_pix)] = m[2]
            
        elif which_model == "d1":

            # data m[I, Q, U]
            m = dust_data + cmb

            data_m[int((2*ii)*N_pix) : int((2*ii+1)*N_pix)] = m[1]
            data_m[int((2*ii+1)*N_pix) : int(((2*ii+2))*N_pix)] = m[2]     

        elif which_model == "d1 and s1":

            # data m[I, Q, U]
            m = synch_data + dust_data + cmb

            data_m[int((2*ii)*N_pix) : int((2*ii+1)*N_pix)] = m[1]
            data_m[int((2*ii+1)*N_pix) : int(((2*ii+2))*N_pix)] = m[2]

        else:
            print("正しいモデル入れてね")

    return data_m.T

# InputとCleanのパワースペクトルを計算

def calc_spectrum(input_Q, input_U, clean_Q, clean_U, nside):
    
    N_pix = hp.nside2npix(nside)
    
    I = np.zeros(N_pix)
    
    Clean_map = [I, clean_Q, clean_U]
    
    Input_map = [I, input_Q, input_U]
    
    l = np.arange(0, 2*nside+1, 1)
    
    Clean_map_cl = hp.sphtfunc.anafast(Clean_map, lmax = 2*nside)
    
    Input_map_cl = hp.sphtfunc.anafast(Input_map, lmax = 2*nside)
    
    Clean_map_dl = Clean_map_cl * l * (l + 1) / (2 * math.pi)
    
    Input_map_dl = Input_map_cl * l * (l + 1) / (2 * math.pi)

    return Clean_map_cl, Clean_map_dl, Input_map_cl, Input_map_dl, l

# InputとCleanの差分のパワースペクトルを計算

def diff_calc_spectrum(input_Q, input_U, clean_Q, clean_U, nside):
    
    N_pix = hp.nside2npix(nside)
    
    I = np.zeros(N_pix)
    
    res_map = [I, clean_Q - input_Q, clean_U - input_U]
    
    l = np.arange(0, 3*nside, 1)
    
    res_map_cl = hp.sphtfunc.anafast(res_map)

    res_map_dl = res_map_cl * l * (l + 1) / (2 * math.pi)

    return res_map_cl, res_map_dl, l
    

# making input map for each frequency.

def make_input_map(cmb, freq_band, which_model, nside):

    band = int(len(freq_band))

    N_pix = hp.nside2npix(nside)

    Input_map = []

    if which_model == "s1":
        
        for ii, freq in enumerate(freq_band):

            # data m[I, Q, U]
            map = fg.fg_make_file(freq, nside)[0] + cmb

            Input_map.append(map)
            
    elif which_model == "d1":

        for ii in freq_band:

            # data m[I, Q, U]
            map = fg.fg_make_file(freq, nside)[1] + cmb

            Input_map.append(map)

    elif which_model == "d1 and s1":

        for ii in freq_band:

            # data m[I, Q, U]
            map = fg.fg_make_file(freq, nside)[0] + fg.fg_make_file(k, nside)[1] + cmb

            Input_map.append(map)

    return Input_map


def calc_clean_map(data_m, D, nside):
    
    #CMB_map = bslash((D.T @ bslash(N, D)), D.T @ bslash (N, m))
    # N + M +3 = Nfreq ===> just number of frequency bands (Synch + Dust)
    # N + 2 = Nfreq ===> just number of frequency bands (Dust)
    # M + 2 = Nfreq ===> just number of frequency bands (Synch)

    N_pix = hp.nside2npix(nside)
    
    clean_map = bslash(D, data_m)

    CMB_map_Q = clean_map[N_pix*0:N_pix*1]
    CMB_map_U = clean_map[N_pix*1:N_pix*2]

    return CMB_map_Q, CMB_map_U
