import numpy as np
import healpy as hp
import fg_make_file as fg
from numpy.linalg import solve as bslash
import math

# ストークスパラメータ Q, U のデータ m　配列を作成

def make_data_m(cmb, freq_band, which_model, nside):

    band = int(len(freq_band))

    N_pix = hp.nside2npix(nside)

    data_m = np.zeros(int(2*N_pix*band))

    if which_model == "s1":
        
        for ii, freq in enumerate(freq_band):

            # data m[I, Q, U]
            m = fg.fg_make_file(freq, nside)[0] + cmb
            
            data_m[int((2*ii)*N_pix) : int((2*ii+1)*N_pix)] = m[1]
            data_m[int((2*ii+1)*N_pix) : int(((2*ii+2))*N_pix)] = m[2]
            
    elif which_model == "d1":

        for ii, freq in enumerate(freq_band):

            # data m[I, Q, U]
            m = fg.fg_make_file(freq, nside)[1] + cmb

            data_m[int((2*ii)*N_pix) : int((2*ii+1)*N_pix)] = m[1]
            data_m[int((2*ii+1)*N_pix) : int(((2*ii+2))*N_pix)] = m[2]     

    elif which_model == "d1 and s1":

        for ii, freq in enumerate(freq_band):

            # data m[I, Q, U]
            m = fg.fg_make_file(freq, nside)[0] + fg.fg_make_file(freq, nside)[1] + cmb

            data_m[int((2*ii)*N_pix) : int((2*ii+1)*N_pix)] = m[1]
            data_m[int((2*ii+1)*N_pix) : int(((2*ii+2))*N_pix)] = m[2]

    else:
        print("正しいモデル入れてね")

    return data_m
    
# Dマトリクスの計算

def D_element(freq, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,):
    
    
    # x = (f / T) * (h / k)
    x = (freq / 2.725) / (2.083661912 * 10**10)
    
    g_freq = ((np.exp(x) - 1)**2) / (np.exp(x) * x**2) * 1000.0
    
    s = g_freq * (freq/freq_bs)**(beta_s)
    
    ss = g_freq * (freq/freq_bs)**(beta_s) * np.log(freq/freq_bs)

    x_d = (freq / T_d) / (2.083661912 * 10**10)

    x_bd = (freq_bd / T_d) / (2.083661912 * 10**10)

    d = g_freq * (freq/freq_bd)**(beta_d + 1) * ((np.exp(x_bd)-1)/(np.exp(x_d)-1))

    dd = d * np.log(freq/freq_bd)

    ddd = d * (((x_d*np.exp(x_d))/(np.exp(x_d)-1)) - (x_bd*np.exp(x_bd))/(np.exp(x_bd)-1))/T_d
    
    return s, ss, d, dd, ddd



def  calc_D_matrix(freq_band, which_model, nside):

    N_pix = hp.nside2npix(nside)

    diag = np.identity(2 * N_pix)

    D_blocks = []

    for freq in freq_band:
        
        s, ss, d, dd, ddd = D_element(freq * 10**9)
        
        if which_model == "s1":
        
            D_blocks.append([diag, s * diag, ss * diag])
            
            D_matrix = np.block(D_blocks)
        
        elif which_model == "d1":

            D_blocks.append([diag, d * diag, dd * diag, ddd * diag])
            
            D_matrix = np.block(D_blocks)

        elif which_model == "d1 and s1":

            D_blocks.append([diag, s * diag, ss * diag, d * diag, dd * diag, ddd * diag])
            
            D_matrix = np.block(D_blocks)

    return D_matrix


# parameter input ver

def  calc_D_matrix_param(freq_band, which_model, nside, freq_bs, beta_s, freq_bd, beta_d, T_d):

    N_pix = hp.nside2npix(nside)

    diag = np.identity(2 * N_pix)

    D_blocks = []

    for freq in freq_band:
        
        s, ss, d, dd, ddd = D_element(freq * 10**9, freq_bs, beta_s, freq_bd, beta_d, T_d)
        
        if which_model == "s1":
        
            D_blocks.append([diag, s * diag, ss * diag])
            
            D_matrix = np.block(D_blocks)
        
        elif which_model == "d1":

            D_blocks.append([diag, d * diag, dd * diag, ddd * diag])
            
            D_matrix = np.block(D_blocks)

        elif which_model == "d1 and s1":

            D_blocks.append([diag, s * diag, ss * diag, d * diag, dd * diag, ddd * diag])
            
            D_matrix = np.block(D_blocks)

    return D_matrix
    
    
# クリーンマプの計算

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


# InputとCleanのパワースペクトルを計算

def calc_spectrum(input_Q, input_U, clean_Q, clean_U, nside):
    
    N_pix = hp.nside2npix(nside)
    
    I = np.zeros(N_pix)
    
    Clean_map = [I, clean_Q, clean_U]
    
    Input_map = [I, input_Q, input_U]
    
    l = np.arange(0, 3*nside, 1)
    
    Clean_map_cl = hp.sphtfunc.anafast(Clean_map)
    
    Input_map_cl = hp.sphtfunc.anafast(Input_map)
    
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

