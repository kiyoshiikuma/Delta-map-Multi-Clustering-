import sys
sys.path.append('/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/Python_delta_map/function')

import numpy as np
import healpy as hp
import quaternionic
import spherical
import fg_make_file as fg
from numpy.linalg import solve as bslash


# 前景放射要素計算　角周波数

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
    
    return s, ss, d, dd, ddd, g_freq


def d_s_vec_calc(freq, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,):

    s, ss, d, dd, ddd, g_freq = D_element(freq * 10**9, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,)

    s_vec = np.array([s, ss]).reshape(-1, 1)

    d_vec = np.array([d, dd, ddd]).reshape(-1, 1)

    return s_vec, d_vec, g_freq


def Amat_calc(freq_band, cmb_freq, which_model, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,):

    A_s = []
    A_d = []
    
    for freq in freq_band:
        
        s_vec, d_vec, g_freq = d_s_vec_calc(freq, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,)

        if cmb_freq != freq:

            if which_model == "s1":

                A_s.append(s_vec)
            
                #print(s_vec)
        
            elif which_model == "d1":
            
                A_d.append(d_vec)

                #print(d_vec)

            elif which_model == "d1 and s1":

                A_s.append(s_vec)
                A_d.append(d_vec)

        else:
            pass

    if which_model == "s1":
        
        A = np.concatenate(A_s).reshape(-1, 2)

    elif which_model == "d1":

        A = np.concatenate(A_d).reshape(-1, 3)

    elif which_model == "d1 and s1":

        A_d = np.concatenate(A_d).reshape(-1, 3)
        A_s = np.concatenate(A_s).reshape(-1, 2)

        #print(A_d.T)
        #print(A_s.T)

        #A = np.vstack([A_d.T, A_s.T])
        A = np.hstack([A_d, A_s])

    else:
        print("正しいモデルを入力してね")
    
    return A.T

def calc_x(freq_band, cmb_freq, nside, cmb_map, which_model, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,):

    data_m = make_input_map(cmb_map, freq_band, which_model, nside)
    Amat = Amat_calc(freq_band, cmb_freq, which_model, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,)
    s_cmb, d_cmb, g_cmb = d_s_vec_calc(cmb_freq, freq_bs = 23*10**9, beta_s = -3., freq_bd = 353*10**9, beta_d = 1.5, T_d =20.1,)

    if which_model == "s1":
        vec_cmb = np.vstack([s_cmb])

    elif which_model == "d1":
        vec_cmb = np.vstack([d_cmb])

    elif which_model == "d1 and s1":
        vec_cmb = np.vstack([d_cmb, s_cmb])

    from numpy.linalg import solve as bslash

    # N + M +3 = Nfreq ===> just number of frequency bands (Synch + Dust)
    # N + 2 = Nfreq ===> just number of frequency bands (Dust)
    # M + 2 = Nfreq ===> just number of frequency bands (Synch)

    alpha_i = - bslash(Amat, vec_cmb)

    Q = 0 * data_m[0][1]
    U = 0 * data_m[0][1]

    #ここのsumと1の要素を入れる順は間違えないで
    sum_alpha = sum(alpha_i)

    freq_list = np.array(freq_band)
    cmb_index = list(np.where(freq_list == cmb_freq)[0])

    alpha_i = alpha_i.T
    alpha_i = alpha_i.tolist()
    alpha_i = alpha_i[0]
    alpha_i.insert(int(cmb_index[0]), 1)
    #print(alpha_i)
    
    for i, freq  in enumerate(freq_band):
        
        Q += alpha_i[i]*data_m[i][1]
        U += alpha_i[i]*data_m[i][2]

    Q = Q / (1 + sum_alpha)
    U = U / (1 + sum_alpha)

    x = np.concatenate([Q, U], 0)
    
    return Q, U, x

#=========================================================================================================#

#==================================フリーパラメータ============================================#



def d_s_vec_calc_pra(freq, freq_bs, beta_s, freq_bd, beta_d, T_d):

    s, ss, d, dd, ddd, g_freq = D_element(freq * 10**9, freq_bs, beta_s, freq_bd, beta_d, T_d)

    s_vec = np.array([s, ss]).reshape(-1, 1)

    d_vec = np.array([d, dd, ddd]).reshape(-1, 1)

    return s_vec, d_vec, g_freq


def Amat_calc_pra(freq_band, cmb_freq, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d):

    A_s = []
    A_d = []
    
    for freq in freq_band:
        
        s_vec, d_vec, g_freq = d_s_vec_calc_pra(freq, freq_bs, beta_s, freq_bd, beta_d, T_d)

        if cmb_freq != freq:

            if which_model == "s1":

                A_s.append(s_vec)
            
                #print(s_vec)
        
            elif which_model == "d1":
            
                A_d.append(d_vec)

                #print(d_vec)

            elif which_model == "d1 and s1":

                A_s.append(s_vec)
                A_d.append(d_vec)

        else:
            pass

    if which_model == "s1":
        
        A = np.concatenate(A_s).reshape(-1, 2)

    elif which_model == "d1":

        A = np.concatenate(A_d).reshape(-1, 3)

    elif which_model == "d1 and s1":

        A_d = np.concatenate(A_d).reshape(-1, 3)
        A_s = np.concatenate(A_s).reshape(-1, 2)

        #print(A_d.T)
        #print(A_s.T)

        #A = np.vstack([A_d.T, A_s.T])
        A = np.hstack([A_d, A_s])

    else:
        print("正しいモデルを入力してね")
    
    return A.T

def calc_x_pra(freq_band, cmb_freq, nside, cmb_map, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d):

    data_m = make_input_map(cmb_map, freq_band, which_model, nside)
    Amat = Amat_calc_pra(freq_band, cmb_freq, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)
    s_cmb, d_cmb, g_cmb = d_s_vec_calc_pra(cmb_freq, freq_bs, beta_s, freq_bd, beta_d, T_d)

    if which_model == "s1":
        vec_cmb = np.vstack([s_cmb])

    elif which_model == "d1":
        vec_cmb = np.vstack([d_cmb])

    elif which_model == "d1 and s1":
        vec_cmb = np.vstack([d_cmb, s_cmb])

    from numpy.linalg import solve as bslash

    # N + M +3 = Nfreq ===> just number of frequency bands (Synch + Dust)
    # N + 2 = Nfreq ===> just number of frequency bands (Dust)
    # M + 2 = Nfreq ===> just number of frequency bands (Synch)

    alpha_i = - bslash(Amat, vec_cmb)

    Q = 0 * data_m[0][1]
    U = 0 * data_m[0][1]

    #ここのsumと1の要素を入れる順は間違えないで
    sum_alpha = sum(alpha_i)

    freq_list = np.array(freq_band)
    cmb_index = list(np.where(freq_list == cmb_freq)[0])

    alpha_i = alpha_i.T
    alpha_i = alpha_i.tolist()
    alpha_i = alpha_i[0]
    alpha_i.insert(int(cmb_index[0]), 1)
    #print(alpha_i)
    
    for i, freq  in enumerate(freq_band):
        
        Q += alpha_i[i]*data_m[i][1]
        U += alpha_i[i]*data_m[i][2]

    Q = Q / (1 + sum_alpha)
    U = U / (1 + sum_alpha)

    x = np.concatenate([Q, U], 0)
    
    return Q, U, x

    
#=========================================================================================================#



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

        for freq in freq_band:

            # data m[I, Q, U]
            map = fg.fg_make_file(freq, nside)[1] + cmb

            Input_map.append(map)

    elif which_model == "d1 and s1":

        for freq in freq_band:

            # data m[I, Q, U]
            map = fg.fg_make_file(freq, nside)[0] + fg.fg_make_file(freq, nside)[1] + cmb

            Input_map.append(map)

    return Input_map

#=========球面調和関数計算====================#

def Cal_sYlm(nside):
    #npix = 12*nside**2
    npix = hp.nside2npix(nside) 
    lmax = 2*nside
    Y_all = []
    Y__all = []
    W_all = []
    X_all = []
    
    for pix in range(npix):
        # thita, phi
        theta, phi = hp.pixelfunc.pix2ang(nside, pix, nest=False)
        
        # quaternionic R
        R = quaternionic.array.from_spherical_coordinates(theta, phi)
        
        # wigner
        wigner = spherical.Wigner(lmax)
        
        Y = wigner.sYlm(2, R)
        Y_ = wigner.sYlm(-2, R)
        W = -1*(Y+Y_)/2
        X = -1j*(Y-Y_)/2
        #Y_all.append(Y)
        #Y__all.append(Y_)
        W_all.append(W)
        X_all.append(X)
    return np.array(W_all), np.array(X_all)


def Cal_matrix(F1lm, F2lm, G1lm, G2lm, clEE, clBB):
    npix = len(F1lm)
    nside = hp.pixelfunc.npix2nside(npix)
    lmax = 2*nside

    # beam, pixwinも考えない場合(= 1) theta_pixも適当---------------------------------------#
    theta_pix = 1
    wl = hp.sphtfunc.pixwin(nside, pol=True, lmax=lmax)
    wl_P = wl[1] * 0 + 1
    beam = hp.sphtfunc.gauss_beam(theta_pix, pol=True, lmax=lmax) * 0 + 1
    beam_EE = beam[:,1] * 0 + 1
    beam_BB = beam[:,2] * 0 + 1
    #----------------------------------------------------------------------------#

    # nestへの変換はコメントアウト
    matrix = np.zeros((npix,npix))*1.0j
    for i in range(npix):
        #i_nest = hp.pixelfunc.ring2nest(nside,i)
        for j in range(npix):
            #j_nest = hp.pixelfunc.ring2nest(nside,j)
            comp = Cal_sum_l(F1lm, F2lm, i, j, clEE, wl_P, beam_EE, nside) + Cal_sum_l(G1lm, G2lm, i, j, clBB, wl_P, beam_BB, nside)
            real = float(format(comp.real, '.13f'))
            imag = float(format(comp.imag, '.13f'))
            #print(real, imag)
            
            # nestへの変換はコメントアウト
            #matrix[i_nest,j_nest] = complex(real, imag)
            matrix[i,j] = complex(real, imag)
            
    return matrix
    

def Cal_sum_m(Flm, Glm, npix1, npix2, l, nside):
    Wl = []
    lmax = 2*nside
    # wigner
    wigner = spherical.Wigner(lmax)
    for m in range(-l,l+1):
        f = Flm[npix1][wigner.Yindex(l,m)]*(Glm[npix2][wigner.Yindex(l,m)].conjugate())
        Wl.append(f)
    return sum(Wl)


def Cal_sum_l(Flm, Glm, npix1, npix2, cl, wl, beam_pol, nside):
    mat_comp = []
    lmax = 2*nside
    for l in range(2,lmax+1):
        #fact = l*(l+1)/(2*np.pi)
        fact = 1
        Wl = Cal_sum_m(Flm, Glm, npix1, npix2, l, nside)
        f = wl[l]**2*beam_pol[l]**2*cl[l]/fact*Wl
        mat_comp.append(f)
    return sum(mat_comp)


def Cal_cov(Wlm, Xlm, clEE, clBB):
    npix = len(Wlm)
    C_QQ = Cal_matrix(Wlm, Wlm, Xlm, Xlm, clEE, clBB)
    C_QU = Cal_matrix(-Wlm, Xlm, Xlm, Wlm, clEE, clBB)
    C_UQ = Cal_matrix(-Xlm, Wlm, Wlm, Xlm, clEE, clBB)
    C_UU = Cal_matrix(Xlm, Xlm, Wlm, Wlm, clEE, clBB)
    
    C = np.zeros((2*npix,2*npix))
    C[0:npix, 0:npix] = C_QQ
    C[0:npix, npix:2*npix] = C_QU
    C[npix:2*npix, 0:npix] = C_UQ
    C[npix:2*npix, npix:2*npix] = C_UU
    
    return C

def Cal_cov_cmb(clEE, clBB, nside):

    Wlm, Xlm = Cal_sYlm(nside) 
    npix = len(Wlm)
    C_QQ = Cal_matrix(Wlm, Wlm, Xlm, Xlm, clEE, clBB)
    C_QU = Cal_matrix(-Wlm, Xlm, Xlm, Wlm, clEE, clBB)
    C_UQ = Cal_matrix(-Xlm, Wlm, Wlm, Xlm, clEE, clBB)
    C_UU = Cal_matrix(Xlm, Xlm, Wlm, Wlm, clEE, clBB)
    
    C = np.zeros((2*npix,2*npix))
    C[0:npix, 0:npix] = C_QQ
    C[0:npix, npix:2*npix] = C_QU
    C[npix:2*npix, 0:npix] = C_UQ
    C[npix:2*npix, npix:2*npix] = C_UU
    
    return C


def Cal_cov_mat(cov_mat_scal, cov_mat_tens, r):

    cov_mat = cov_mat_scal + r * cov_mat_tens
    
    return cov_mat


def log_det(A):
    try:
        logdet = np.sum(np.log(np.linalg.eigvalsh(A)))
        return logdet
    except np.linalg.LinAlgError:
        raise ValueError("Matrix is numerically singular!")

def smootinhg_map(input_map, input_fwhm, nside):
    
    alm = hp.sphtfunc.map2alm(input_map, lmax = 2*nside)

    fwhm = 2200 * (4/nside) ^ 2
    
    smoothed_map = hp.sphtfunc.alm2map(alm, nside, lmax=2*nside, pixwin=true, verbose=false, fwhm=input_fwhm*pi/10800.)

    return smoothed_map
