import healpy as hp
import numpy as np

def noise_calc_simle(pol_sen, nside, noise_seed):

    resol_pixel = hp.nside2resol(nside, arcmin = True) 

    sigma = pol_sen / resol_pixel

    pixel = hp.nside2npix(nside)
    
    # average, σ, N
    np.random.seed(noise_seed[0])
    
    noise_map_1 = np.random.normal(0.0, sigma, pixel)

    np.random.seed(noise_seed[1])
    
    noise_map_2 = np.random.normal(0.0, sigma, pixel)

    np.random.seed(noise_seed[2])
    
    noise_map_3 = np.random.normal(0.0, sigma, pixel)

    noise_map = [noise_map_1, noise_map_2, noise_map_3]

    return noise_map, sigma

def noise_map_make(freq_band, nside, noise_seed_array):

    noise_map_array = []
    sigma_array = []

    for i, freq in enumerate(freq_band):

        noise_seed = noise_seed_array[i]

        noise_map, sigma = noise_calc(freq, nside, noise_seed)

        noise_map_array.append(noise_map)

        sigma_array.append(sigma)

    return noise_map_array, sigma_array

    
def noise_calc(freq, nside, noise_seed):

    resol_pixel = hp.nside2resol(nside, arcmin = True) 

    if 40 == freq:
        pol_sen = 37.42 

    elif 50 == freq:
        pol_sen = 33.46

    elif 60 == freq:
        pol_sen = 21.31

    elif 68 == freq:
        pol_sen = 16.87

    elif 78 == freq:
        pol_sen = 12.07

    elif 89 == freq:
        pol_sen = 11.30

    elif 100 == freq:
        pol_sen = 6.56

    elif 119 == freq:
        pol_sen = 4.58

    elif 140 == freq:
        pol_sen = 4.79

    elif 166 == freq:
        pol_sen = 5.57

    elif 195 == freq:
        pol_sen = 5.85

    elif 235 == freq:
        pol_sen = 10.79

    elif 280 == freq:
        pol_sen = 13.80

    elif 337 == freq:
        pol_sen = 21.95

    elif 402 == freq:
        pol_sen = 47.45

    else:
        print("正しい周波数をいれてね")

    sigma = pol_sen / resol_pixel

    pixel = hp.nside2npix(nside)
    
    # average, σ, N
    np.random.seed(noise_seed[0])
    
    noise_map_1 = np.random.normal(0.0, sigma, pixel)

    np.random.seed(noise_seed[1])
    
    noise_map_2 = np.random.normal(0.0, sigma, pixel)

    np.random.seed(noise_seed[2])
    
    noise_map_3 = np.random.normal(0.0, sigma, pixel)

    noise_map = [noise_map_1, noise_map_2, noise_map_3]

    return noise_map, sigma

def noise_map_make(freq_band, nside, noise_seed_array):

    noise_map_array = []
    sigma_array = []

    for i, freq in enumerate(freq_band):

        noise_seed = noise_seed_array[i]

        noise_map, sigma = noise_calc(freq, nside, noise_seed)

        noise_map_array.append(noise_map)

        sigma_array.append(sigma)

    return noise_map_array, sigma_array



        

        

        

    

    
