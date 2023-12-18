using LinearAlgebra
using Healpix


function D_element(freq, freq_bs = 23*10^9, beta_s = -3., freq_bd = 353*10^9, beta_d = 1.5, T_d =20.1,)

    # x = (f / T) * (h / k)
    x = (freq / 2.725) / (2.083661912 * 10^10)
    
    g_freq = ((exp(x) - 1)^2) / (exp(x) * x^2) * 1000.0
    
    s = g_freq * (freq/freq_bs)^(beta_s)
    
    ss = g_freq * (freq/freq_bs)^(beta_s) * log(freq/freq_bs)

    x_d = (freq / T_d) / (2.083661912 * 10^10)

    x_bd = (freq_bd / T_d) / (2.083661912 * 10^10)

    d = g_freq * (freq / freq_bd)^(beta_d + 1) * ((exp(x_bd)-1)/(exp(x_d)-1))

    dd = d * log(freq / freq_bd)

    ddd = d * (((x_d * exp(x_d)) / (exp(x_d)-1)) - (x_bd * exp(x_bd)) / (exp(x_bd)-1)) / T_d
    
    return s, ss, d, dd, ddd
    
end


function calc_D_matrix(freq_band, which_model, nside)
    
    N_pix = nside2npix(nside)
    diag = Matrix{Float64}(I, 2 * N_pix, 2 * N_pix)

    D_blocks = []

    for freq in freq_band
        
        s, ss, d, dd, ddd = D_element(freq * 1e9)

        if which_model == "s1"
            
            push!(D_blocks, [diag s * diag ss * diag])
            
        elseif which_model == "d1"
            
            push!(D_blocks, [diag d * diag dd * diag ddd * diag])
            
        elseif which_model == "d1 and s1"
            
            push!(D_blocks, [diag s * diag ss * diag d * diag dd * diag ddd * diag])
            
        else
            println("正しいモデルを入力してね")
            
            return
        end
    end

    D_matrix = vcat(D_blocks...)
    
    return D_matrix
end


function calc_D_matrix_param(freq_band, which_model, nside, freq_bs, beta_s, freq_bd, beta_d, T_d)

    N_pix = nside2npix(nside)

    diag = Matrix{Float64}(I, 2 * N_pix, 2 * N_pix)

    D_blocks = []

    for freq in freq_band
        
        s, ss, d, dd, ddd = D_element(freq * 10^9, freq_bs, beta_s, freq_bd, beta_d, T_d)
        
        if which_model == "s1"
            
            push!(D_blocks, [diag s * diag ss * diag])
            
        elseif which_model == "d1"
            
            push!(D_blocks, [diag d * diag dd * diag ddd * diag])
            
        elseif which_model == "d1 and s1"
            
            push!(D_blocks, [diag s * diag ss * diag d * diag dd * diag ddd * diag])
            
        else
            println("正しいモデルを入力してね")
            
            return
        end
    end

    D_matrix = vcat(D_blocks...)
    
    return D_matrix
    
end

# クリーンマプの計算

function calc_clean_map(data_m, D, nside)
    
    #CMB_map = bslash((D.T @ bslash(N, D)), D.T @ bslash (N, m))
    # N + M +3 = Nfreq ===> just number of frequency bands (Synch + Dust)
    # N + 2 = Nfreq ===> just number of frequency bands (Dust)
    # M + 2 = Nfreq ===> just number of frequency bands (Synch)

    N_pix = nside2npix(nside)
    
    clean_map = D \ data_m

    CMB_map_Q = clean_map[N_pix*0 + 1 : N_pix*1]
    CMB_map_U = clean_map[N_pix*1 + 1 : N_pix*2]

    return CMB_map_Q, CMB_map_U

end


function calc_spectrum(input_Q, input_U, clean_Q, clean_U, nside)
  
    N_pix = Healpix.nside2npix(nside)

    II = zeros(Float64, N_pix)
    
    Clean_map = [II, clean_Q, clean_U]
    
    Input_map = [II, input_Q, input_U]

    # ell is 1 - 3 nside -1
    ell = 1:1:3*nside - 1
    
    Clean_map_cl = hp.sphtfunc.anafast(Clean_map)
    
    Input_map_cl = hp.sphtfunc.anafast(Input_map)
    
    Clean_map_dl = Clean_map_cl * ell * (ell + 1) / (2 * math.pi)
    
    Input_map_dl = Input_map_cl * ell * (ell + 1) / (2 * math.pi)

    return Clean_map_cl, Clean_map_dl, Input_map_cl, Input_map_dl, ell

end





