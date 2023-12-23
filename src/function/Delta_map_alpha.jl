using LinearAlgebra
using Healpix


include("masked_map.jl")


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
    
    return s, ss, d, dd, ddd, g_freq
    
end


function d_s_vec_calc_pra(freq, freq_bs, beta_s, freq_bd, beta_d, T_d)

    s, ss, d, dd, ddd, g_freq = D_element(freq * 10^9, freq_bs, beta_s, freq_bd, beta_d, T_d)

    s_vec = [s, ss]

    d_vec = [d, dd, ddd]

    return s_vec, d_vec, g_freq

end


function Amat_calc_pra(freq_band, cmb_freq, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)
    
    A_s = []
    A_d = []
    
    for freq in freq_band
        s_vec, d_vec, g_freq = d_s_vec_calc_pra(freq, freq_bs, beta_s, freq_bd, beta_d, T_d)

        if cmb_freq != freq
            
            if which_model == "s1"
                
                push!(A_s, s_vec)
                
            elseif which_model == "d1"
                
                push!(A_d, d_vec)
                
            elseif which_model == "d1 and s1"
                
                push!(A_s, s_vec)
                push!(A_d, d_vec)
                
            end
            
        else
            # pass
        end
    end

    if which_model == "s1"
        A = hcat(A_s...)
        
    elseif which_model == "d1"
        A = hcat(A_d...)
        
    elseif which_model == "d1 and s1"
        A_d = hcat(A_d...)
        A_s = hcat(A_s...)
        A = vcat(A_d, A_s)
        
    else
        println("正しいモデルを入力してね")
    end
    
    return A
end


function calc_x_pra(freq_band, cmb_freq, nside, cmb_map, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)
    
    Q_map, U_map = make_input_map(cmb_map, freq_band, which_model, nside)
    
    Amat = Amat_calc_pra(freq_band, cmb_freq, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)
    
    s_cmb, d_cmb, g_cmb = d_s_vec_calc_pra(cmb_freq, freq_bs, beta_s, freq_bd, beta_d, T_d)

    if which_model == "s1"
        vec_cmb = s_cmb
        
    elseif which_model == "d1"
        vec_cmb = d_cmb
        
    elseif which_model == "d1 and s1"
        vec_cmb = [d_cmb; s_cmb]
        
    end

    alpha_i = - Amat \ vec_cmb

    Q = zeros(size(Q_map[1]))
    U = zeros(size(U_map[1]))

    sum_alpha = sum(alpha_i)

    freq_list = collect(freq_band)
    cmb_index = findall(x -> x == cmb_freq, freq_list)[1]

    insert!(alpha_i, cmb_index, 1.)
    
    for (i, freq) in enumerate(freq_band)
        Q .+= alpha_i[i] * Q_map[i]
        U .+= alpha_i[i] * U_map[i]
    end

    Q /= 1 + sum_alpha
    U /= 1 + sum_alpha

    x = vcat(Q, U)

    return Q, U, x
end



function make_input_map(cmb, freq_band, which_model, nside)
    
    band = length(freq_band)
    
    N_pix = Healpix.nside2npix(nside)
    
    data_m = zeros(Float64, 2 * N_pix * band)

    Q_map = []
    U_map = []
    
    for (ii, freq) in enumerate(freq_band)
        
        # file name
        dir = "../../make_map/fg_map_make/fg_map/"
        nside_name = "nside_"
        nside_n = string(nside)

        fg_synch = "Synch_"
        GHz = "_GHz_"
        
        # Synchrotron name
        Synch_name = string(dir, fg_synch, freq, GHz, nside_name, nside_n)
        
        fg_dust = "Dust_"
        
        # Dust name
        Dust_name = string(dir, fg_dust, freq, GHz, nside_name, nside_n)
        
        #synch_data = np.load(Synch_name)
        synch_data_I = Healpix.readMapFromFITS(Synch_name, 1, Float64)
        synch_data_Q = Healpix.readMapFromFITS(Synch_name, 2, Float64)
        synch_data_U = Healpix.readMapFromFITS(Synch_name, 3, Float64)
        
        
        #dust_data = np.load(Dust_name)
        Dust_data_I = Healpix.readMapFromFITS(Dust_name, 1, Float64)
        Dust_data_Q = Healpix.readMapFromFITS(Dust_name, 2, Float64)
        Dust_data_U = Healpix.readMapFromFITS(Dust_name, 3, Float64)
        
        if which_model == "s1"
            
            m_Q = synch_data_Q + cmb[2, :]
            m_U = synch_data_U + cmb[3, :]

            push!(Q_map, m_Q)
            push!(U_map, m_U)
            
        elseif which_model == "d1"

            # data m[I, Q, U]
            m_Q = Dust_data_Q + cmb[2, :]
            m_U = Dust_data_U + cmb[3, :]

            push!(Q_map, m_Q)
            push!(U_map, m_U)

        elseif which_model == "d1 and s1"

            # data m[I, Q, U]
            m_Q = synch_data_Q + Dust_data_Q + cmb[2, :]
            m_U = synch_data_U + Dust_data_U + cmb[3, :]

            push!(Q_map, m_Q)
            push!(U_map, m_U)

        else
            print("正しいモデル入れてね")

        end
        
    end

    return Q_map, U_map

end


# ======================maskありのinput_map=========================================#
# =================================================================================#

function masked_make_input_map(cmb, freq_band, which_model, nside)
    
    band = length(freq_band)
    
    N_pix = Healpix.nside2npix(nside)
    
    data_m = zeros(Float64, 2 * N_pix * band)

    Q_map = []
    U_map = []
    
    for (ii, freq) in enumerate(freq_band)
        
        # file name
        dir = "../../make_map/fg_map_make/fg_map/"
        nside_name = "nside_"
        nside_n = string(nside)

        fg_synch = "Synch_"
        GHz = "_GHz_"
        
        # Synchrotron name
        Synch_name = string(dir, fg_synch, freq, GHz, nside_name, nside_n)
        
        fg_dust = "Dust_"
        
        # Dust name
        Dust_name = string(dir, fg_dust, freq, GHz, nside_name, nside_n)
        
        #synch_data = np.load(Synch_name)
        synch_data_I = Healpix.readMapFromFITS(Synch_name, 1, Float64)
        synch_data_Q = Healpix.readMapFromFITS(Synch_name, 2, Float64)
        synch_data_U = Healpix.readMapFromFITS(Synch_name, 3, Float64)
        
        
        #dust_data = np.load(Dust_name)
        Dust_data_I = Healpix.readMapFromFITS(Dust_name, 1, Float64)
        Dust_data_Q = Healpix.readMapFromFITS(Dust_name, 2, Float64)
        Dust_data_U = Healpix.readMapFromFITS(Dust_name, 3, Float64)

        if which_model == "s1"

            # ここでマスクをすればいい
            m_Q = masked_data_calc(synch_data_Q + cmb[2, :])
            m_U = masked_data_calc(synch_data_U + cmb[3, :])

            push!(Q_map, m_Q)
            push!(U_map, m_U)
            
        elseif which_model == "d1"

            # data m[I, Q, U]
            m_Q = masked_data_calc(Dust_data_Q + cmb[2, :])
            m_U = masked_data_calc(Dust_data_U + cmb[3, :])

            push!(Q_map, m_Q)
            push!(U_map, m_U)

        elseif which_model == "d1 and s1"

            # data m[I, Q, U]
            m_Q = masked_data_calc(synch_data_Q + Dust_data_Q + cmb[2, :])
            m_U = masked_data_calc(synch_data_U + Dust_data_U + cmb[3, :])

            push!(Q_map, m_Q)
            push!(U_map, m_U)

        else
            print("正しいモデル入れてね")

        end
        
    end

    return Q_map, U_map

end


function masked_calc_x_pra(freq_band, cmb_freq, nside, cmb_map, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)
    
    Q_map, U_map = masked_make_input_map(cmb_map, freq_band, which_model, nside)
    
    Amat = Amat_calc_pra(freq_band, cmb_freq, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)
    
    s_cmb, d_cmb, g_cmb = d_s_vec_calc_pra(cmb_freq, freq_bs, beta_s, freq_bd, beta_d, T_d)

    if which_model == "s1"
        vec_cmb = s_cmb
        
    elseif which_model == "d1"
        vec_cmb = d_cmb
        
    elseif which_model == "d1 and s1"
        vec_cmb = [d_cmb; s_cmb]
        
    end

    alpha_i = - Amat \ vec_cmb

    Q = zeros(size(Q_map[1]))
    U = zeros(size(U_map[1]))

    sum_alpha = sum(alpha_i)

    freq_list = collect(freq_band)
    cmb_index = findall(x -> x == cmb_freq, freq_list)[1]

    insert!(alpha_i, cmb_index, 1.)
    
    for (i, freq) in enumerate(freq_band)
        Q .+= alpha_i[i] * Q_map[i]
        U .+= alpha_i[i] * U_map[i]
    end

    Q /= 1 + sum_alpha
    U /= 1 + sum_alpha

    x = vcat(Q, U)

    return Q, U, x
end


#===========logdetの計算==========#

function log_det(A)
    try
        # 固有値に分解、対数とって足し合わせ
        logdet = sum(log.(eigvals(Symmetric(A))))
        return logdet
    catch
        throw(ArgumentError("Matrix is numerically singular!"))
    end
end


function cholesky_logdet(A)
    # ¥コレスキー分解
    L = cholesky(A).L
    
    # 行列の対数行列式の2倍
    logdet = 2 * sum(log.(diag(L)))

    return logdet
end
