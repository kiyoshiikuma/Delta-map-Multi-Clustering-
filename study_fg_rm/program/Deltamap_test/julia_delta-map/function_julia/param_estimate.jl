using Optim
using PyCall

include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/function_julia/make_data_m.jl")
include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/function_julia/Delta_map.jl")
include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/function_julia/Delta_map_alpha.jl")
include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/function_julia/cov_mat_calc.jl")
include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/function_julia/calc_noise.jl");


# delta_mapを用いてパラメータ推定を行う

function minimize_function_chi_sq(beta_s, r, cov_mat_scal, cov_mat_tens, x_map, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)
    
    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, 140, nside, cmb_data, "s1", freq_bs, beta_s, freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    return x_with_noise' / cov_mat * x_with_noise

end

function minimize_function_like(r, beta_s, cov_mat_scal, cov_mat_tens, x_map, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, 140, nside, cmb_data, "s1", freq_bs, beta_s, freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    det_C = cholesky_logdet(cov_mat)
    
    return x_with_noise' / cov_mat * x_with_noise + det_C

end


function r_estimate(r_ini, beta_ini, num_iterations, cov_mat_scal, cov_mat_tens, x_map, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)
    
    #num_iterations = 10

    r_array = []
    beta_array = []

    # 関数を用意

    # chi二乗
    minimize_function_chi_sq_new(beta_s, r) = minimize_function_chi_sq(beta_s, r, cov_mat_scal, cov_mat_tens, x_map, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    # likelihood
    minimize_function_like_new(r, beta_s) = minimize_function_like(r, beta_s, cov_mat_scal, cov_mat_tens, x_map, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)



    for i in 1:num_iterations

        if i == 1

            #======================chi二乗　を最小=================================#
            # 固定するr
            fixed_r = r_ini

            # r固定した関数を定義
            objective_function_p(beta_s) = minimize_function_chi_sq_new(beta_s, fixed_r)

            # 最適化オプションを設定（1変数 → Brent）
            opt = optimize(objective_function_p, -4.0, -2.0, Brent())

            # 結果
            optimum_beta_s = Optim.minimizer(opt)
            #println("Optimum beta_s: ", optimum_beta_s)

            #optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(beta_array, optimum_beta_s)

            #=========================likelihood 最小====================================#

            # betaを更新
            fixed_beta_s = optimum_beta_s

            # beta固定した関数を定義
            objective_function_r(r) = minimize_function_like_new(r, fixed_beta_s)

            # 最適化オプションを設定（1変数 → Brent）
            opt = optimize(objective_function_r, 0., 1., Brent())
        
            # 結果
            optimum_r = Optim.minimizer(opt)
            #println("Optimum r: ", optimum_r)

            #optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(r_array, optimum_r)
        
        else

            #======================chi二乗　を最小=================================#
        
            # r更新
            fixed_r = r_array[i-1]

            # r固定した関数を定義
            objective_function_pp(beta_s) = minimize_function_chi_sq_new(beta_s, fixed_r)

            # 最適化オプションを設定（1変数 → Brent）
            opt = optimize(objective_function_pp, -4.0, -2.0, Brent())

            # 結果
            optimum_beta_s = Optim.minimizer(opt)
            #println("Optimum beta_s: ", optimum_beta_s)

        
            #optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(beta_array, optimum_beta_s)

            #=========================likelihood 最小====================================#

            # beta更新
            fixed_beta_s = optimum_beta_s

            # beta固定した関数を定義
            objective_function_rr(r) = minimize_function_like_new(r, fixed_beta_s)

            # 最適化オプションを設定（1変数 → Brent）
            opt = optimize(objective_function_rr, 0., 1., Brent())
        
            # 結果
            optimum_r = Optim.minimizer(opt)
            #println("Optimum r: ", optimum_r)

            #optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(r_array, optimum_r)
    
        end
    
    end

    return r_array, beta_array
    
end