using Optim
using PyCall
using Healpix
using NPZ

@pyimport sys
pushfirst!(PyVector(pyimport("sys")["path"]), "../../src/function")

@pyimport smoothing_map as sm
@pyimport cmb_make_file as cmb_make
@pyimport scipy.optimize as scipy_optimize

include("make_data_m.jl")
include("Delta_map.jl")
include("Delta_map_alpha.jl") 
include("cov_mat_calc.jl")
include("calc_noise.jl");

# delta_mapを用いてパラメータ推定を行う
# 最小化する関数を下のように作成

function minimize_function_chi_sq(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)
    
    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    return x_with_noise' / cov_mat * x_with_noise

end

function minimize_function_like(r, beta_s, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    det_C = cholesky_logdet(cov_mat)
    
    return x_with_noise' / cov_mat * x_with_noise + det_C

end

function minimize_function_chi_sq_s1(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)
    
    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, beta_s[1], freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    return x_with_noise' / cov_mat * x_with_noise

end


function minimize_function_chi_sq_d1_and_d1_s1(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)
    
    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    return x_with_noise' / cov_mat * x_with_noise

end

function minimize_function_like_free(r, beta_s, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    cov_mat = calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r[1]) + art_noise_cov_mat

    Q_1, U_1, x_map = calc_x_pra(freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, beta_s, freq_bd, beta_d, T_d)

    x_with_noise = x_map + noise_map
    
    det_C = cholesky_logdet(cov_mat)
    
    return x_with_noise' / cov_mat * x_with_noise + det_C

end


#　各ファイルを読み込むから、r_like_estimateより早く動く

function r_like_estimate_fast(r_ini, num_iterations, accuracy_r, accuracy_like, random_seed_cmb, seed_syn, nside, r_input, freq_band, which_model, cmb_freq, cov_mat_scal, cov_mat_tens)

    #================CMBを読む,  r と nside, seedを指定=================================#

    # file name
    dir =  "../../make_map/cmb_map_make/cmb_map/"

    # r_0.01_nside_4_seed_1みたいな感じ"

    r_name = "r_"
    r_n = string(r_input)
    
    nside_name = "_nside_"
    nside_n = string(nside)

    # cmbファイルは seed が 0 ~ 999なのでそれを指定する(1~1000に直せばいいんだけどさ)
    seed_name = "_seed_"
    seed_n = string(seed_syn)
    
    # Synchrotron name
    cmb_name = string(dir, r_name, r_n, nside_name, nside_n, seed_name, seed_n)
    
    #synch_data = np.load(Synch_name)
    cmb_I = Healpix.readMapFromFITS(cmb_name, 1, Float64)
    cmb_Q = Healpix.readMapFromFITS(cmb_name, 2, Float64)
    cmb_U = Healpix.readMapFromFITS(cmb_name, 3, Float64)

    cmb_data = [cmb_I'; cmb_Q'; cmb_U']
    
    #============================cmbパワースペクトルを読む=====================================#
    #cl_scal, cl_tens = cmb_make.cmb_cl_calc(nside, random_seed_cmb, seed_syn);

    #パワースペクトルの計算

    #cl_scal_EE = cl_scal[:,2][1:2*nside+1]
    #cl_scal_BB = cl_scal[:,3][1:2*nside+1]
    #cl_tens_EE = cl_tens[:,2][1:2*nside+1]
    #cl_tens_BB = cl_tens[:,3][1:2*nside+1];

    # 共分散行列の計算
    #cov_mat_scal = calc_cmb_cov_mat(cl_scal_EE, cl_scal_BB, nside)
    #cov_mat_tens = calc_cmb_cov_mat(cl_tens_EE, cl_tens_BB, nside);

    # 変数に指定しておくほうがいいかも 何回も読み込むのはめんどくさい,,,,

    #cov_mat_scal = npzread("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/covariance_matrix/cov_mat_scal.npy")
    #cov_mat_tens = npzread("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/covariance_matrix/cov_mat_tens.npy");
    
    #============================パラメータ指定(重要じゃない、てか実質使用されない)====================================#
    
    # 各パラメータの指定
    beta_s, freq_bs, freq_bd, beta_d, T_d = -3, 23*10^9, 353*10^9, 1.5, 20.1
    
    #============================ノイズ作成=====================================#
    
    #　人工ノイズの作成
    noise_seed = [1, 2, 3]
    pol_sen = 0.2 # μK

    art_noise_cov_mat =  calc_noise_cov_mat(pol_sen, nside, noise_seed);
    art_noise_map, sigma = calc_noise_map(pol_sen, nside, noise_seed);

    noise_map = [art_noise_map[1]; art_noise_map[2]];

    #============================パラメータ推定=====================================#
    
    r_array = []
    beta_array = []
    like_array = []

    # 関数を用意

    # chi二乗
    minimize_function_chi_sq_new(beta_s, r) = minimize_function_chi_sq(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    # likelihood
    minimize_function_like_new(r, beta_s) = minimize_function_like(r, beta_s, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

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

            optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(r_array, optimum_r)
            push!(like_array, optimum_value)
        
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

            optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(r_array, optimum_r)
            push!(like_array, optimum_value)

            # 条件
            if abs(-r_array[i] + r_array[i-1]) <= accuracy_r && abs(like_array[i] - like_array[i-1]) <= accuracy_like
            
                break
            end
    
        end
    
    end

    return r_array, beta_array, r_array[end], beta_array[end]
    
end


#=============================scipy minimize========================================#

function minimize_r_like_estimate_fast(r_ini, num_iterations, accuracy_r, accuracy_like, random_seed_cmb, seed_syn, nside, r_input, freq_band, which_model, cmb_freq, cov_mat_scal, cov_mat_tens)

    #================CMBを読む,  r と nside, seedを指定=================================#

    # file name
    dir = "../../make_map/cmb_map_make/cmb_map/"

    # r_0.01_nside_4_seed_1みたいな感じ"

    r_name = "r_"
    r_n = string(r_input)
    
    nside_name = "_nside_"
    nside_n = string(nside)

    # cmbファイルは seed が 0 ~ 999なのでそれを指定する(1~1000に直せばいいんだけどさ)
    seed_name = "_seed_"
    seed_n = string(seed_syn)
    
    # Synchrotron name
    cmb_name = string(dir, r_name, r_n, nside_name, nside_n, seed_name, seed_n)
    
    #synch_data = np.load(Synch_name)
    cmb_I = Healpix.readMapFromFITS(cmb_name, 1, Float64)
    cmb_Q = Healpix.readMapFromFITS(cmb_name, 2, Float64)
    cmb_U = Healpix.readMapFromFITS(cmb_name, 3, Float64)

    cmb_data = [cmb_I'; cmb_Q'; cmb_U']
    
    #============================cmbパワースペクトルを読む=====================================#
    #cl_scal, cl_tens = cmb_make.cmb_cl_calc(nside, random_seed_cmb, seed_syn);

    #パワースペクトルの計算

    #cl_scal_EE = cl_scal[:,2][1:2*nside+1]
    #cl_scal_BB = cl_scal[:,3][1:2*nside+1]
    #cl_tens_EE = cl_tens[:,2][1:2*nside+1]
    #cl_tens_BB = cl_tens[:,3][1:2*nside+1];

    # 共分散行列の計算
    #cov_mat_scal = calc_cmb_cov_mat(cl_scal_EE, cl_scal_BB, nside)
    #cov_mat_tens = calc_cmb_cov_mat(cl_tens_EE, cl_tens_BB, nside);

    # 変数に指定しておくほうがいいかも 何回も読み込むのはめんどくさい,,,,

    #cov_mat_scal = npzread("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/covariance_matrix/cov_mat_scal.npy")
    #cov_mat_tens = npzread("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/covariance_matrix/cov_mat_tens.npy");
    
    #============================パラメータ指定(重要じゃない、てか実質使用されない)====================================#
    
    # 各パラメータの指定
    beta_s, freq_bs, freq_bd, beta_d, T_d = -3, 23*10^9, 353*10^9, 1.5, 20.1
    
    #============================ノイズ作成=====================================#
    
    #　人工ノイズの作成
    noise_seed = [1, 2, 3]
    pol_sen = 0.2 # μK

    art_noise_cov_mat =  calc_noise_cov_mat(pol_sen, nside, noise_seed);
    art_noise_map, sigma = calc_noise_map(pol_sen, nside, noise_seed);

    noise_map = [art_noise_map[1]; art_noise_map[2]];

    #============================パラメータ推定=====================================#
    
    r_array = []
    beta_s_array = []
    beta_d_array = []
    T_d_array = []
    like_array = []

    # 関数を用意する

    # chi二乗
    minimize_function_chi_sq_s1_new(beta_s, r, beta_d, T_d) = minimize_function_chi_sq_s1(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    minimize_function_chi_sq_d1_and_d1_s1_new(beta_s, r, beta_d, T_d) = minimize_function_chi_sq_d1_and_d1_s1(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    # likelihood
    minimize_function_like_free_new(r, beta_s, beta_d, T_d) = minimize_function_like_free(r, beta_s, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    for i in 1:num_iterations

        if i == 1

            #======================chi二乗　を最小=================================#
            # 固定する r
            fixed_r = r_ini

            if which_model == "s1"

                fixed_beta_d, fixed_T_d = 1.5, 20.1

                # r固定した関数を定義
                objective_function_s1(beta_s, r) = minimize_function_chi_sq_s1_new(beta_s, r, fixed_beta_d, fixed_T_d)

                objective_function_s1_new(beta_s) = objective_function_s1(beta_s, fixed_r)

                # beta の初期値と範囲
                x0 = [-3.]

                bounds = [(-4., -2.)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_s1_new, x0, bounds=bounds)

                optimum_beta_s = result["x"][1]
                push!(beta_s_array, optimum_beta_s)

                #=========================likelihood 最小====================================#

                # betaを更新
                fixed_beta_s = optimum_beta_s

                # beta固定した関数を定義
                objective_function_r_s1(r) = minimize_function_like_free_new(r, fixed_beta_s, fixed_beta_d, fixed_T_d)

                # r の初期値と範囲
                x0 = [0.]

                bounds = [(0, 1)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_r_s1, x0, bounds=bounds)

                optimum_r = result["x"][1]
                optimum_value = result["fun"]

                push!(r_array, optimum_r)
                push!(like_array, optimum_value)

            #======================chi二乗　を最小=================================#
            elseif which_model == "d1"

                fixed_beta_s = -3
                
                # r固定した関数を定義
                objective_function_d1(r, beta_d, T_d) = minimize_function_chi_sq_d1_and_d1_s1_new(fixed_beta_s, r, beta_d, T_d)

                objective_function_d1_new(vars) = objective_function_d1(fixed_r, vars[1], vars[2])
                
                # 初期値 beta_d, T_d
                initial_values = [1.5, 20]

                bounds = [(1., 2.), (5, 40)]
                
                # scipy.optimizeのminimize関数を使用
                result = scipy_optimize.minimize(objective_function_d1_new, initial_values, bounds = bounds)

                # 結果
                optimum_beta_d = result["x"][1]
                optimum_T_d = result["x"][2]

                push!(beta_d_array, optimum_beta_d)
                push!(T_d_array, optimum_T_d)

                #=========================likelihood 最小====================================#

                # beta_d, T_d を更新
                fixed_beta_d = optimum_beta_d
                fixed_T_d = optimum_T_d
                
                # beta固定した関数を定義
                objective_function_r_d1(r) = minimize_function_like_free_new(r, fixed_beta_s, fixed_beta_d, fixed_T_d)

                x0 = [0.]

                bounds = [(0, 1)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_r_d1, x0, bounds=bounds)

                optimum_r = result["x"][1]
                optimum_value = result["fun"]

                push!(r_array, optimum_r)
                push!(like_array, optimum_value)
                
            #======================chi二乗　を最小=================================#
            elseif which_model == "d1 and s1"

                # 初期値
                initial_values = [-3, 1.5, 20]

                bounds = [(-4., -2.), (1., 2.), (5, 40)]

                # r固定した関数を定義

                objective_function_d1_s1(vars) = minimize_function_chi_sq_d1_and_d1_s1_new(vars[1], fixed_r, vars[2], vars[3])

                # scipy.optimizeのminimize関数を使用
                result = scipy_optimize.minimize(objective_function_d1_s1, initial_values, bounds = bounds)

                # 結果
                optimum_beta_s = result["x"][1]
                optimum_beta_d = result["x"][2]
                optimum_T_d = result["x"][3]

                push!(beta_s_array, optimum_beta_s)
                push!(beta_d_array, optimum_beta_d)
                push!(T_d_array, optimum_T_d)

                #=========================likelihood 最小====================================#

                # beta_s, beta_d, T_d を更新
                fixed_beta_s = optimum_beta_s
                fixed_beta_d = optimum_beta_d
                fixed_T_d = optimum_T_d
                
                # beta固定した関数を定義
                objective_function_r_d1_s1(r) = minimize_function_like_free_new(r, fixed_beta_s, fixed_beta_d, fixed_T_d)

                x0 = [0.]

                bounds = [(0, 1)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_r_d1_s1, x0, bounds=bounds)

                optimum_r = result["x"][1]
                optimum_value = result["fun"]

                push!(r_array, optimum_r)
                push!(like_array, optimum_value)
            
            else
                Print("正しいモデルを入れてね")

            end

        else

            #======================chi二乗　を最小=================================#
            # r更新
            fixed_r = r_array[i-1]
            
            if which_model == "s1"
                
                fixed_beta_d, fixed_T_d = 1.5, 20.1
                
                # r固定した関数を定義
                objective_function_s1_newnew(beta_s, r) = minimize_function_chi_sq_s1_new(beta_s, r, fixed_beta_d, fixed_T_d)

                objective_function_s1_newnewnew(beta_s) = objective_function_s1_newnew(beta_s, fixed_r)

                fixed_beta_s = r_array[i-1]

                # r の初期値と範囲
                x0 = [-3.]

                bounds = [(-4., -2.)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_s1_newnewnew, x0, bounds=bounds)

                optimum_beta_s = result["x"][1]
                push!(beta_s_array, optimum_beta_s)

                #=========================likelihood 最小====================================#

                # betaを更新
                fixed_beta_s = optimum_beta_s

                # beta固定した関数を定義
                objective_function_r_s1_new(r) = minimize_function_like_free_new(r, fixed_beta_s, fixed_beta_d, fixed_T_d)

                x0 = [0.]

                bounds = [(0, 1)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_r_s1_new, x0, bounds=bounds)

                optimum_r = result["x"][1]
                optimum_value = result["fun"]

                push!(r_array, optimum_r)
                push!(like_array, optimum_value)
            
            #======================chi二乗　を最小=================================#
            
            elseif which_model == "d1"

                fixed_beta_s = -3

                bounds = [(1., 2.), (5, 40)]
                
                # r固定した関数を定義
                objective_function_d1_newnew(r, beta_d, T_d) = minimize_function_chi_sq_d1_and_d1_s1_new(fixed_beta_s, r, beta_d, T_d)

                objective_function_d1_newnewnew(vars) = objective_function_d1_newnew(fixed_r, vars[1], vars[2])
                
                # 初期値
                initial_values = [1.5, 20]

                # scipy.optimizeのminimize関数を使用
                result = scipy_optimize.minimize(objective_function_d1_newnewnew, initial_values, bounds = bounds)

                # 結果
                optimum_beta_d = result["x"][1]
                optimum_T_d = result["x"][2]

                push!(beta_d_array, optimum_beta_d)
                push!(T_d_array, optimum_T_d)

                #=========================likelihood 最小====================================#

                # beta_d, T_d を更新
                fixed_beta_d = optimum_beta_d
                fixed_T_d = optimum_T_d
                
                # beta固定した関数を定義
                objective_function_r_d1_new(r) = minimize_function_like_free_new(r, fixed_beta_s, fixed_beta_d, fixed_T_d)

                x0 = [0.]

                bounds = [(0, 1)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_r_d1_new, x0, bounds=bounds)

                optimum_r = result["x"][1]
                optimum_value = result["fun"]

                push!(r_array, optimum_r)
                push!(like_array, optimum_value)
                

            #======================chi二乗　を最小=================================#
            elseif which_model == "d1 and s1"

                # 初期値
                initial_values = [-3, 1.5, 20]

                bounds = [(-4., -2.), (1., 2.), (5, 40)]

                # r固定した関数を定義

                objective_function_d1_s1_new(vars) = minimize_function_chi_sq_d1_and_d1_s1_new(vars[1], fixed_r, vars[2], vars[3])

                # scipy.optimizeのminimize関数を使用
                result = scipy_optimize.minimize(objective_function_d1_s1_new, initial_values, bounds = bounds)

                # 結果
                optimum_beta_s = result["x"][1]
                optimum_beta_d = result["x"][2]
                optimum_T_d = result["x"][3]

                push!(beta_s_array, optimum_beta_s)
                push!(beta_d_array, optimum_beta_d)
                push!(T_d_array, optimum_T_d)

                #=========================likelihood 最小====================================#

                # beta_s, beta_d, T_d を更新
                fixed_beta_s = optimum_beta_s
                fixed_beta_d = optimum_beta_d
                fixed_T_d = optimum_T_d
                
                # beta固定した関数を定義
                objective_function_r_d1_s1_new(r) = minimize_function_like_free_new(r, fixed_beta_s, fixed_beta_d, fixed_T_d)

                x0 = [0.]

                bounds = [(0, 1)]

                # 最適化の実行
                result = scipy_optimize.minimize(objective_function_r_d1_s1_new, x0, bounds=bounds)

                optimum_r = result["x"][1]
                optimum_value = result["fun"]

                push!(r_array, optimum_r)
                push!(like_array, optimum_value)
            
            else
                Print("正しいモデルを入れてね")

            end
            
            # 条件
            if abs(-r_array[i] + r_array[i-1]) <= accuracy_r && abs(like_array[i] - like_array[i-1]) <= accuracy_like
            
                break
            end
    
        end
    
    end

    if which_model == "s1"
        beta_s_estimate = beta_s_array[end]
        beta_d_estimate = beta_d
        T_d_estimate = T_d

    elseif which_model == "d1"
        beta_s_estimate = beta_s
        beta_d_estimate = beta_d_array[end]
        T_d_estimate = T_d_array[end]

    elseif which_model == "d1 and s1"
        beta_s_estimate = beta_s_array[end]
        beta_d_estimate = beta_d_array[end]
        T_d_estimate = T_d_array[end]

    else
        print("正しいモデルを入れて")

    end
    
    return r_array, beta_s_array, beta_d_array, T_d_array, r_array[end], beta_s_estimate, beta_d_estimate, T_d_estimate
    
end



function r_like_estimate(r_ini, beta_ini, num_iterations, accuracy_r, accuracy_like, random_seed_cmb, seed_syn, nside, r_input, freq_band, which_model, cmb_freq)

    # CMBをここで作成 r_inputでrを指定
    
    #seed_syn =5123

    cmb_data = cmb_make.cmb_make_file(nside, r_input, random_seed_cmb, seed_syn);
    cl_scal, cl_tens = cmb_make.cmb_cl_calc(nside, random_seed_cmb, seed_syn);

    #パワースペクトルの計算

    cl_scal_EE = cl_scal[:,2][1:2*nside+1]
    cl_scal_BB = cl_scal[:,3][1:2*nside+1]
    cl_tens_EE = cl_tens[:,2][1:2*nside+1]
    cl_tens_BB = cl_tens[:,3][1:2*nside+1];

    # 共分散行列の計算
    cov_mat_scal = calc_cmb_cov_mat(cl_scal_EE, cl_scal_BB, nside)
    cov_mat_tens = calc_cmb_cov_mat(cl_tens_EE, cl_tens_BB, nside);

    # 各パラメータの指定
    beta_s, freq_bs, freq_bd, beta_d, T_d = -3, 23*10^9, 353*10^9, 1.5, 20.1

    #　人工ノイズの作成
    noise_seed = [1, 2, 3]
    pol_sen = 0.2 # μK

    art_noise_cov_mat =  calc_noise_cov_mat(pol_sen, nside, noise_seed);
    art_noise_map, sigma = calc_noise_map(pol_sen, nside, noise_seed);

    noise_map = [art_noise_map[1]; art_noise_map[2]];

    # ここから、パラメータ推定
    
    r_array = []
    beta_array = []
    like_array = []

    # 関数を用意

    # chi二乗
    minimize_function_chi_sq_new(beta_s, r) = minimize_function_chi_sq(beta_s, r, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

    # likelihood
    minimize_function_like_new(r, beta_s) = minimize_function_like(r, beta_s, cov_mat_scal, cov_mat_tens, noise_map, art_noise_cov_mat, freq_band, cmb_freq, nside, cmb_data, which_model, freq_bs, freq_bd, beta_d, T_d)

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

            optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(r_array, optimum_r)
            push!(like_array, optimum_value)
        
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

            optimum_value = Optim.minimum(opt)
            #println("Minimum value: ", optimum_value)

            push!(r_array, optimum_r)
            push!(like_array, optimum_value)

            # 条件
            accuracy_r
            if abs(-r_array[i] + r_array[i-1]) <= accuracy_r && abs(like_array[i] - like_array[i-1]) <= accuracy_like
            
                break
            end

        end
    
    end

    return r_array, beta_array, r_array[end], beta_array[end]
    
end