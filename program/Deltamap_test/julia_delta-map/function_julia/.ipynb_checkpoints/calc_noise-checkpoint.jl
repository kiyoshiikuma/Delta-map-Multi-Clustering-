using Healpix
using StatsPlots, Distributions, Random
using LinearAlgebra

# ノイズのマップの計算

function calc_noise_map(pol_sen, nside, noise_seed)

    # radian　で取得　→ arcmin
    resol_pixel = nside2resol(nside) * (60) * (180 / pi)

    sigma = pol_sen / resol_pixel

    pixel = nside2npix(nside)

    # average, σ, N
    # Normal(μ,σ)
    d = Normal(0,sigma)

    # case 1
    Random.seed!(noise_seed[1])

    noise_map_1 = rand(d, pixel) 

    # case 2
    Random.seed!(noise_seed[2])

    noise_map_2 = rand(d, pixel) 
    
    # case 3
    Random.seed!(noise_seed[3])

    noise_map_3 = rand(d, pixel) 

    # noise_map I, Q. U
    noise_map = [noise_map_1, noise_map_2, noise_map_3]

    return noise_map, sigma

end


function calc_noise_cov_mat(pol_sen, nside, noise_seed)

    # radian　で取得　→ arcmin
    resol_pixel = nside2resol(nside) * (60) * (180 / pi)

    sigma = pol_sen / resol_pixel

    pixel = nside2npix(nside)
    
    noise_cov_mat = Matrix{Float64}(I, 2*pixel, 2*pixel) .* sigma^2

    return noise_cov_mat

end

    