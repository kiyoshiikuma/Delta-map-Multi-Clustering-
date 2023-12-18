using WignerD
using Healpix
using Base.Threads
using Printf

#　スピン2の球面調和関数を求める。
# 構造は save_Ylm [npix, ell, m]
# save_Ylm[pix,l+1, lmax+1+m]
# スピン2より l = 0, 1の要素は0なので l+1だから、 l は 2からなので、配列の index 3 からつまっていく、、、、
# lmax+1+m　も mの条件を守りながら詰められていく

# (ex)
# nside = 4
# npix = 192 by l = 2*nside + 1 = 9 by -l < m < l

function calc_Yslm(nside, s)
    
    lmax= 2*nside
    npix = nside2npix(nside)
    res = Resolution(nside)
    
    save_Ylm = zeros(ComplexF64, npix, lmax+1, 2*lmax+1);
    
    @views for pix in 1:npix
        ang = pix2angRing(res, pix)
        for l in 2:lmax
            @inbounds @threads for m in -l:l
                #WignerD.wignerDjmn(l,m,n,φ,θ,ψ)
                save_Ylm[pix,l+1, lmax+1+m] = (-1)^m*sqrt((2*l+1)/(4pi))*WignerD.wignerDjmn(l,-m,s,ang[2],ang[1],0)
            end    
        end
    end

    return save_Ylm
end

# WXの定義より

function calc_WX(nside)

    Ylm = calc_Yslm(nside, 2)
    Ylm_ = calc_Yslm(nside, -2)
    
    Wlm = -1*(Ylm + Ylm_)/2
    Xlm = -1im*(Ylm - Ylm_)/2

    return Wlm, Xlm
    
end

# mでの総和をとる

function calc_sum_m(Flm, Glm, npix1, npix2, l)
    
    # Σ_m W(npix_i) * W^*(npix_j)
    WW_XX = sum(Flm[npix1,l,:] .* conj(Glm[npix2,l,:]))

    return WW_XX

end

# lでの総和をとる ellの範囲に注意、足し合わせは ell =2 から行われる
# sスムージング関数の ell に注意

function calc_sum_l(Flm, Glm, npix1, npix2, cl, wl, beam_pol, nside)
    
    mat_comp = []
    
    lmax = 2*nside
    
    for l in 3 : lmax + 1

        # dl → cl
        #fact = l * (l+1) / (2 * pi)
        fact = 1

        # WW_XX は -lからlだけど 3から、lmax + 1 まで詰まっている
        WW_XX = calc_sum_m(Flm, Glm, npix1, npix2, l)

        # Σ_l Cl * wl^2 * W(npix_i, npix_j)
        Cov_mat = wl[l]^2 * beam_pol[l]^2 * cl[l] * WW_XX / fact

        #=============================================================
        # 各成分を13桁まで保持
        real_part = parse(Float64, @sprintf("%.13f", real(Cov_mat)))
        imag_part = parse(Float64, @sprintf("%.13f", imag(Cov_mat)))
        
        # 13桁の方
        push!(mat_comp, real_part)
        ==============================================================#
    
        # im 部分は小さいので切り落とす
        push!(mat_comp, real(Cov_mat))
        #push!(mat_comp, Cov_mat)

    end 
        
    return sum(mat_comp)

end

# 共分散行列の各ブロック行列を計算する C_QQ, C_QU, C_UQ, C_UU

function calc_cov_mat(F1lm, F2lm, G1lm, G2lm, clEE, clBB, nside)
    
    npix = nside2npix(nside)
    lmax = 2*nside

    # beam, pixwinも考えない場合(= 1) fwhmも適当---------------------------------------#
    fwhm = 0
    
    wl = pixwin(nside, pol=true)

    # 0 < ell < lmax より　julia　の　indexとしては　　1 ~ lmax + 1　　
    wl_P = wl[2][1:lmax + 1] * 0 .+ 1
    
    beam = gaussbeam(fwhm, Int(lmax + 1), pol=true) * 0 .+ 1
    
    #beam_EE = beam[:,2] * 0 .+ 1
    #beam_BB = beam[:,3] * 0 .+ 1

    beam_EE = beam * 0 .+ 1
    beam_BB = beam * 0 .+ 1
    
    #----------------------------------------------------------------------------#

    matrix = zeros((npix,npix))*1.0im
    for i in 1 : npix
        
        @inbounds @threads for j in 1 : npix
            
            comp = calc_sum_l(F1lm, F2lm, i, j, clEE, wl_P, beam_EE, nside) + calc_sum_l(G1lm, G2lm, i, j, clBB, wl_P, beam_BB, nside)

            # i = 1の時　(n1, n1),　(n1, n2)　, (n1, n3)...,(n1, n_npix)　 
            # i = 2の時　(n2, n1),　(n2, n2)　, (n2, n3)...,(n2, n_npix)　 
            # :
            # :
            # i = n_npixの時　(n_npix, n1),　(n_npix, n2)　, (n_npix, n3)...,(n_npix, n_npix)　 
            
            matrix[i,j] = comp
        end
    end
    
    return matrix

end

# signal 共分散行列を計算する C

function calc_cov_mat(Wlm, Xlm, clEE, clBB, nside)

    npix = nside2npix(nside)
    C_QQ = calc_cov_mat(Wlm, Wlm, Xlm, Xlm, clEE, clBB, nside)
    C_QU = calc_cov_mat(-Wlm, Xlm, Xlm, Wlm, clEE, clBB, nside)
    C_UQ = calc_cov_mat(-Xlm, Wlm, Wlm, Xlm, clEE, clBB, nside)
    C_UU = calc_cov_mat(Xlm, Xlm, Wlm, Wlm, clEE, clBB, nside)
    
    C = zeros((2*npix,2*npix))
    C[1:npix, 1:npix] = C_QQ
    C[1:npix, npix+1:2*npix] = C_QU
    C[npix+1:2*npix, 1:npix] = C_UQ
    C[npix+1:2*npix, npix+1:2*npix] = C_UU
    
    return C
    
end

# CMB の signal 共分散行列を計算する 

function calc_cmb_cov_mat(clEE, clBB, nside)

    Wlm, Xlm = calc_WX(nside)
        
    npix = nside2npix(nside)

    C_QQ = calc_cov_mat(Wlm, Wlm, Xlm, Xlm, clEE, clBB, nside)
    C_QU = calc_cov_mat(-Wlm, Xlm, Xlm, Wlm, clEE, clBB, nside)
    C_UQ = calc_cov_mat(-Xlm, Wlm, Wlm, Xlm, clEE, clBB, nside)
    C_UU = calc_cov_mat(Xlm, Xlm, Wlm, Wlm, clEE, clBB, nside)
    
    C = zeros((2*npix,2*npix))
    C[1:npix, 1:npix] = C_QQ
    C[1:npix, npix+1:2*npix] = C_QU
    C[npix+1:2*npix, 1:npix] = C_UQ
    C[npix+1:2*npix, npix+1:2*npix] = C_UU
    
    return C

end

function calc_all_cov_mat(cov_mat_scal, cov_mat_tens, r)

    cov_mat = cov_mat_scal + r * cov_mat_tens
    
    return cov_mat

end