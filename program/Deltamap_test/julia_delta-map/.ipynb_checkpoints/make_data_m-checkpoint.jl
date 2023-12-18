import Healpix

function make_data_m(cmb, freq_band, which_model, nside)
    
    band = length(freq_band)

    N_pix = Healpix.nside2npix(nside)

    data_m = zeros(Float64, 2 * N_pix * band)

    for (ii, freq) in enumerate(freq_band)

        # file name
        dir = "/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/fg_map_data/"
        nside_name = "nside_"
        nside_n = string(nside)

        fg_synch = "Synch_"
        GHz = "_GHz_"

        # Synchrotron name
        Synch_name = string(dir, fg_synch, freq, GHz, nside_name, nside_n)

        fg_dust = "Dust_"

        # Dust name
        Dust_name = string(dir, fg_dust, freq, GHz, nside_name, nside_n)

        #==
        # file name
        dir = "/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_fg_rm/program/Deltamap_test/julia_delta-map/fg_map_data/"
        #dir = "fg_map_data/"
        nside_name = "nside_"
        nside_n = str(nside)
        #npy =".npy"

        fg_synch = "Synch_"
        GHz = "_GHz_"

        # Synchrotron name
        Synch_name = f"{dir}{fg_synch}{freq}{GHz}{nside_name}{nside_n}{npy}"

        fg_dust = "Dust_"
        GHz = "_GHz"

        # Dust name
        Dust_name = f"{dir}{fg_dust}{freq}{GHz}{nside_name}{nside_n}{npy}"
        ==#

        
        #synch_data = np.load(Synch_name)
        synch_data_I = Healpix.readMapFromFITS(Synch_name, 1, Float64)
        synch_data_Q = Healpix.readMapFromFITS(Synch_name, 2, Float64)
        synch_data_U = Healpix.readMapFromFITS(Synch_name, 3, Float64)
        
        
        #dust_data = np.load(Dust_name)
        Dust_data_I = Healpix.readMapFromFITS(Dust_name, 1, Float64)
        Dust_data_Q = Healpix.readMapFromFITS(Dust_name, 2, Float64)
        Dust_data_U = Healpix.readMapFromFITS(Dust_name, 3, Float64)

        
        if which_model == "s1"
            
            m_Q = synch_data_Q + cmb[2]
            m_U = synch_data_U + cmb[3]
            
            data_m[Int((2*(ii-1))*N_pix)+1 : Int((2*(ii-1)+1)*N_pix)] = m_Q
            data_m[Int((2*(ii-1)+1)*N_pix) : Int(((2*(ii-1)+2))*N_pix)] = m_U
            
        elseif which_model == "d1"

            # data m[I, Q, U]
            m_Q = Dust_data_Q + cmb[2]
            m_U = Dust_data_U + cmb[3]

            data_m[Int((2*ii)*N_pix) : Int((2*ii+1)*N_pix)] = m_Q
            data_m[Int((2*ii+1)*N_pix) : Int(((2*ii+2))*N_pix)] = m_U     

        elseif which_model == "d1 and s1"

            # data m[I, Q, U]
            m_Q = synch_data_Q + Dust_data_Q + cmb[2]
            m_U = synch_data_U + Dust_data_U + cmb[3]

            data_m[Int((2*(ii-1))*N_pix)+1 : Int((2*(ii-1)+1)*N_pix)] = m_Q
            data_m[Int((2*(ii-1)+1)*N_pix)+1 : Int(((2*(ii-1)+2))*N_pix)] = m_U

        else
            print("正しいモデル入れてね")

        end
        
    end

    return data_m'

end