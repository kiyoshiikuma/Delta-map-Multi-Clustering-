# maskをしたmap返す

using Healpix

function masked_data_calc(map)

    mask_map = Healpix.readMapFromFITS("../../mask_p06_Nside4.v2.fits", 1, Float64);

    unmasked_indices = findall(x -> x == 1, mask_map);
    masked_indices = findall(x -> x == 0, mask_map);

    # maskされていない配列を前に持ってくる

    masked_map = [map[i] for i in 1:length(map) if mask_map[i] == 1];

    return masked_map
end

function masked_map_calc(map)

    mask_map = Healpix.readMapFromFITS("../../mask_p06_Nside4.v2.fits", 1, Float64);

    mask_map_I = map[1, :] .* mask_map;
    mask_map_Q = map[2, :] .* mask_map;
    mask_map_U = map[3, :] .* mask_map;

    return [mask_map_I, mask_map_Q, mask_map_U]

end