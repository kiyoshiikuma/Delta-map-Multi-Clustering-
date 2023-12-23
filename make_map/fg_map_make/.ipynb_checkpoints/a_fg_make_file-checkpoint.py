import pysm3
import pysm3.units as u
import numpy as np
import sys, platform, os
import math
import healpy


def fg_make_file(freq, nside):
    
    Synch = pysm3.Sky(nside, preset_strings=["s1"])
    Dust = pysm3.Sky(nside, preset_strings=["d1"])

    map_Synch = Synch.get_emission(freq * u.GHz)
    map_Dust = Dust.get_emission(freq * u.GHz)

    map_Synch = (map_Synch.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(freq*u.GHz))).value
    map_Dust = (map_Dust.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(freq*u.GHz))).value

    # file name
    dir = "fg_map_data/"
    nside_n = "nside_"
    
    fg_synch = "Synch_"
    GHz = "_GHz_"
    
    Synch_name = f"{dir}{fg_synch}{freq}{GHz}{nside_n}{nside}"

    fg_dust = "Dust_"
    
    Dust_name = f"{dir}{fg_dust}{freq}{GHz}{nside_n}{nside}"
    
    np.save(Synch_name, map_Synch)
    np.save(Dust_name, map_Dust)

    healpy.fitsfunc.write_map(Synch_name, map_Synch, overwrite=True)
    healpy.fitsfunc.write_map(Dust_name, map_Dust, overwrite=True)

    return map_Synch, map_Dust

