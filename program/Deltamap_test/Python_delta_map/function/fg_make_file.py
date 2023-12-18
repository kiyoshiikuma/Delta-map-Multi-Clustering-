import pysm3
import pysm3.units as u
import numpy as np
import sys, platform, os
import math

def fg_make_file(freq, nside):
    
    Synch = pysm3.Sky(nside, preset_strings=["s1"])
    Dust = pysm3.Sky(nside, preset_strings=["d1"])

    map_Synch = Synch.get_emission(freq * u.GHz)
    map_Dust = Dust.get_emission(freq * u.GHz)

    map_Synch = (map_Synch.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(freq*u.GHz))).value
    map_Dust = (map_Dust.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(freq*u.GHz))).value

    return map_Synch, map_Dust