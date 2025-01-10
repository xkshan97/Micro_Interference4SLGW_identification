import numpy as np
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import healpy as hp
from healpy.newvisufunc import projview, newprojplot

skypro = hp.read_map('cWB_micro/data/1025263199_3580_0/L1H1V1_1025265244.000_1025265244.000_1025265244.000/skyprobcc.fits')

projview(
    skypro,
    coord=["G", "C"],
    graticule=True,
    graticule_labels=True,
    unit=r"cbar label",
    xlabel="RA",
    ylabel="DEC",
    cb_orientation="vertical",
    latitude_grid_spacing=18,
    projection_type="hammer"
)
