import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from tqdm import tqdm
from astropy.cosmology import Planck18
import astropy.units as u


hdu=fits.open('JWST_catalogue/JADES_SF_mock_r1_v1.2.fits')
name = hdu[1].data.names
galaxy_redshift = hdu[1].data['redshift']
galaxy_SFR_10 = 10**hdu[1].data['SFR_10']
galaxy_SFR_100 = 10**hdu[1].data['SFR_100']
HST_F606W_fnu = hdu[1].data['HST_F606W_fnu']
HST_flux = HST_F606W_fnu * 1.0000000000000002e-32 #转成erg/s/cm^2/Hz
apparent_mag = -2.5 * np.log10(HST_flux) - 48.6 
apparent_mag[np.isinf(apparent_mag)] = 100
index_plot = np.where(apparent_mag < 26)[0]

RA = hdu[1].data['RA']
DEC = hdu[1].data['DEC']

plt.scatter(RA[index_plot], DEC[index_plot], s=np.exp(26 / apparent_mag[index_plot] - 1))