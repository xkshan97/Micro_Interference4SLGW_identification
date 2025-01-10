import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from tqdm import tqdm

def Mstar_V_d(log_v):
    a = 2.5
    b = 5
    logM = a * log_v + b
    return logM

Sigma_v_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/sigma_v.csv', delimiter=',')
lens_z = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')

M_star_set = 10**Mstar_V_d(np.log10(Sigma_v_set))

hdu=fits.open('JWST_catalogue/JADES_SF_mock_r1_v1.2.fits')
name = hdu[1].data.names
galaxy_redshift = hdu[1].data['redshift']
galaxy_SFR_10 = 10**hdu[1].data['SFR_10']
galaxy_SFR_100 = 10**hdu[1].data['SFR_100']
metallicity = 10**hdu[1].data['metallicity'] * 0.01524
mstar = 10**hdu[1].data['mStar']

metallicity_set = []
galaxy_redshift_set = []
mstar_set = []
galaxy_SFR_set = []
for i in tqdm(range(len(lens_z))):
    index_tmp = np.where((galaxy_redshift<lens_z[i]+1 + 0.01)&(galaxy_redshift>lens_z[i]+1-0.01)&(np.log10(mstar)>np.log10(M_star_set[i])-0.01)&(np.log10(mstar)>np.log10(M_star_set[i])+0.01))
    galaxy_redshift_set.append(np.mean(galaxy_redshift[index_tmp]))
    metallicity_set.append(np.mean(np.log10(metallicity[index_tmp])))
    mstar_set.append(np.mean(np.log10(mstar[index_tmp])))
    galaxy_SFR_set.append(np.mean(galaxy_SFR_10[index_tmp]))
    
metallicity_set = np.array(metallicity_set)
galaxy_redshift_set = np.array(galaxy_redshift_set)
mstar_set = np.array(mstar_set) 
galaxy_SFR_set = np.array(galaxy_SFR_set)

plt.hist(10**metallicity_set)

plt.hist(galaxy_redshift_set - lens_z)

plt.hist(mstar_set - np.log10(M_star_set))

plt.hist2d(np.log10(metallicity_0p5_1p5), np.log10(mstar_0p5_1p5))
plt.xlabel('$\log_{10} Z$')
plt.ylabel('$\log_{10} M_\mathrm{star}$')

plt.hist2d(np.log10(metallicity_0p5_1p5), np.log10(galaxy_SFR_10_0p5_1p5))
plt.xlabel('$\log_{10} Z$')
plt.ylabel('$\log_{10} \mathrm{SFR}$[M$_\odot\ \mathrm{yr}^{-1}$]')

plt.hist(np.log10(galaxy_SFR_10_0p5_1p5))


plt.scatter(galaxy_redshift, metallicity) 
plt.semilogy()