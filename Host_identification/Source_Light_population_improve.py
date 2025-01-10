import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from tqdm import tqdm
from astropy.cosmology import Planck18
import astropy.units as u


hdu=fits.open('../Paper4_Source_Rec/JWST_catalogue/JADES_SF_mock_r1_v1.2.fits')
name = hdu[1].data.names
galaxy_redshift = hdu[1].data['redshift']
galaxy_SFR_10 = 10**hdu[1].data['SFR_10']
galaxy_SFR_100 = 10**hdu[1].data['SFR_100']
galaxy_SFR_50 = (galaxy_SFR_10 + galaxy_SFR_100) / 2
SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')

index_set_selected = [170, 120, 46]

for index in index_set_selected:
    z_s = SampleParam[0][int(AcceptLensIndex[index])]
    dis_abs = np.abs(galaxy_redshift - z_s) #星系红移与此引力波红移之间的距离。
    

    #检验一下0.01误差范围内的星等分布
    index_range = np.where(dis_abs < 0.01)[0]
    mag_range = []
    JWST_mag_range = []
    axis_ratio_range = []
    position_angle_range = []
    R_e_range = []
    sersic_n_range = []
    SFR_50_range = []
    SFR_100_range = []
    selected_source_z_s = []

    for i in range(len(index_range)):
        index_tmp = index_range[i]
        selected_source_z_s.append(galaxy_redshift[index_tmp])
        
        HST_F606W_fnu_tmp = hdu[1].data['HST_F606W_fnu'][index_tmp]
        NRC_F200W_fnu_tmp = hdu[1].data['NRC_F200W_fnu'][index_tmp] 


        HST_flux_tmp = HST_F606W_fnu_tmp * 1.0000000000000002e-32 #转成erg/s/cm^2/Hz
        apparent_mag_tmp = -2.5 * np.log10(HST_flux_tmp) - 48.6 
        
        
        JWST_flux_tmp = NRC_F200W_fnu_tmp * 1.0000000000000002e-32 #转成erg/s/cm^2/Hz
        JWST_apparent_mag_tmp = -2.5 * np.log10(JWST_flux_tmp) - 48.6
            
            
        axis_ratio_tmp = hdu[1].data['axis_ratio'][index_tmp]
        position_angle_tmp = hdu[1].data['position_angle'][index_tmp]
        R_e_tmp = hdu[1].data['Re_circ'][index_tmp]
        sersic_n_tmp = hdu[1].data['sersic_n'][index_tmp]
        
        
        mag_range.append(apparent_mag_tmp)
        JWST_mag_range.append(JWST_apparent_mag_tmp)

        axis_ratio_range.append(axis_ratio_tmp)
        position_angle_range.append(position_angle_tmp)
        R_e_range.append(R_e_tmp)
        sersic_n_range.append(sersic_n_tmp)
        SFR_50_range.append(galaxy_SFR_50[index_tmp])
        SFR_100_range.append(galaxy_SFR_100[index_tmp])

    mag_weight_by_SFR_50 = []
    mag_weight_by_SFR_100 = []
    JWST_mag_weight_by_SFR_50 = []
    JWST_mag_weight_by_SFR_100 = []
    for i in tqdm(range(len(index_range))):
        weight_tmp_50 = SFR_50_range[i] / np.min(SFR_50_range)
        for j in range(int(weight_tmp_50)):
            mag_weight_by_SFR_50.append(mag_range[i])
            JWST_mag_weight_by_SFR_50.append(JWST_mag_range[i])
        weight_tmp_100 = SFR_100_range[i] / np.min(SFR_100_range)
        for j in range(int(weight_tmp_100)):
            mag_weight_by_SFR_100.append(mag_range[i])
            JWST_mag_weight_by_SFR_100.append(JWST_mag_range[i])
            
    plt.hist(mag_range, density=True, bins=20, label='unweight')
    plt.hist(JWST_mag_range, density=True, bins=20, label='JWST unweight')
     
    plt.hist(mag_weight_by_SFR_50, density=True, bins = 20, label='weighted by 50 Myr SFR')
    plt.hist(JWST_mag_weight_by_SFR_50, density=True, bins = 20, label='JWST weighted by 50 Myr SFR')
    # plt.hist(mag_weight_by_SFR_100, density=True, bins = 20, label='weighted by 100 Myr SFR')
    plt.xlabel('mag')
    plt.ylabel('pdf')
    plt.grid()
    plt.legend()
    plt.savefig('./SampleResult_GWGalaxy/mag_Res_' + str(index) +'.png', dpi=450)
    plt.close()
    np.savetxt('./SampleResult_GWGalaxy/mag_range_' + str(index) + '.csv', mag_range, delimiter=',') 
    np.savetxt('./SampleResult_GWGalaxy/mag_range_weight_by_SFR_50' + str(index) + '.csv', mag_weight_by_SFR_50, delimiter=',')
    np.savetxt('./SampleResult_GWGalaxy/mag_range_weight_by_SFR_100' + str(index) + '.csv', mag_weight_by_SFR_100, delimiter=',')
    np.savetxt('./SampleResult_GWGalaxy/JWST_mag_range_' + str(index) + '.csv', JWST_mag_range, delimiter=',') 
    np.savetxt('./SampleResult_GWGalaxy/JWST_mag_range_weight_by_SFR_50' + str(index) + '.csv', JWST_mag_weight_by_SFR_50, delimiter=',')
    np.savetxt('./SampleResult_GWGalaxy/JWST_mag_range_weight_by_SFR_100' + str(index) + '.csv', JWST_mag_weight_by_SFR_100, delimiter=',')
     
    

    SFR_50_range = np.array(SFR_50_range)
    index_gtr_1 = np.where(SFR_50_range >= 1)[0]
    
    mag_Res = []
    JWST_mag_Res = []
    for index_tmp in index_gtr_1:
        mag_tmp = mag_range[index_tmp]
        JWST_mag_tmp = JWST_mag_range[index_tmp]
        axis_ratio_tmp = axis_ratio_range[index_tmp]
        position_angle_tmp = position_angle_range[index_tmp]
        angular_diameter_distance_to_source = Planck18.angular_diameter_distance(selected_source_z_s[index_tmp]).to(u.kpc).value
        R_e_tmp = R_e_range[index_tmp] / angular_diameter_distance_to_source * u.rad.to(u.arcsec)
        sersic_n_tmp = sersic_n_range[index_tmp]
        SFR_50_tmp = SFR_50_range[index_tmp]
        mag_Res.append([mag_tmp, axis_ratio_tmp, position_angle_tmp, R_e_tmp, sersic_n_tmp, SFR_50_tmp])
        JWST_mag_Res.append([JWST_mag_tmp, axis_ratio_tmp, position_angle_tmp, R_e_tmp, sersic_n_tmp, SFR_50_tmp])
    np.savetxt('./SampleResult_GWGalaxy/mag_Res_' + str(index) + '.csv', mag_Res, delimiter=',')
    np.savetxt('./SampleResult_GWGalaxy/JWST_mag_Res_' + str(index) + '.csv', JWST_mag_Res, delimiter=',')
    
'''
经过验证，通过谱计算出的星等，跟直接用三个HST波段算出的差不多。所以就直接用HST的三个波段吧。
'''
# hdu_spectrum=fits.open('JWST_catalogue/JADES_SF_mock_r1_v1.2_spec_5A_30um_z_1p5_2.fits')
# spectrum = hdu_spectrum['FULL SED'].data #erg / s / cm ^ 2 / A
# spectrum_WL = hdu_spectrum['FULL SED WL'].data
# spectrum_red = hdu_spectrum['OBJECT PROPERTIES'].data
# for i in range(len(spectrum_red)):
#     if spectrum_red[i][0] == index_min:
#         index_in_spectrum = i
#         print(index_in_spectrum)

# spectrum_WL_red = spectrum_WL * (1 + spectrum_red[index_in_spectrum][1])

# index_in_filter = np.where((spectrum_WL_red < 9300)&(spectrum_WL_red > 5000))
# plt.plot(spectrum_WL, spectrum[index_in_spectrum])
# plt.plot(spectrum_WL_red, spectrum[index_in_spectrum], 'r--')
# plt.plot(spectrum_WL_red[index_in_filter], spectrum[index_in_spectrum][index_in_filter], 'k')
# plt.xlim(0, 10000)
# np.diff(spectrum_WL_red[index_in_filter])
# Sum_flux = 0.8 * np.sum(spectrum[index_in_spectrum][index_in_filter] * np.diff(spectrum_WL_red[index_in_filter])[0])
# f_start = 3 * 10 ** (8) / (5000 * 10 ** (-10))
# f_end = 3 * 10 ** (8) / (9300 * 10 ** (-10))
# flux_per_hz = Sum_flux / (f_start - f_end)
# apparent_mag = -2.5 * np.log10(flux_per_hz) - 48.6 