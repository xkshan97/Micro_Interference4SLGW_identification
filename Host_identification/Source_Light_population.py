import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from tqdm import tqdm



hdu=fits.open('JWST_catalogue/JADES_SF_mock_r1_v1.2.fits')
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
    axis_ratio_range = []
    position_angle_range = []
    R_e_range = []
    sersic_n_range = []
    SFR_50_range = []
    SFR_100_range = []

    for i in range(len(index_range)):
        index_tmp = index_range[i]
        HST_F435W_fnu_tmp = hdu[1].data['HST_F435W_fnu'][index_tmp]
        HST_F606W_fnu_tmp = hdu[1].data['HST_F606W_fnu'][index_tmp]
        HST_F775W_fnu_tmp = hdu[1].data['HST_F775W_fnu'][index_tmp]

        Ave_flux_tmp = (HST_F435W_fnu_tmp + HST_F606W_fnu_tmp + HST_F775W_fnu_tmp) / 3 * 1.0000000000000002e-32 #转成erg/s/cm^2/Hz
        apparent_mag_tmp = -2.5 * np.log10(Ave_flux_tmp) - 48.6 
        axis_ratio_tmp = hdu[1].data['axis_ratio'][index_tmp]
        position_angle_tmp = hdu[1].data['position_angle'][index_tmp]
        R_e_tmp = hdu[1].data['Re_circ'][index_tmp]
        sersic_n_tmp = hdu[1].data['sersic_n'][index_tmp]
        
        mag_range.append(apparent_mag_tmp)
        axis_ratio_range.append(axis_ratio_tmp)
        position_angle_range.append(position_angle_tmp)
        R_e_range.append(R_e_tmp)
        sersic_n_range.append(sersic_n_tmp)
        SFR_50_range.append(galaxy_SFR_50[index_tmp])
        SFR_100_range.append(galaxy_SFR_100[index_tmp])

    mag_weight_by_SFR_50 = []
    mag_weight_by_SFR_100 = []
    for i in tqdm(range(len(index_range))):
        weight_tmp_50 = SFR_50_range[i] / np.min(SFR_50_range)
        for j in range(int(weight_tmp_50)):
            mag_weight_by_SFR_50.append(mag_range[i])
        weight_tmp_100 = SFR_100_range[i] / np.min(SFR_100_range)
        for j in range(int(weight_tmp_100)):
            mag_weight_by_SFR_100.append(mag_range[i])
            
    plt.hist(mag_range, density=True, bins=20, label='unweight')
    plt.hist(mag_weight_by_SFR_50, density=True, bins = 20, label='weighted by 50 Myr SFR')
    plt.hist(mag_weight_by_SFR_100, density=True, bins = 20, label='weighted by 100 Myr SFR')
    plt.xlabel('mag')
    plt.ylabel('pdf')
    plt.grid()
    plt.legend()
    plt.savefig('./SampleResult_GWGalaxy/mag_pdf_' + str(index) + '.png', dpi=450)
    plt.close()
    np.savetxt('./SampleResult_GWGalaxy/mag_range_' + str(index) + '.csv', mag_range, delimiter=',') 
    np.savetxt('./SampleResult_GWGalaxy/mag_range_weight_by_SFR_50' + str(index) + '.csv', mag_weight_by_SFR_50, delimiter=',')
    np.savetxt('./SampleResult_GWGalaxy/mag_range_weight_by_SFR_100' + str(index) + '.csv', mag_weight_by_SFR_100, delimiter=',')
    
    index_min_mag = np.where(np.array(mag_range) == np.min(mag_range))[0][0]
    mag_min = mag_range[index_min_mag]
    axis_ratio_min = axis_ratio_range[index_min_mag]
    position_angle_min = position_angle_range[index_min_mag]
    R_e_min = R_e_range[index_min_mag]
    sersic_n_min = sersic_n_range[index_min_mag]
    mag_min_Res = [mag_min, axis_ratio_min, position_angle_min, R_e_min, sersic_n_min, np.sum(np.array(mag_weight_by_SFR_50)<25.5) / len(mag_weight_by_SFR_50), np.sum(np.array(mag_weight_by_SFR_100)<25.5) / len(mag_weight_by_SFR_100)]
    np.savetxt('./SampleResult_GWGalaxy/mag_min_Res_' + str(index) + '.csv', mag_min_Res, delimiter=',')
    print('min mag: ', mag_min)
    print('min axis ratio: ', axis_ratio_min)
    print('min position angle: ', position_angle_min)
    print('min R_e: ', R_e_min)
    print('min sersic n: ', sersic_n_min)

    try:
        index_max_mag = np.where((np.array(mag_range)<25.5)&(np.array(mag_range)>25.4))[0][0]
        mag_max = mag_range[index_max_mag]
        axis_ratio_max = axis_ratio_range[index_max_mag]
        position_angle_max = position_angle_range[index_max_mag]
        R_e_max = R_e_range[index_max_mag]
        sersic_n_max = sersic_n_range[index_max_mag]
        mag_max_Res = [mag_max, axis_ratio_max, position_angle_max, R_e_max, sersic_n_max, np.sum(np.array(mag_weight_by_SFR_50)<25.5) / len(mag_weight_by_SFR_50), np.sum(np.array(mag_weight_by_SFR_100)<25.5) / len(mag_weight_by_SFR_100)]
        np.savetxt('./SampleResult_GWGalaxy/mag_max_Res_' + str(index) + '.csv', mag_max_Res, delimiter=',')
        print('max mag: ', mag_max)
        print('max axis ratio: ', axis_ratio_max)
        print('max position angle: ', position_angle_max)
        print('max R_e: ', R_e_max)
        print('max sersic n: ', sersic_n_max)
    except IndexError:
        print('minimum magnitude greater than 25.5')
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