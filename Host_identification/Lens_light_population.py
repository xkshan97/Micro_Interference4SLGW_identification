import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck18
import astropy.units as u

'''
#根据collett 2015的描述，引用了Hyde & Bernardi 2009年的Fundamental plane的文章
#使用方法：用log10_R_e随机生成一个log10_R_e，这个值是sersic轮廓的有效半径
#根据collett 2015的描述，各波段的log10_R_e一样，所以各波段只需要随机生成一次。
#然后根据生成的log10_R_e带入Apparent_magnitude_r/i/z程序，得到三个波段的星等。三个波段的星等通过选取不同的a,b,c来体现出差异。
def log10_R_e():
    #according to the average value and std in Collett reffed Hyde & Bernardi 2009
    #and assuming the effective radius Re is the same in all observed bands.
    return np.random.normal(0.39, 0.0552)

def Apparent_magitude_r(log10_R_e_, sigma_v, z_l):
    #r band 视星等，系数a，b，c来自于 Hyde & Bernardi 2009 表2，direct的方法。
    a = 1.1701
    b = 0.3029
    c = -8.0858
    mu_e = (log10_R_e_ - a * np.log10(sigma_v) - c) / b
    R_e_ = 10 ** log10_R_e_
    DA_l = Planck18.angular_diameter_distance(z_l)
    r_e = (R_e_ * u.kpc / DA_l).to(1) * (u.rad).to(u.arcsec)
    return mu_e - 5 * np.log10(r_e) - 2.5 * np.log10(2 * np.pi) + 10 * np.log10(1 + z_l)

def Apparent_magitude_i(log10_R_e_, sigma_v, z_l):
    #r band 视星等，系数a，b，c来自于 Hyde & Bernardi 2009 表2，direct的方法。
    a = 1.1990
    b = 0.3036
    c = -8.0481
    mu_e = (log10_R_e_ - a * np.log10(sigma_v) - c) / b
    R_e_ = 10 ** log10_R_e_
    DA_l = Planck18.angular_diameter_distance(z_l)
    r_e = (R_e_ * u.kpc / DA_l).to(1) * (u.rad).to(u.arcsec)
    return mu_e - 5 * np.log10(r_e) - 2.5 * np.log10(2 * np.pi) + 10 * np.log10(1 + z_l)

def Apparent_magitude_z(log10_R_e_, sigma_v, z_l):
    #r band 视星等，系数a，b，c来自于 Hyde & Bernardi 2009 表2，direct的方法。
    a = 1.2340
    b = 0.3139
    c = -8.2161
    mu_e = (log10_R_e_ - a * np.log10(sigma_v) - c) / b
    R_e_ = 10 ** log10_R_e_
    DA_l = Planck18.angular_diameter_distance(z_l)
    r_e = (R_e_ * u.kpc / DA_l).to(1) * (u.rad).to(u.arcsec)
    return mu_e - 5 * np.log10(r_e) - 2.5 * np.log10(2 * np.pi) + 10 * np.log10(1 + z_l)

'''




def Sample_u_R2mag_i(sigma_v, z_l):
    #参考的Goldstein 2019(SN)，但是均值是按照Wempe 2021修改的(GW)。
    #下面的程序针对的i band。
    #得到的是AB星等，但是我假设了没有颜色，所以AB星等就是V波段视星等MV。
    u_s = 19.40
    Q = 0.75
    sigma_u = 0.6
    R_s = 0.465
    sigma_R = 0.241
    V_s = 2.201
    sigma_V = 0.110
    rho_Ru = 0.753
    rho_Vu = -0.001
    rho_RV = 0.542
    u_s_c = u_s - Q * z_l
    V = np.log10(sigma_v)
    mean = (u_s_c + (V - V_s)/sigma_V * sigma_u * rho_Vu, R_s + (V - V_s)/sigma_V * sigma_R * rho_RV)
    cov = [[sigma_u ** 2 * (1 - rho_Vu ** 2), sigma_R * sigma_u * (rho_Ru - rho_RV * rho_Vu)], 
           [sigma_R * sigma_u * (rho_Ru - rho_RV * rho_Vu), sigma_R ** 2 * (1 - rho_RV ** 2)]]
    u_, R_ = np.random.multivariate_normal(mean, cov)
    R_e = 10 ** R_ * 1 / 0.7 * u.kpc
    Dl_l = Planck18.luminosity_distance(z_l)
    mag_i = u_ - 5 * np.log10((R_e / Dl_l).to(1) * (u.rad).to(u.arcsec)) - 2.5 * np.log10(2 * np.pi) + 10 * np.log10(1 + z_l)
    return mag_i, R_e.value


Generate_SLGW_Galaxy_mag_R_e = 'False'
if Generate_SLGW_Galaxy_mag_R_e:
    np.random.seed(1000000)
    lens_z = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')
    Sigma_v_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/sigma_v.csv', delimiter=',')

    index_set_selected = [170, 120, 46]

    for index in index_set_selected:
        z_l = lens_z[index]
        sigma_v = Sigma_v_set[index]
        Sample_u_R2mag_i_Res_tmp = Sample_u_R2mag_i(sigma_v, z_l)
        Sample_u_R2mag_i_Res = [Sample_u_R2mag_i_Res_tmp[0].value, Sample_u_R2mag_i_Res_tmp[1]]
        np.savetxt('./SampleResult_GWGalaxy/GW_Lens_Galaxy_' + str(index) + '.csv', Sample_u_R2mag_i_Res, delimiter=',') 