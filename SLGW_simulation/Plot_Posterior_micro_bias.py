import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import matplotlib.lines as mlines
from tqdm import tqdm
import json
import corner
from astropy.cosmology import Planck15
from getdist import plots, MCSamples
import getdist
import matplotlib.lines as mlines
import IPython

#参数估计结果

SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')


AcceptLensIndex = np.loadtxt('./SampleResult/AcceptLensIndex.csv', delimiter=',')
imagenum = np.loadtxt('./SampleResult/imagenumber.csv', delimiter=',')
SNR_network = np.loadtxt('./Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')
kappa = np.loadtxt('./SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('./SampleResult/gamma.csv', delimiter=',')
timedelay = np.loadtxt('./SampleResult/timedelayDays.csv', delimiter=',')
lens_z_set = np.loadtxt('./SampleResult/lens_z.csv', delimiter=',')
magnification = np.loadtxt('./SampleResult/magnification.csv', delimiter=',')
Total_inject_time = np.loadtxt('./Lensed_SampleResult_Micro/Total_inject_time.csv', delimiter=',')

trigger_time_H1 = np.loadtxt('./Lensed_SampleResult/trigger_time_H1.csv', delimiter=',')
trigger_time_L1 = np.loadtxt('./Lensed_SampleResult/trigger_time_L1.csv', delimiter=',')
trigger_time_V1 = np.loadtxt('./Lensed_SampleResult/trigger_time_V1.csv', delimiter=',')


dt_h1_from_earth_center = np.loadtxt('./Lensed_SampleResult/dt_h1_from_earth_center.csv', delimiter=',')
dt_l1_from_earth_center = np.loadtxt('./Lensed_SampleResult/dt_l1_from_earth_center.csv', delimiter=',')      
dt_v1_from_earth_center = np.loadtxt('./Lensed_SampleResult/dt_v1_from_earth_center.csv', delimiter=',')      

sample_rate = 4096

IndexInAccept = np.loadtxt('../Paper5_cWB_Modify/Total_Data_start_time_and_snr_and_match_and_outdir_name/IndexInAccept.csv', delimiter=',')
SNR_selected = np.loadtxt('../Paper5_cWB_Modify/Total_Data_start_time_and_snr_and_match_and_outdir_name/SNR_selected.csv', delimiter=',')
labels=['q', '\mathcal{M}', 'a_1', 'a_2', '\iota', 'POL', '\mathrm{ra}', '\mathrm{\sin{(dec)}}']
names = labels

for i_ in [120]:
    i_ = int(i_)
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = int(imagenum[i_])
    save_index = 0
    IndexSumimage = int(np.sum(imagenum[0:i_]))
    

    for j_ in range(imag_num):
        
        if SNR_network[IndexSumimage] > 12:
            
            if SNR_network[IndexSumimage] in SNR_selected:
                td_second = timedelay[IndexSumimage] * 24 * 3600
                SampleParam9_td = SampleParam[9][i] + td_second
                # print(SNR_network[IndexSumimage])
                sqrt_mag_abs = np.sqrt(np.abs(magnification[IndexSumimage]))
                if 1 - kappa[IndexSumimage] - gamma[IndexSumimage] > 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] > 0:    
                    """================================="""
                    Lens_Type = "Minimum"
                elif (1 - kappa[IndexSumimage] - gamma[IndexSumimage]) * (1 - kappa[IndexSumimage] + gamma[IndexSumimage]) < 0:
                    Lens_Type = "Saddle"
                    
                try:
                    outdir = outdir = outdir = 'outdir_Lensed_' + Lens_Type + '_with_micro/outdir' + str(i) + "_" + str(save_index)
                    
                    #=======================
                    #injection parameters
                    luminosity_distance = Planck15.luminosity_distance(SampleParam[0][i]).value
                    Mz = (SampleParam[0][i] + 1) * (SampleParam[1][i] * SampleParam[2][i]) ** (3/5) / (SampleParam[1][i] + SampleParam[2][i]) ** (1/5)
                    injection_parameters = np.array([SampleParam[2][i]/SampleParam[1][i], Mz, SampleParam[3][i], SampleParam[4][i],SampleParam[5][i], SampleParam[6][i],SampleParam[7][i], np.sin(SampleParam[8][i])])
                    #注意geocent_time改成了H1_time。
                    
                    
                    
                    


                    file_post = outdir + '/H1_L1_V1_result.json'

                    with open(file_post,"r") as f:
                        data_post = json.load(f)

                    
                    #最大似然点
                    j = np.where(data_post['posterior']['content']['log_likelihood'] == np.max(data_post['posterior']['content']['log_likelihood']))[0][0]
                    q_max = data_post['posterior']['content']['mass_ratio'][j]
                    Mz_o_max = data_post['posterior']['content']['chirp_mass'][j]
                    luminosity_dis_max = data_post['posterior']['content']['luminosity_distance'][j]
                    dec_max = data_post['posterior']['content']['dec'][j]
                    ra_max = data_post['posterior']['content']['ra'][j]
                    theta_jn_max = data_post['posterior']['content']['theta_jn'][j]
                    psi_max = data_post['posterior']['content']['psi'][j]
                    phase_max = data_post['posterior']['content']['phase'][j]
                    a_1_max = data_post['posterior']['content']['a_1'][j]
                    a_2_max = data_post['posterior']['content']['a_2'][j]
                    tilt_1_max = data_post['posterior']['content']['tilt_1'][j]
                    tilt_2_max = data_post['posterior']['content']['tilt_2'][j]
                    phi_12_max = data_post['posterior']['content']['phi_12'][j]
                    phi_jl_max = data_post['posterior']['content']['phi_jl'][j]
                    geocent_time_max = data_post['posterior']['content']['geocent_time'][j]
                    max_likelihood_parameters = np.array([q_max, Mz_o_max, a_1_max, a_1_max, theta_jn_max, psi_max, ra_max, dec_max, luminosity_dis_max, geocent_time_max]) 
                    SNR_tmp = np.sqrt(data_post['posterior']['content']['H1_matched_filter_snr'][j]['real']**2 + data_post['posterior']['content']['L1_matched_filter_snr'][j]['real']**2 + data_post['posterior']['content']['V1_matched_filter_snr'][j]['real']**2)
                    
                    
                    
                    nest_samples = data_post['posterior']['content']
                    mass_ratio_res = nest_samples['mass_ratio']
                    chirp_mass_res = nest_samples['chirp_mass']
                    ra_res = nest_samples['ra']
                    dec_res = nest_samples['dec']
                    theta_jn_res = nest_samples['theta_jn']
                    psi_res = nest_samples['psi']
                    a_1_res = nest_samples['a_1']
                    a_2_res = nest_samples['a_2']
                    luminosity_distance_res = nest_samples['luminosity_distance']
                    geocent_time = nest_samples['geocent_time']


                    samples_res = []
                    for _ in tqdm(range(len(chirp_mass_res))):
                        
                        # m1 = chirp_mass[i]*(1 + mass_ratio[i])**(1/5)/mass_ratio[i]**(3/5)
                        # m2 = mass_ratio[i]*m1
                        samples_res.append([mass_ratio_res[_], chirp_mass_res[_], a_1_res[_], a_2_res[_], theta_jn_res[_], psi_res[_], ra_res[_], np.sin(dec_res[_])]) 

                    samples_res = np.array(samples_res)
                    
                    samples_final = MCSamples(samples=np.array(samples_res), names = names, labels = labels)
                    
                    samples_final.updateSettings({'contours': [0.6826, 0.9544]})

                    g = plots.get_subplot_plotter()
                    g.settings.legend_fontsize = 10
                    g.settings.axes_fontsize = 15
                    g.settings.axes_labelsize = 15
                    g.settings.axis_tick_x_rotation = 45
                    g.triangle_plot(samples_final, colors=['grey'],filled=True,legend_loc='upper right',
                                line_args=[{'color':'grey'}], markers={'q': injection_parameters[0], '\mathcal{M}': injection_parameters[1], 'a_1': injection_parameters[2], 'a_2': injection_parameters[3], '\iota': injection_parameters[4], 'POL': injection_parameters[5], '\mathrm{ra}': injection_parameters[6], '\mathrm{\sin{(dec)}}': injection_parameters[7]}
                            ,marker_args={'color':'black','lw':20,'lw':1} ,contour_lws=[1])
                    plt.savefig('Final_Bias_show/corner_' + str(i) + str(save_index) +'.png', dpi=450)
                    # plt.close()

                except FileNotFoundError:
                    pass
            save_index += 1
        IndexSumimage += 1



