import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import matplotlib.lines as mlines
from tqdm import tqdm
import json
import corner
from astropy.cosmology import Planck15

#参数估计结果
spin1x = np.loadtxt('./Sim_GW_Data_precessing/spin1x.csv', delimiter=',')      
spin2x = np.loadtxt('./Sim_GW_Data_precessing/spin2x.csv', delimiter=',')      
spin1y = np.loadtxt('./Sim_GW_Data_precessing/spin1y.csv', delimiter=',')      
spin2y = np.loadtxt('./Sim_GW_Data_precessing/spin2y.csv', delimiter=',')      
spin1z = np.loadtxt('./Sim_GW_Data_precessing/spin1z.csv', delimiter=',')      
spin2z = np.loadtxt('./Sim_GW_Data_precessing/spin2z.csv', delimiter=',')

SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')

trigger_time_H1 = np.loadtxt('./Sim_GW_Data_precessing/trigger_time_H1.csv', delimiter=',')
trigger_time_L1 = np.loadtxt('./Sim_GW_Data_precessing/trigger_time_L1.csv', delimiter=',')


dt_h1_from_earth_center = np.loadtxt('./Sim_GW_Data_precessing/dt_h1_from_earth_center.csv', delimiter=',')
dt_l1_from_earth_center = np.loadtxt('./Sim_GW_Data_precessing/dt_l1_from_earth_center.csv', delimiter=',')      
dt_v1_from_earth_center = np.loadtxt('./Sim_GW_Data_precessing/dt_v1_from_earth_center.csv', delimiter=',')      

sample_rate = 4096

SNR_select_index = np.loadtxt('../Paper5_cWB_Modify/Data_start_time_and_snr_and_match_and_outdir_name/SNR_select_index_unlens.csv', delimiter=',')

for i in SNR_select_index:
    i = int(i)
    try:
        outdir = 'outdir_unlensed_precessing/outdir' + str(i)
        
        #=======================
        #injection parameters
        luminosity_distance = Planck15.luminosity_distance(SampleParam[0][i]).value
        Mz = (SampleParam[0][i] + 1) * (SampleParam[1][i] * SampleParam[2][i]) ** (3/5) / (SampleParam[1][i] + SampleParam[2][i]) ** (1/5)
        injection_parameters = np.array([SampleParam[2][i]/SampleParam[1][i], Mz, SampleParam[3][i], SampleParam[4][i], spin1x[i], spin1y[i], spin1z[i], spin2x[i], spin2y[i], spin2z[i], SampleParam[5][i], SampleParam[6][i],SampleParam[7][i], SampleParam[8][i], luminosity_distance, trigger_time_H1[i] - dt_h1_from_earth_center[i]])
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
        spin_1x_max = data_post['posterior']['content']['spin_1x'][j]
        spin_1y_max = data_post['posterior']['content']['spin_1y'][j]
        spin_1z_max = data_post['posterior']['content']['spin_1z'][j]
        spin_2x_max = data_post['posterior']['content']['spin_2x'][j]
        spin_2y_max = data_post['posterior']['content']['spin_2y'][j]
        spin_2z_max = data_post['posterior']['content']['spin_2z'][j]
        tilt_1_max = data_post['posterior']['content']['tilt_1'][j]
        tilt_2_max = data_post['posterior']['content']['tilt_2'][j]
        phi_12_max = data_post['posterior']['content']['phi_12'][j]
        phi_jl_max = data_post['posterior']['content']['phi_jl'][j]
        geocent_time_max = data_post['posterior']['content']['geocent_time'][j]
        iota_max = data_post['posterior']['content']['iota'][j]
        max_likelihood_parameters = np.array([q_max, Mz_o_max, a_1_max, a_1_max, spin_1x_max, spin_1y_max, spin_1z_max, spin_2x_max, spin_2y_max, spin_2z_max, iota_max, psi_max, ra_max, dec_max, luminosity_dis_max, geocent_time_max]) 
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
        spin_1x_res = nest_samples['spin_1x']
        spin_1y_res = nest_samples['spin_1y']
        spin_1z_res = nest_samples['spin_1z']
        spin_2x_res = nest_samples['spin_2x']
        spin_2y_res = nest_samples['spin_2y']
        spin_2z_res = nest_samples['spin_2z']
        luminosity_distance_res = nest_samples['luminosity_distance']
        geocent_time = nest_samples['geocent_time']
        iota_res = nest_samples['iota']


        samples_res = []
        for _ in tqdm(range(len(chirp_mass_res))):
            
            # m1 = chirp_mass[i]*(1 + mass_ratio[i])**(1/5)/mass_ratio[i]**(3/5)
            # m2 = mass_ratio[i]*m1
            samples_res.append([mass_ratio_res[_], chirp_mass_res[_], a_1_res[_], a_2_res[_], spin_1x_res[_], spin_1y_res[_], spin_1z_res[_], spin_2x_res[_], spin_2y_res[_], spin_2z_res[_], iota_res[_], psi_res[_], ra_res[_], dec_res[_], luminosity_distance_res[_], geocent_time[_]]) 

        samples_res = np.array(samples_res)
        figure = corner.corner(samples_res, labels=[r"$q$", r"$\mathcal{M}$", r"$a_1$", r"$a_2$", r"$a_{1x}$", r"$a_{1y}$", r"$a_{1z}$", r"$a_{2x}$", r"$a_{2y}$", r"$a_{2z}$", r"$\iota$", r"POL", r"RA", r"DEC", r"D_L", r"t_c"],
                            show_titles=True, label_kwargs={"fontsize": 15},smooth=True, smooth1d=True, title_fmt='.2f', quantiles=[0.0223, 0.977],
                            )

        corner.overplot_lines(figure,injection_parameters, color='r')
        corner.overplot_points(figure, injection_parameters[None], marker="s", color='r')

        corner.overplot_lines(figure,max_likelihood_parameters, color='b')
        corner.overplot_points(figure, max_likelihood_parameters[None], marker="s", color='b')


        red_line = mlines.Line2D([], [], color='r', label='Injected')
        blue_line = mlines.Line2D([], [], color='b', label='Max likelihood')
        yellow_line = mlines.Line2D([], [], color='y', label='SNR = ' + str(SNR_tmp))
        plt.legend(handles=[red_line, blue_line, yellow_line],  bbox_to_anchor=(0.5, 8), loc = 'center', prop={'size': 20})
        plt.tight_layout()
        figure.savefig('IntermediaPlot_corner_unlens_precessing/corner_' + str(i) + '.png', dpi=450)
        plt.close()

    except FileNotFoundError:
        pass



