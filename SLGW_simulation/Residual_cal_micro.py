#画一下unlens的重构的和输入的波形。
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import matplotlib.pyplot as plt
import json
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
from astropy.cosmology import Planck15
import scipy
from scipy.interpolate import interp1d
from tqdm import tqdm
from  multiprocessing import Process,Pool
import time
import pycbc.psd
from gwpy.timeseries import TimeSeries


sample_RATE = '4KHZ'
signal_duration = 4096
G = 6.67 * 10**(-11)
c = 3 * 10**(8)
M_sun = 2 * 10**(30)
apx = 'IMRPhenomPv2'
f_lower = 10
duration = 4096
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 0.1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)


SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')

two_image = np.loadtxt('./Lensed_SampleResult/two_image_index.csv', delimiter=',')

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


apx = 'IMRPhenomPv2'

def Cal_residual(data, signal):
    return np.sum((data - signal)**2) / len(data)

def Cal_residual_freq(data, signal):
    return np.abs(np.sum((data - signal)*np.conjugate(data - signal)))

def Cal_residual_cross_freq(data1, signal1, data2, signal2):
    return np.abs(np.sum((data1 - signal1)*np.conjugate(data2 - signal2)))
# Residual_res = []
# Residual_res_max = []
# Residual_res_noise = []
Residual_res_freq = []
Residual_res_noise_freq = []
Residual_res_cross_freq = []
Residual_res_noise_cross_freq = []
Residual_SNR = []

IndexInAccept = np.loadtxt('../Paper5_cWB_Modify/Total_Data_start_time_and_snr_and_match_and_outdir_name/IndexInAccept.csv', delimiter=',')
SNR_selected = np.loadtxt('../Paper5_cWB_Modify/Total_Data_start_time_and_snr_and_match_and_outdir_name/SNR_selected.csv', delimiter=',')

for i_ in tqdm(IndexInAccept):
    i_ = int(i_)
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = int(imagenum[i_])
    save_index = 0
    IndexSumimage = int(np.sum(imagenum[0:i_]))
    

    for j_ in range(imag_num):
        
        if SNR_network[IndexSumimage] > 12:
            
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
                outdir = outdir = 'outdir_Lensed_' + Lens_Type + '_with_micro/outdir' + str(i) + "_" + str(save_index)
                file_post = outdir + '/H1_L1_V1_result.json'

                with open(file_post,"r") as f:
                    data_post = json.load(f)

                
                match_micro = []
                #最大似然点
                j = np.where(data_post['posterior']['content']['log_likelihood'] == np.max(data_post['posterior']['content']['log_likelihood']))[0][0]
                q = data_post['posterior']['content']['mass_ratio'][j]
                Mz_o = data_post['posterior']['content']['chirp_mass'][j]
                luminosity_dis = data_post['posterior']['content']['luminosity_distance'][j]
                dec = data_post['posterior']['content']['dec'][j]
                ra = data_post['posterior']['content']['ra'][j]
                theta_jn = data_post['posterior']['content']['theta_jn'][j]
                psi = data_post['posterior']['content']['psi'][j]
                phase = data_post['posterior']['content']['phase'][j]
                a_1 = data_post['posterior']['content']['a_1'][j]
                a_2 = data_post['posterior']['content']['a_2'][j]
                spin_1x = data_post['posterior']['content']['spin_1x'][j]
                spin_1y = data_post['posterior']['content']['spin_1y'][j]
                spin_1z = data_post['posterior']['content']['spin_1z'][j]
                spin_2x = data_post['posterior']['content']['spin_2x'][j]
                spin_2y = data_post['posterior']['content']['spin_2y'][j]
                spin_2z = data_post['posterior']['content']['spin_2z'][j]
                tilt_1 = data_post['posterior']['content']['tilt_1'][j]
                tilt_2 = data_post['posterior']['content']['tilt_2'][j]
                phi_12 = data_post['posterior']['content']['phi_12'][j]
                phi_jl = data_post['posterior']['content']['phi_jl'][j]
                geocent_time = data_post['posterior']['content']['geocent_time'][j]
                m1_o = Mz_o * (1 + q)**(1/5) / q ** (3/5)
                m2_o = q * m1_o
                SNR_tmp = np.sqrt(data_post['posterior']['content']['H1_matched_filter_snr'][j]['real']**2 + data_post['posterior']['content']['L1_matched_filter_snr'][j]['real']**2 + data_post['posterior']['content']['V1_matched_filter_snr'][j]['real']**2)

                hp, hc = get_td_waveform(approximant=apx,
                                        mass1=m1_o, mass2=m2_o,
                                        spin1x = spin_1x, spin1y = spin_1y,
                                        spin1z = spin_1z, spin2x = spin_2x,
                                        spin2y = spin_2y, spin2z = spin_2z,
                                        inclination = theta_jn,
                                        distance = luminosity_dis,
                                        delta_t=1.0/sample_rate, f_lower=25, f_ref = 50, coa_phase = phase)


                det_h1 = Detector('H1')
                det_l1 = Detector('L1')
                
                
                
                hp.start_time += geocent_time
                hc.start_time += geocent_time
                hp_h1, hc_h1 = hp.copy(), hc.copy()
                hp_l1, hc_l1 = hp.copy(), hc.copy()
                # 从地心到探测器的时间
                dt_h1_re = det_h1.time_delay_from_earth_center(ra, dec, hp.start_time.gpsSeconds + hp.start_time.gpsNanoSeconds / 10 ** 9)
                dt_l1_re = det_l1.time_delay_from_earth_center(ra, dec, hp.start_time.gpsSeconds + hp.start_time.gpsNanoSeconds / 10 ** 9)
                            
                hp_h1.start_time += dt_h1_re
                hc_h1.start_time += dt_h1_re

                hp_l1.start_time += dt_l1_re
                hc_l1.start_time += dt_l1_re


                # We get back the fp and fc antenna pattern weights.
                fp_h1_re, fc_h1_re = det_h1.antenna_pattern(ra, dec, psi, geocent_time)
                fp_l1_re, fc_l1_re = det_l1.antenna_pattern(ra, dec, psi, geocent_time)
                
                signal_h1 = fp_h1_re * hp_h1 + fc_h1_re * hc_h1
                signal_l1 = fp_l1_re * hp_l1 + fc_l1_re * hc_l1


                #得到频域的信号

                
                tilde_signal_h1 = scipy.fft.fft(signal_h1.data) * signal_h1.delta_t
                tilde_signal_l1 = scipy.fft.fft(signal_l1.data) * signal_l1.delta_t

                freq_h1 = scipy.fft.fftfreq(len(signal_h1), signal_h1.delta_t)
                freq_l1 = scipy.fft.fftfreq(len(signal_l1), signal_l1.delta_t)


                Bilby_func_h1 = interp1d(freq_h1, tilde_signal_h1)
                Bilby_func_l1 = interp1d(freq_l1, tilde_signal_l1)


                
                
                #输入的
                luminosity_dis = Planck15.luminosity_distance(SampleParam[0][i]).value
                hp_orign, hc_orign = get_td_waveform(approximant=apx,
                                                mass1=SampleParam[1][i] * (1 + SampleParam[0][i]), mass2=SampleParam[2][i] * (1 + SampleParam[0][i]),
                                                spin1z = SampleParam[3][i], spin2z = SampleParam[4][i],
                                                inclination = SampleParam[5][i], 
                                                distance = luminosity_dis,
                                                delta_t=1.0/sample_rate, f_lower=f_lower, f_ref = 50, coa_phase = 0)
                            
                            
                hp_20, hc_20 = get_td_waveform(approximant=apx,
                                                mass1=SampleParam[1][i] * (1 + SampleParam[0][i]), mass2=SampleParam[2][i] * (1 + SampleParam[0][i]),
                                                spin1z = SampleParam[3][i], spin2z = SampleParam[4][i],
                                                inclination = SampleParam[5][i], 
                                                distance = luminosity_dis,
                                                delta_t=1.0/sample_rate, f_lower=20, f_ref = 50, coa_phase = 0)
                            
                    
                #把时间调成0试试，这样做意味着警报时间和开始时间对齐了

                # merger_index = np.where(hp.sample_times <= 0)[0][-1]
                hp_orign.start_time = 0
                hc_orign.start_time = 0
                
                

                det_h1 = Detector('H1')
                det_l1 = Detector('L1')
                


                #不同探测器的波形
                hp_h1, hc_h1 = hp_orign.copy(), hc_orign.copy()
                hp_l1, hc_l1 = hp_orign.copy(), hc_orign.copy()
            
                
                #得到不同探测器的时间
                start_time_signal = int(SampleParam9_td) - (hp_orign.duration - hp_20.duration) #保证30Hz的开始时间是在start time上
                start_time = int(SampleParam9_td)
                
                declination = SampleParam[8][i]
                right_ascension = SampleParam[7][i]
                polarization = SampleParam[6][i]
                
                #校准不同探测的波形时间
                hp_h1.start_time += start_time_signal
                hc_h1.start_time += start_time_signal
                hp_l1.start_time += start_time_signal
                hc_l1.start_time += start_time_signal
                
                
                
                # We get back the fp and fc antenna pattern weights.
                fp_h1, fc_h1 = det_h1.antenna_pattern(right_ascension, declination, polarization, start_time)
                fp_l1, fc_l1 = det_l1.antenna_pattern(right_ascension, declination, polarization, start_time)
                
                # 从地心到探测器的时间
                dt_h1 = det_h1.time_delay_from_earth_center(right_ascension, declination, start_time)
                dt_l1 = det_l1.time_delay_from_earth_center(right_ascension, declination, start_time)
                
            
                
                #校准不同探测的波形时间
                hp_h1.start_time += dt_h1 
                hc_h1.start_time += dt_h1
                hp_l1.start_time += dt_l1 
                hc_l1.start_time += dt_l1
                

                ## Apply the factors to get the detector frame strain
                signal_h1_inj = fp_h1 * hp_h1 + fc_h1 * hc_h1
                signal_l1_inj = fp_l1 * hp_l1 + fc_l1 * hc_l1
                
                
                tilde_signal_h1_inj = scipy.fft.fft(signal_h1_inj.data) * signal_h1_inj.delta_t

                freq_h1_inj = scipy.fft.fftfreq(len(signal_h1_inj), signal_h1_inj.delta_t)

                tilde_signal_l1_inj = scipy.fft.fft(signal_l1_inj.data) * signal_l1_inj.delta_t

                freq_l1_inj = scipy.fft.fftfreq(len(signal_l1_inj), signal_l1_inj.delta_t)
                
                """透镜的放大"""
                
                
                Ff_inter = np.loadtxt("../Paper4_MicroLens_Modify/Ffabs" + Lens_Type + "_0_300/Total" + str(i) + "_" + str(save_index) + ".csv", delimiter=',')
                ThetaF_inter = np.loadtxt("../Paper4_MicroLens_Modify/ThetaF" + Lens_Type + "_0_300/Total" + str(i) + "_" + str(save_index) + ".csv", delimiter=',')
                    


                tilde_signal_h1_inj_lensed = tilde_signal_h1_inj * Ff_inter * np.exp(complex(0,1) * ThetaF_inter)
                tilde_signal_l1_inj_lensed = tilde_signal_l1_inj * Ff_inter * np.exp(complex(0,1) * ThetaF_inter)
                
                tilde_signal_h1_inj_lensed_func = interp1d(freq_h1_inj, tilde_signal_h1_inj_lensed)
                tilde_signal_l1_inj_lensed_func = interp1d(freq_l1_inj, tilde_signal_l1_inj_lensed)
                 
                signal_h1_inj_lensed = np.real(scipy.fft.ifft(tilde_signal_h1_inj_lensed) / signal_h1_inj.delta_t)
                signal_l1_inj_lensed = np.real(scipy.fft.ifft(tilde_signal_l1_inj_lensed) / signal_l1_inj.delta_t)
                
                index_plot = np.where((freq_h1 > f_lower)&(freq_h1 < 1024))
                
                plt.title('SNR = ' + str(SNR_tmp) + ' Type = ' + Lens_Type)
                plt.loglog(freq_h1[index_plot], np.abs(tilde_signal_h1[index_plot]), label = 'Bilby')
                plt.loglog(freq_h1_inj[index_plot], np.abs(tilde_signal_h1_inj_lensed[index_plot]), '--', label='inj micro')
                plt.xlim(1, 500)
                plt.ylim(10**(-26), 10**(-23))
                plt.legend()
                plt.savefig('IntermediaPlot_waveform_micro/H1_Micro_Tilde_res_vs_inj_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()

                plt.title('SNR = ' + str(SNR_tmp) + ' Type = ' + Lens_Type)
                plt.loglog(freq_l1[index_plot], np.abs(tilde_signal_l1[index_plot]), label = 'Bilby')
                plt.loglog(freq_l1_inj[index_plot], np.abs(tilde_signal_l1_inj_lensed[index_plot]), '--', label='inj micro')
                plt.xlim(1, 500)
                plt.ylim(10**(-26), 10**(-23))
                plt.legend()
                plt.savefig('IntermediaPlot_waveform_micro/L1_Micro_Tilde_res_vs_inj_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()

                plt.title('SNR = ' + str(SNR_tmp) + ' Type = ' + Lens_Type)
                plt.plot(signal_h1.sample_times, signal_h1, label='Bilby')
                plt.plot(signal_h1_inj.sample_times, signal_h1_inj_lensed, '--', label='inj micro')
                plt.xlim(geocent_time - 1, geocent_time + 0.1)
                plt.legend()
                plt.savefig('IntermediaPlot_waveform_micro/H1_Time_res_vs_inj_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()

                plt.title('SNR = ' + str(SNR_tmp) + ' Type = ' + Lens_Type)
                plt.plot(signal_l1.sample_times, signal_l1, label='Bilby')
                plt.plot(signal_l1_inj.sample_times, signal_l1_inj_lensed, '--', label='inj micro')
                plt.xlim(geocent_time - 1, geocent_time + 0.1)
                plt.legend()
                plt.savefig('IntermediaPlot_waveform_micro/L1_Time_res_vs_inj_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()
                
                
                """下面开始计算residual"""
                """读取数据"""
                fname_list = ['./Total_Sim_GW_Data_Lensed_' + Lens_Type + '_with_micro/H1/H-H1_GWOSC_4KHZ_R1-'+str(int(SampleParam9_td - signal_duration/2))+'-4096_' + str(i) + "_" + str(save_index) + '.gwf', './Total_Sim_GW_Data_Lensed_' + Lens_Type + '_with_micro/L1/L-L1_GWOSC_4KHZ_R1-'+str(int(SampleParam9_td - signal_duration/2))+'-4096_' + str(i) + "_" + str(save_index) + '.gwf']
                fchannel_list = ['H1:GWOSC_4KHZ_R1_STRAIN', 'L1:GWOSC_4KHZ_R1_STRAIN']
                data_h1 = TimeSeries.read(fname_list[0], channel = fchannel_list[0])
                data_l1 = TimeSeries.read(fname_list[1], channel = fchannel_list[1])
                
                # data_h1 = data_h1.highpass(20)
                # data_l1 = data_l1.highpass(20)
                # data_h1 = data_h1.lowpass(300)
                # data_l1 = data_l1.lowpass(300)
                data_start_index_h1 = np.where(np.array(data_h1.times) >= signal_h1.sample_times.data[0])[0][0]
                data_end_index_h1 = np.where(np.array(data_h1.times) <= signal_h1.sample_times.data[-1])[0][-1]

                data_start_index_l1 = np.where(np.array(data_l1.times) >= signal_l1.sample_times.data[0])[0][0]
                data_end_index_l1 = np.where(np.array(data_l1.times) <= signal_l1.sample_times.data[-1])[0][-1]
                
                signal_rec_inter_h1 = interp1d(signal_h1.sample_times, signal_h1)
                signal_rec_inter_l1 = interp1d(signal_l1.sample_times, signal_l1)
                # signal_rec_inter = interp1d(signal_h1_inj.sample_times, signal_h1_inj_lensed)
                data_segment_h1 = data_h1[data_start_index_h1:data_end_index_h1]
                data_segment_l1 = data_l1[data_start_index_l1:data_end_index_l1]

                '''得到数据中的频域信号'''
                tilde_data_segment_h1 = scipy.fft.fft(data_segment_h1.value) * data_segment_h1.dt

                freq_data_segment_h1 = scipy.fft.fftfreq(len(data_segment_h1.value), data_segment_h1.dt)


                data_segment_func_h1 = interp1d(freq_data_segment_h1, tilde_data_segment_h1)

                tilde_data_segment_l1 = scipy.fft.fft(data_segment_l1.value) * data_segment_l1.dt

                freq_data_segment_l1 = scipy.fft.fftfreq(len(data_segment_l1.value), data_segment_l1.dt)


                data_segment_func_l1 = interp1d(freq_data_segment_l1, tilde_data_segment_l1)
                
                
                
                
                
                '''得到纯噪声的频域信号'''
                tilde_noise_h1 = scipy.fft.fft(data_h1[data_start_index_h1+40960:data_end_index_h1+40960].value) * data_h1[data_start_index_h1+40960:data_end_index_h1+40960].dt

                freq_noise_h1 = scipy.fft.fftfreq(len(data_h1[data_start_index_h1+40960:data_end_index_h1+40960].value), data_h1[data_start_index_h1+40960:data_end_index_h1+40960].dt)


                noise_func_h1 = interp1d(freq_noise_h1, tilde_noise_h1)

                tilde_noise_l1 = scipy.fft.fft(data_l1[data_start_index_l1+40960:data_end_index_l1+40960].value) * data_l1[data_start_index_l1+40960:data_end_index_l1+40960].dt

                freq_noise_l1 = scipy.fft.fftfreq(len(data_l1[data_start_index_l1+40960:data_end_index_l1+40960].value), data_l1[data_start_index_l1+40960:data_end_index_l1+40960].dt)


                noise_func_l1 = interp1d(freq_noise_l1, tilde_noise_l1)
                
                 
                plt.title('SNR = ' + str(SNR_tmp) + ' Type = ' + Lens_Type)
                plt.plot(data_segment_h1.times, data_segment_h1.value)
                # plt.plot(data[data_start_index:data_end_index].times, data[data_start_index-40960:data_end_index-40960].value)
                plt.plot(signal_h1.sample_times, signal_h1, label='Bilby')
                plt.plot(data_segment_h1.times, signal_rec_inter_h1(data_segment_h1.times), '--' ,label='interp1d')
                plt.xlim(signal_h1.sample_times[0] + signal_h1.duration * 5 / 10, signal_h1.sample_times[0] + signal_h1.duration * 8 / 10)
                plt.savefig('IntermediaPlot_waveform_micro/H1_Read_Time_res_vs_data_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()

                plt.title('SNR = ' + str(SNR_tmp) + ' Type = ' + Lens_Type)
                plt.plot(data_segment_l1.times, data_segment_l1.value)
                # plt.plot(data[data_start_index:data_end_index].times, data[data_start_index-40960:data_end_index-40960].value)
                plt.plot(signal_l1.sample_times, signal_l1, label='Bilby')
                plt.plot(data_segment_l1.times, signal_rec_inter_l1(data_segment_l1.times), '--' ,label='interp1d')
                plt.xlim(signal_l1.sample_times[0] + signal_l1.duration * 5 / 10, signal_l1.sample_times[0] + signal_l1.duration * 8 / 10)
                plt.savefig('IntermediaPlot_waveform_micro/L1_Read_Time_res_vs_data_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()
                
                plt.plot(data_segment_h1.times, data_segment_h1.value - signal_rec_inter_h1(data_segment_h1.times))
                # plt.plot(data_segment.times, data_segment.value)
                plt.plot(data_h1[data_start_index_h1:data_end_index_h1].times, data_h1[data_start_index_h1-40960:data_end_index_h1-40960].value)
                # plt.plot(signal_h1.sample_times, signal_h1, label='Bilby')
                # plt.plot(data_segment.times, signal_rec_inter(data_segment.times), '--' ,label='interp1d')
                # plt.xlim(signal_h1.sample_times[0] + signal_h1.duration * 6 / 10, signal_h1.sample_times[0] + signal_h1.duration * 10 / 10)
                plt.savefig('IntermediaPlot_waveform_micro/H1_Residual_' + str(i) + '.png', dpi=450)
                plt.close()

                
                plt.plot(data_segment_l1.times, data_segment_l1.value - signal_rec_inter_l1(data_segment_l1.times))
                # plt.plot(data_segment.times, data_segment.value)
                plt.plot(data_l1[data_start_index_l1:data_end_index_l1].times, data_l1[data_start_index_l1-40960:data_end_index_l1-40960].value)
                # plt.plot(signal_h1.sample_times, signal_h1, label='Bilby')
                # plt.plot(data_segment.times, signal_rec_inter(data_segment.times), '--' ,label='interp1d')
                # plt.xlim(signal_h1.sample_times[0] + signal_h1.duration * 6 / 10, signal_h1.sample_times[0] + signal_h1.duration * 10 / 10)
                plt.savefig('IntermediaPlot_waveform_micro/L1_Residual_' + str(i) + '.png', dpi=450)
                plt.close()
                
                
                freq_used_in_freq_residual_cal = np.arange(25, 250)
                plt.semilogy(freq_used_in_freq_residual_cal, np.abs(data_segment_func_h1(freq_used_in_freq_residual_cal)), label='data')
                
                plt.semilogy(freq_used_in_freq_residual_cal, np.abs(tilde_signal_h1_inj_lensed_func(freq_used_in_freq_residual_cal)), label='inj')

                plt.semilogy(freq_used_in_freq_residual_cal, np.abs(Bilby_func_h1(freq_used_in_freq_residual_cal)), label='bilby')
                plt.legend()
                # plt.ylim(10**(-25), 10**(-20))
                plt.ylim(10**(-26), 10**(-23))
                plt.savefig('IntermediaPlot_waveform_micro/H1_freq_domain_abs_' + str(i) + '.png', dpi=450)
                plt.close()

                plt.semilogy(freq_used_in_freq_residual_cal, np.abs(data_segment_func_l1(freq_used_in_freq_residual_cal)), label='data')
                
                plt.semilogy(freq_used_in_freq_residual_cal, np.abs(tilde_signal_l1_inj_lensed_func(freq_used_in_freq_residual_cal)), label='inj')

                plt.semilogy(freq_used_in_freq_residual_cal, np.abs(Bilby_func_l1(freq_used_in_freq_residual_cal)), label='bilby')
                plt.legend()
                # plt.ylim(10**(-25), 10**(-20))
                plt.ylim(10**(-26), 10**(-23))
                plt.savefig('IntermediaPlot_waveform_micro/L1_freq_domain_abs_' + str(i) + '.png', dpi=450)
                plt.close()

                
                
                # Residual_res.append(Cal_residual(data_segment.value, signal_rec_inter(data_segment.times)))
                # Residual_res_noise.append(Cal_residual(data[data_start_index+40960:data_end_index+40960].value, 0))
                # Residual_res_max.append(Cal_residual(data_segment.value, 0))
                
                Residual_res_freq.append(Cal_residual_freq(data_segment_func_h1(freq_used_in_freq_residual_cal), Bilby_func_h1(freq_used_in_freq_residual_cal)) + Cal_residual_freq(data_segment_func_h1(-freq_used_in_freq_residual_cal), Bilby_func_h1(-freq_used_in_freq_residual_cal)))
                Residual_res_cross_freq.append(Cal_residual_cross_freq(data_segment_func_h1(freq_used_in_freq_residual_cal), Bilby_func_h1(freq_used_in_freq_residual_cal), data_segment_func_l1(freq_used_in_freq_residual_cal), Bilby_func_l1(freq_used_in_freq_residual_cal)) + Cal_residual_cross_freq(data_segment_func_h1(-freq_used_in_freq_residual_cal), Bilby_func_h1(-freq_used_in_freq_residual_cal), data_segment_func_l1(-freq_used_in_freq_residual_cal), Bilby_func_l1(-freq_used_in_freq_residual_cal)))
                Residual_res_noise_freq.append(Cal_residual_freq(noise_func_h1(freq_used_in_freq_residual_cal), 0) + Cal_residual_freq(noise_func_h1(-freq_used_in_freq_residual_cal), 0))
                Residual_res_noise_cross_freq.append(Cal_residual_cross_freq(noise_func_h1(freq_used_in_freq_residual_cal), 0, noise_func_l1(freq_used_in_freq_residual_cal), 0) + Cal_residual_cross_freq(noise_func_h1(-freq_used_in_freq_residual_cal), 0, noise_func_l1(-freq_used_in_freq_residual_cal), 0))
                Residual_SNR.append(SNR_tmp)
                 
            except FileNotFoundError:
                pass
            save_index += 1
        IndexSumimage += 1


# np.savetxt('IntermediaPlot_waveform_micro/Residual_res.csv', Residual_res, delimiter=',')
# np.savetxt('IntermediaPlot_waveform_micro/Residual_res_max.csv', Residual_res_max, delimiter=',')
# np.savetxt('IntermediaPlot_waveform_micro/Residual_res_noise.csv', Residual_res_noise, delimiter=',')

np.savetxt('IntermediaPlot_waveform_micro/Residual_res_freq.csv', Residual_res_freq, delimiter=',')
np.savetxt('IntermediaPlot_waveform_micro/Residual_res_noise_freq.csv', Residual_res_noise_freq, delimiter=',')
np.savetxt('IntermediaPlot_waveform_micro/Residual_res_cross_freq.csv', Residual_res_cross_freq, delimiter=',')
np.savetxt('IntermediaPlot_waveform_micro/Residual_res_noise_cross_freq.csv', Residual_res_noise_cross_freq, delimiter=',')
np.savetxt('IntermediaPlot_waveform_micro/Residual_SNR.csv', Residual_SNR, delimiter=',')


# Residual_res_macro = np.loadtxt('IntermediaPlot_waveform_macro/Residual_res.csv', delimiter=',')
# Residual_res = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res.csv', delimiter=',')
# Residual_res_max = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res_max.csv', delimiter=',')
# Residual_res_noise = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res_noise.csv', delimiter=',')
Residual_SNR = np.loadtxt('IntermediaPlot_waveform_micro/Residual_SNR.csv', delimiter=',')
Residual_SNR_macro = np.loadtxt('IntermediaPlot_waveform_macro/Residual_SNR.csv', delimiter=',')

Residual_res_macro_freq = np.loadtxt('IntermediaPlot_waveform_macro/Residual_res_freq.csv', delimiter=',')
Residual_res_macro_cross_freq = np.loadtxt('IntermediaPlot_waveform_macro/Residual_res_cross_freq.csv', delimiter=',')
Residual_res_freq = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res_freq.csv', delimiter=',')
Residual_res_noise_freq = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res_noise_freq.csv', delimiter=',')
Residual_res_cross_freq = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res_cross_freq.csv', delimiter=',')
Residual_res_noise_cross_freq = np.loadtxt('IntermediaPlot_waveform_micro/Residual_res_noise_cross_freq.csv', delimiter=',')


# plt.hist(np.log10(np.abs(Residual_res)), label='micro', histtype='step', density=True, bins=10)
# plt.hist(np.log10(np.abs(Residual_res_macro)), label='macro', histtype='step', density=True, bins=10)
# plt.hist(np.log10(np.abs(Residual_res_noise)), label='noise', histtype='step', density=True, bins=10)
# plt.legend()
# plt.savefig('IntermediaPlot_waveform_micro/statistic_Residual.png', dpi=450)
# plt.close()

# plt.hist(np.log10(np.abs(Residual_res_max)), label='micro max', histtype='step', density=True, bins=10)
# plt.hist(np.log10(np.abs(Residual_res_macro)), label='macro', histtype='step', density=True, bins=10)
# plt.hist(np.log10(np.abs(Residual_res_noise)), label='noise', histtype='step', density=True, bins=10)
# plt.legend()
# plt.savefig('IntermediaPlot_waveform_micro/max_statistic_Residual.png', dpi=450)
# plt.close()



plt.hist(np.log10(np.abs(Residual_res_freq)), label='micro', histtype='step', density=True)
plt.hist(np.log10(np.abs(Residual_res_macro_freq)), label='macro', histtype='step', density=True)
plt.hist(np.log10(np.abs(Residual_res_noise_freq)), label='noise', histtype='step', density=True)
plt.legend()




plt.hist(np.log10(np.abs(Residual_res_cross_freq)), label='micro', histtype='step', density=True)
plt.hist(np.log10(np.abs(Residual_res_macro_cross_freq)), label='macro', histtype='step', density=True)
plt.hist(np.log10(np.abs(Residual_res_noise_cross_freq)), label='noise', histtype='step', density=True)
plt.legend()