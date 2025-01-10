'''
修改成10Hz了，与原先的版本不同，我这里保存的文件名字是i而不是i_，与unlens和micro total一致了。
'''

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import matplotlib.pyplot as pp
import pycbc.noise
import pycbc.psd
import pylab
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
from astropy.cosmology import Planck15
from tqdm import tqdm
from pycbc import frame
import pycbc.filter
from pycbc.filter import highpass_fir, matched_filter
from pycbc.psd import welch, interpolate
from  multiprocessing import Process,Pool
from scipy.interpolate import interp1d
import scipy
import sys
import matplotlib.pyplot as plt
import gc
# Generate a PSD using an analytic expression for 
# the full design Advanced LIGO noise curve
f_lower = 10
duration = 4096
sample_rate = 4096
sample_RATE = '4KHZ'
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 0.1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)


SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

apx = 'IMRPhenomPv2'

# Now, let's generate noise that has the same spectrum
htilde_l1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=1000)
htilde_h1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=2000)
htilde_v1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=3000)


AcceptLensIndex = np.loadtxt('./SampleResult/AcceptLensIndex.csv', delimiter=',')

magnification = np.loadtxt('./SampleResult/magnification.csv', delimiter=',')
timedelay = np.loadtxt('./SampleResult/timedelayDays.csv', delimiter=',')
kappa = np.loadtxt('./SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('./SampleResult/gamma.csv', delimiter=',')
imagenum = np.loadtxt('./SampleResult/imagenumber.csv', delimiter=',')









def SNR_cal(IndexStart, IndexEnd):
    Error = []
    H1_snr_set = []
    L1_snr_set = []
    V1_snr_set = []
  
    trigger_time_H1 = []
    trigger_time_L1 = []
    trigger_time_V1 = []
  
    gw_duration = [] #10Hz的持续时间
    gw_duration_20 = [] #20Hz的持续时间
    
    dt_h1_from_earth_center = []
    dt_l1_from_earth_center = []
    dt_v1_from_earth_center = []
    
    for i_ in tqdm(range(IndexStart, IndexEnd)):
        
        # print(i_)
        i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
        imag_num = imagenum[i_]
        
        sum_num = int(np.sum(imagenum[0:i_]))
        save_index = 0
        for j_ in range(int(imag_num)):
            
            td_second = timedelay[sum_num] * 24 * 3600
            SampleParam9_td = SampleParam[9][i] + td_second
            
            sqrt_mag_abs = np.sqrt(np.abs(magnification[sum_num]))
            if 1 - kappa[sum_num] - gamma[sum_num] > 0 and 1 - kappa[sum_num] + gamma[sum_num] > 0:    
                """================================="""
                Lens_Type = "Minimum"
            elif (1 - kappa[sum_num] - gamma[sum_num]) * (1 - kappa[sum_num] + gamma[sum_num]) < 0:
                Lens_Type = "Saddle"
            elif 1 - kappa[sum_num] - gamma[sum_num] < 0 and 1 - kappa[sum_num] + gamma[sum_num] < 0:
                Lens_Type = "Maximum"
                    
            sum_num += 1
            """透镜的放大"""
            if Lens_Type == "Minimum": 
                Ff_abs = [sqrt_mag_abs] * len(lens_freq)
                ThetaF = [0]*len(lens_freq) 
                
            elif Lens_Type == "Saddle":
                Ff_abs = [sqrt_mag_abs] * len(lens_freq)
                ThetaF1 = [np.pi/2]*(len(lens_freq)//2)
                ThetaF2 = [-np.pi/2]*(len(lens_freq)//2) 
                ThetaF = ThetaF1 + ThetaF2
                
            elif Lens_Type == "Maximum":
                Ff_abs = [sqrt_mag_abs] * len(lens_freq)
                ThetaF = [np.pi]*len(lens_freq)


            Ff_inter = interp1d(lens_freq, Ff_abs)
            ThetaF_inter = interp1d(lens_freq, ThetaF)
            """================================"""



            # Equivelantly in the time domain
            hoft_h1 = htilde_h1.to_timeseries()
            hoft_l1 = htilde_l1.to_timeseries()
            hoft_v1 = htilde_v1.to_timeseries()
            hoft_h1.start_time = 0
            hoft_l1.start_time = 0
            hoft_v1.start_time = 0
         
            """对时域的噪音hoft用welch方法重新估计psd，并与原始的psd对比"""

            # Use Welch's method with 4s segments
            psd_h1_estimated = interpolate(welch(hoft_h1), 1.0 / hoft_h1.duration)
            psd_l1_estimated = interpolate(welch(hoft_l1), 1.0 / hoft_l1.duration)
            psd_v1_estimated = interpolate(welch(hoft_v1), 1.0 / hoft_v1.duration)
        
            # pp.plot(hoft_l1_h1.sample_times, hoft_l1_h1)
            # print('start time = ', hoft_l1_h1.start_time)
            luminosity_dis = Planck15.luminosity_distance(SampleParam[0][i]).value
                    # NOTE: Inclination runs from 0 to pi, with poles at 0 and pi
                    #       coa_phase runs from 0 to 2 pi.
            
            
            hp_orign, hc_orign = get_td_waveform(approximant=apx,
                                    mass1=SampleParam[1][i] * (1 + SampleParam[0][i]), mass2=SampleParam[2][i] * (1 + SampleParam[0][i]),
                                    spin1z = SampleParam[3][i], spin2z = SampleParam[4][i],
                                    inclination = SampleParam[5][i],
                                    distance = luminosity_dis,
                                    delta_t=1.0/sample_rate, f_lower=f_lower, f_ref = 50, coa_phase=0)
            
            hp_20, hc_20 = get_td_waveform(approximant=apx,
                                        mass1=SampleParam[1][i] * (1 + SampleParam[0][i]), mass2=SampleParam[2][i] * (1 + SampleParam[0][i]),
                                        spin1z = SampleParam[3][i], spin2z = SampleParam[4][i],
                                        inclination = SampleParam[5][i],
                                        distance = luminosity_dis,
                                        delta_t=1.0/sample_rate, f_lower=20, f_ref = 50, coa_phase=0)
                    
            
            
            merger_index = np.where(hp_orign.sample_times <= 0)[0][-1]
            if(merger_index == np.where(hc_orign.sample_times <= 0)[0][-1]):
                pass
            else:
                print("Error !!!!!!!!!!!!!!, GW index i = ", i)
            
            #把时间调成0试试，这样做意味着警报时间和开始时间对齐了
            hp_orign.start_time = 0
            hc_orign.start_time = 0
            
            hp_20.start_time = 0
            hc_20.start_time = 0
            
            gw_duration.append(hp_orign.duration)
            gw_duration_20.append(hp_20.duration)
            
            
            det_h1 = Detector('H1')
            det_l1 = Detector('L1')
            det_v1 = Detector('V1')
            
            
            #不同探测器的波形
            hp_h1, hc_h1 = hp_orign.copy(), hc_orign.copy()
            hp_l1, hc_l1 = hp_orign.copy(), hc_orign.copy()
            hp_v1, hc_v1 = hp_orign.copy(), hc_orign.copy()

            
            #得到不同探测器的时间
            start_time_signal = int(SampleParam9_td) - (hp_orign.duration - hp_20.duration) #保证20Hz的开始时间是在start time上
            start_time = int(SampleParam9_td)
            

         
            declination = SampleParam[8][i]
            right_ascension = SampleParam[7][i]
            polarization = SampleParam[6][i]
            
            #校准不同探测的波形时间
            hp_h1.start_time += start_time_signal
            hc_h1.start_time += start_time_signal
            hp_l1.start_time += start_time_signal
            hc_l1.start_time += start_time_signal
            hp_v1.start_time += start_time_signal
            hc_v1.start_time += start_time_signal
        
            #校准不同探测器噪声时间
            hoft_h1.start_time += start_time
            hoft_l1.start_time += start_time
            hoft_v1.start_time += start_time
        
            hoft_h1.start_time -= duration // 2
            hoft_l1.start_time -= duration // 2
            hoft_v1.start_time -= duration // 2
         
            # We get back the fp and fc antenna pattern weights.
            fp_h1, fc_h1 = det_h1.antenna_pattern(right_ascension, declination, polarization, start_time)
            fp_l1, fc_l1 = det_l1.antenna_pattern(right_ascension, declination, polarization, start_time)
            fp_v1, fc_v1 = det_v1.antenna_pattern(right_ascension, declination, polarization, start_time)
            # 从地心到探测器的时间
            dt_h1 = det_h1.time_delay_from_earth_center(right_ascension, declination, start_time)
            dt_l1 = det_l1.time_delay_from_earth_center(right_ascension, declination, start_time)
            dt_v1 = det_v1.time_delay_from_earth_center(right_ascension, declination, start_time)
            
            dt_h1_from_earth_center.append(dt_h1)
            dt_l1_from_earth_center.append(dt_l1)
            dt_v1_from_earth_center.append(dt_v1)
            
            #校准不同探测的波形时间
            hp_h1.start_time += dt_h1 
            hc_h1.start_time += dt_h1
            hp_l1.start_time += dt_l1
            hc_l1.start_time += dt_l1
            hp_v1.start_time += dt_v1
            hc_v1.start_time += dt_v1

            ## Apply the factors to get the detector frame strain
            signal_h1 = fp_h1 * hp_h1 + fc_h1 * hc_h1
            signal_l1 = fp_l1 * hp_l1 + fc_l1 * hc_l1
            signal_v1 = fp_v1 * hp_v1 + fc_v1 * hc_v1
            
            ##对时域信号做傅立叶变换
            tilde_signal_h1 = scipy.fft.fft(signal_h1.data) * signal_h1.delta_t
            tilde_signal_l1 = scipy.fft.fft(signal_l1.data) * signal_l1.delta_t
            tilde_signal_v1 = scipy.fft.fft(signal_v1.data) * signal_v1.delta_t
            
            freq_h1 = scipy.fft.fftfreq(len(signal_h1), signal_h1.delta_t)
            freq_l1 = scipy.fft.fftfreq(len(signal_l1), signal_l1.delta_t)
            freq_v1 = scipy.fft.fftfreq(len(signal_v1), signal_v1.delta_t)
            
            
            #注入透镜的振幅（新版本）
            tilde_signal_h1_lensed = tilde_signal_h1 * Ff_inter(freq_h1) * np.exp(complex(0,1) * ThetaF_inter(freq_h1))
            tilde_signal_l1_lensed = tilde_signal_l1 * Ff_inter(freq_l1) * np.exp(complex(0,1) * ThetaF_inter(freq_l1))
            tilde_signal_v1_lensed = tilde_signal_v1 * Ff_inter(freq_v1) * np.exp(complex(0,1) * ThetaF_inter(freq_v1))
           
            #画一下加入透镜之前和之后在频域上的不同
            # pp.loglog(freq_h1, np.abs(tilde_signal_h1))
            # pp.loglog(freq_h1, np.abs(tilde_signal_h1_lensed), label='lensed')
            # pp.legend()
            
            
            #得到透镜以后的时域
            signal_h1_lensed = np.real(scipy.fft.ifft(tilde_signal_h1_lensed) / signal_h1.delta_t)
            signal_l1_lensed = np.real(scipy.fft.ifft(tilde_signal_l1_lensed) / signal_l1.delta_t)
            signal_v1_lensed = np.real(scipy.fft.ifft(tilde_signal_v1_lensed) / signal_v1.delta_t)

            #画一下加入透镜之前和之后时域上的不同
            # pp.plot(signal_h1.sample_times, signal_h1)
            # pp.plot(signal_h1.sample_times, signal_h1_lensed, '--', label='lensed')
            # pp.legend()
            # pp.xlim(end_time_h1+hp_h1.duration*3/4, end_time_h1+hp_h1.duration)
            
            
            signal_h1_func = interp1d(signal_h1.sample_times, signal_h1_lensed)
            signal_l1_func = interp1d(signal_l1.sample_times, signal_l1_lensed)
            signal_v1_func = interp1d(signal_v1.sample_times, signal_v1_lensed)
            
            #储存merger时间，也就是hp等于0的位置处所对应的时间。
            trigger_time_H1.append(signal_h1.sample_times[merger_index])
            trigger_time_L1.append(signal_l1.sample_times[merger_index])
            trigger_time_V1.append(signal_v1.sample_times[merger_index])
            
            signal_with_noise_h1 = hoft_h1.copy()
            signal_with_noise_l1 = hoft_l1.copy()
            signal_with_noise_v1 = hoft_v1.copy()
            
            #信号周围的开始时间和结束时间
            signal_with_noise_h1_index_start = np.where(signal_with_noise_h1.sample_times.data>=signal_h1.sample_times[0])[0][0]
            signal_with_noise_l1_index_start = np.where(signal_with_noise_l1.sample_times.data>=signal_l1.sample_times[0])[0][0]
            signal_with_noise_v1_index_start = np.where(signal_with_noise_v1.sample_times.data>=signal_v1.sample_times[0])[0][0]

            signal_with_noise_h1_index_end = np.where(signal_with_noise_h1.sample_times.data<=signal_h1.sample_times[-1])[0][-1]
            signal_with_noise_l1_index_end = np.where(signal_with_noise_l1.sample_times.data<=signal_l1.sample_times[-1])[0][-1]
            signal_with_noise_v1_index_end = np.where(signal_with_noise_v1.sample_times.data<=signal_v1.sample_times[-1])[0][-1]
            
            
            signal_with_noise_h1.data[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end] += signal_h1_func(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end])
            signal_with_noise_l1.data[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end] += signal_l1_func(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end])
            signal_with_noise_v1.data[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end] += signal_v1_func(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end])
            
            
            #===============================
            """matched filter SNR， 依然用无透镜的模板来匹配"""
            # Calculate the complex Signal-to-noise. This is a complex vector
            # because the signal could have any phase.
            signal_h1.resize(len(signal_with_noise_h1))
            signal_l1.resize(len(signal_with_noise_l1))
            signal_v1.resize(len(signal_with_noise_v1))
          
            snr_h1 = matched_filter(signal_h1, signal_with_noise_h1, psd=psd_h1_estimated, low_frequency_cutoff=f_lower)
            snr_l1 = matched_filter(signal_l1, signal_with_noise_l1, psd=psd_l1_estimated, low_frequency_cutoff=f_lower)
            snr_v1 = matched_filter(signal_v1, signal_with_noise_v1, psd=psd_v1_estimated, low_frequency_cutoff=f_lower)
           

            # Remove regions corrupted by filter wraparound
            snr_h1 = snr_h1[len(snr_h1) // 4: len(snr_h1) * 3 // 4]
            snr_l1 = snr_l1[len(snr_l1) // 4: len(snr_l1) * 3 // 4]
            snr_v1 = snr_v1[len(snr_v1) // 4: len(snr_v1) * 3 // 4]
    
           
            H1_snr_set.append(np.max(abs(snr_h1.data)))
            L1_snr_set.append(np.max(abs(snr_l1.data)))
            V1_snr_set.append(np.max(abs(snr_v1.data)))
            if(np.sqrt(H1_snr_set[-1]**2 + L1_snr_set[-1]**2 + V1_snr_set[-1]**2)>=12):
                
                plt.figure(i)
                plt.plot(signal_h1.sample_times, signal_h1.data, 'r', label = 'H1')
                plt.plot(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end], signal_h1_func(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end]), 'r--')
                plt.plot(signal_l1.sample_times, signal_l1.data, 'b', label = 'L1')
                plt.plot(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end], signal_l1_func(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end]), 'b--')
                plt.plot(signal_v1.sample_times, signal_v1.data, 'k', label = 'V1')
                plt.plot(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end], signal_v1_func(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end]), 'k--')
                plt.xlim(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end][len(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end])//2], signal_with_noise_h1.sample_times[signal_with_noise_h1_index_end] + 0.05)
                plt.legend()
                plt.xlabel('t[s]')
                plt.ylabel('strain')
                plt.savefig('./IntermediaPlot_Gen_inter_Macro/waveform_' + str(i) + '_' +str(save_index) +'.png', dpi=450)
                plt.close()
                
                frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/H1/H-H1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_h1.start_time.gpsSeconds) + '-4096_' + str(i) + '_'  + str(save_index) + '.gwf','H1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_h1)
                frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/L1/L-L1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_l1.start_time.gpsSeconds) + '-4096_' + str(i) + '_'  + str(save_index) + '.gwf','L1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_l1)
                frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/V1/V-V1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_v1.start_time.gpsSeconds) + '-4096_' + str(i) + '_'  + str(save_index) + '.gwf','V1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_v1)
                
                save_index += 1
            del Ff_abs, ThetaF, hoft_h1, hoft_l1, hoft_v1, psd_h1_estimated, psd_l1_estimated, psd_v1_estimated, hp_orign, hc_orign, hp_20, hc_20, hp_h1, hc_h1, hp_l1, hc_l1, hp_v1, hc_v1, signal_h1, signal_l1, signal_v1, tilde_signal_h1, tilde_signal_l1, tilde_signal_v1, freq_h1, freq_l1, freq_v1, tilde_signal_h1_lensed, tilde_signal_l1_lensed, tilde_signal_v1_lensed, signal_h1_lensed, signal_l1_lensed, signal_v1_lensed, signal_with_noise_h1, signal_with_noise_l1, signal_with_noise_v1
            gc.collect()
                

            
    return H1_snr_set, L1_snr_set, V1_snr_set, trigger_time_H1, trigger_time_L1, trigger_time_V1, gw_duration, gw_duration_20, Error, dt_h1_from_earth_center, dt_l1_from_earth_center, dt_v1_from_earth_center
        
        # frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/H1/' + str(i) + '.gwf','H1',signal_with_noise_h1)
        # frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/L1/' + str(i) + '.gwf','L1',signal_with_noise_l1)
        # frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/V1/' + str(i) + '.gwf','V1',signal_with_noise_v1)
        # frame.write_frame('./Sim_GW_Data_Lensed_' + Lens_Type + '_only_macro' + '/K1/' + str(i) + '.gwf','K1',signal_with_noise_k1)



if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = 256
    res_z = [0 for i in range(threadcount)]
    length = len(AcceptLensIndex)
    
     
    percount = length//threadcount

    print("percount = ", percount)
   
    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
        if i < threadcount - 1:
            res = pool.apply_async(func=SNR_cal, args=(percount * i, percount * (i + 1),))
            res_z[i] = res
            # print("0",i)
        else:
            # print("1",i)
            res = pool.apply_async(func=SNR_cal, args=(percount * i, length,)) 
            res_z[i] = res

    pool.close()
    pool.join() 
    # Result1 = list(res_z[0].get())
    # Result1 = list(res_z[0].get())
    Result1, Result2, Result3, Result4, Result5, Result6, Result7, Result8, Result9, Result10, Result11, Result12  = list(res_z[0].get()[0]), list(res_z[0].get()[1]), list(res_z[0].get()[2]), list(res_z[0].get()[3]), list(res_z[0].get()[4]), list(res_z[0].get()[5]), list(res_z[0].get()[6]), list(res_z[0].get()[7]), list(res_z[0].get()[8]), list(res_z[0].get()[9]), list(res_z[0].get()[10]), list(res_z[0].get()[11])
    for i in range(threadcount-1):
        Result1.extend(list(res_z[i+1].get()[0]))
        Result2.extend(list(res_z[i+1].get()[1]))
        Result3.extend(list(res_z[i+1].get()[2]))
        Result4.extend(list(res_z[i+1].get()[3]))
        Result5.extend(list(res_z[i+1].get()[4]))
        Result6.extend(list(res_z[i+1].get()[5]))
        Result7.extend(list(res_z[i+1].get()[6]))
        Result8.extend(list(res_z[i+1].get()[7]))
        Result9.extend(list(res_z[i+1].get()[8]))
        Result10.extend(list(res_z[i+1].get()[9]))
        Result11.extend(list(res_z[i+1].get()[10]))
        Result12.extend(list(res_z[i+1].get()[11]))
      
      

SNR_network = np.sqrt(np.array(Result1)**2 + np.array(Result2)**2 + np.array(Result3)**2)       
Index_gtr_12 = np.where(SNR_network>=12)[0]



np.savetxt('./Lensed_SampleResult/H1_snr.csv', Result1, delimiter=',')
np.savetxt('./Lensed_SampleResult/L1_snr.csv', Result2, delimiter=',')
np.savetxt('./Lensed_SampleResult/V1_snr.csv', Result3, delimiter=',')
np.savetxt('./Lensed_SampleResult/trigger_time_H1.csv', Result4, delimiter=',')
np.savetxt('./Lensed_SampleResult/trigger_time_L1.csv', Result5, delimiter=',')
np.savetxt('./Lensed_SampleResult/trigger_time_V1.csv', Result6, delimiter=',')
np.savetxt('./Lensed_SampleResult/gw_duration.csv', Result7, delimiter=',')
np.savetxt('./Lensed_SampleResult/gw_duration_20.csv', Result8, delimiter=',')
np.savetxt('./Lensed_SampleResult/Error_index.csv', Result9, delimiter=',')
np.savetxt('./Lensed_SampleResult/dt_h1_from_earth_center.csv', Result10, delimiter=',')
np.savetxt('./Lensed_SampleResult/dt_l1_from_earth_center.csv', Result11, delimiter=',')      
np.savetxt('./Lensed_SampleResult/dt_v1_from_earth_center.csv', Result12, delimiter=',')      

np.savetxt('./Lensed_SampleResult/SNR_network_only_macro.csv', SNR_network, delimiter=',')


one_image = []
two_image = []
three_image = []
four_image = []
sum_num = 0
for i_ in tqdm(range(len(AcceptLensIndex))):
    start_num = sum_num
    for j_ in range(int(imagenum[i_])):
        sum_num += 1
    end_num = sum_num
    index_tmp = np.where((Index_gtr_12>=start_num)&(Index_gtr_12<end_num))[0]
    lens_num = len(index_tmp)
    if lens_num == 1:
        one_image.append(i_)
    elif lens_num == 2:
        two_image.append(i_)
    elif lens_num == 3:
        three_image.append(i_)
    elif lens_num == 4:
        four_image.append(i_)
   
np.savetxt('./Lensed_SampleResult/two_image_index.csv', two_image, delimiter=',')
np.savetxt('./Lensed_SampleResult/three_image_index.csv', three_image, delimiter=',')
np.savetxt('./Lensed_SampleResult/four_image_index.csv', four_image, delimiter=',') 


#二像、三像、四像的信息。
SNR_network = np.loadtxt('./Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')
imagenum = np.loadtxt('./SampleResult/imagenumber.csv', delimiter=',')

#找到像的index。
two_image = np.loadtxt('./Lensed_SampleResult/two_image_index.csv', delimiter=',')
three_image = np.loadtxt('./Lensed_SampleResult/three_image_index.csv', delimiter=',')
four_image = np.loadtxt('./Lensed_SampleResult/four_image_index.csv', delimiter=',') 
try:
    length_four = len(four_image)
except TypeError:
    four_image = np.array([four_image])
try:
    length_three = len(three_image)
except TypeError:
    three_image = np.array([three_image])
#两个像
two_kappa = []
two_gamma = []
two_mag = []
two_td = []

for i__ in tqdm(range(len(two_image))):
    i_ = int(two_image[i__])
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = imagenum[i_]
    sum_num = int(np.sum(imagenum[0:i_]))
    
    for j_ in range(int(imag_num)):
        td_second = timedelay[sum_num] * 24 * 3600 
        if SNR_network[sum_num] >= 12:
            two_kappa.append(kappa[sum_num])
            two_gamma.append(gamma[sum_num])
            two_mag.append(magnification[sum_num])
            two_td.append(td_second)
            
            

            
            
        sum_num += 1
    

np.savetxt('./Lensed_SampleResult/two_image_kappa.csv', two_kappa, delimiter=',')
np.savetxt('./Lensed_SampleResult/two_image_gamma.csv', two_gamma, delimiter=',')
np.savetxt('./Lensed_SampleResult/two_image_mag.csv', two_mag, delimiter=',') 
np.savetxt('./Lensed_SampleResult/two_image_td.csv', two_td, delimiter=',') 

        






#三个像
three_kappa = []
three_gamma = []
three_mag = []
three_td = []

for i__ in tqdm(range(len(three_image))):
    i_ = int(three_image[i__])
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = imagenum[i_]
    sum_num = int(np.sum(imagenum[0:i_]))
    
    for j_ in range(int(imag_num)):
        
        td_second = timedelay[sum_num] * 24 * 3600
        
       
        if SNR_network[sum_num] >= 12:
            
            three_kappa.append(kappa[sum_num])
            three_gamma.append(gamma[sum_num])
            three_mag.append(magnification[sum_num])
            three_td.append(td_second)
            
            
            
            
        sum_num += 1
    

np.savetxt('./Lensed_SampleResult/three_image_kappa.csv', three_kappa, delimiter=',')
np.savetxt('./Lensed_SampleResult/three_image_gamma.csv', three_gamma, delimiter=',')
np.savetxt('./Lensed_SampleResult/three_image_mag.csv', three_mag, delimiter=',') 
np.savetxt('./Lensed_SampleResult/three_image_td.csv', three_td, delimiter=',') 

#四个像
four_kappa = []
four_gamma = []
four_mag = []
four_td = []

try:
    length_four = len(four_image)
except TypeError:
    four_image = np.array([four_image])
for i__ in tqdm(range(len(four_image))):
    i_ = int(four_image[i__])
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = imagenum[i_]
    sum_num = int(np.sum(imagenum[0:i_]))
    
    for j_ in range(int(imag_num)):
        
        td_second = timedelay[sum_num] * 24 * 3600
        
        if SNR_network[sum_num] >= 12:
            
            four_kappa.append(kappa[sum_num])
            four_gamma.append(gamma[sum_num])
            four_mag.append(magnification[sum_num])
            four_td.append(td_second)
            
            
            
            
        sum_num += 1
    

np.savetxt('./Lensed_SampleResult/four_image_kappa.csv', four_kappa, delimiter=',')
np.savetxt('./Lensed_SampleResult/four_image_gamma.csv', four_gamma, delimiter=',')
np.savetxt('./Lensed_SampleResult/four_image_mag.csv', four_mag, delimiter=',') 
np.savetxt('./Lensed_SampleResult/four_image_td.csv', four_td, delimiter=',') 
    