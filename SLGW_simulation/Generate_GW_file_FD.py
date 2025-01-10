'''
修改成10Hz了，然后运行的文件没有全部运行，而是只运行了1068个
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

# Generate a PSD using an analytic expression for 
# the full design Advanced LIGO noise curve
f_lower = 10
duration = 4096
sample_rate = 4096
sample_RATE = '4KHZ'
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)

SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

apx = 'IMRPhenomPv2'
"""=================================="""
"""画图"""
'''
# Let's take a look at the spectrum
pylab.loglog(psd_l1_h1.sample_frequencies, psd_l1_h1**0.5, label = 'L1 & H1')
pylab.xlim(20, 1024) 
pylab.legend()
pylab.xlabel('Frequency (Hz)')
pylab.ylabel('Strain^2 / Hz')
pylab.grid()
'''
"""=================================="""

# Now, let's generate noise that has the same spectrum
htilde_l1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=1000)
htilde_h1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=2000)
htilde_v1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=3000)
# pylab.loglog(htilde_l1.sample_frequencies, np.abs(htilde_l1))
# pylab.loglog(psd_l1_h1.sample_frequencies, psd_l1_h1**0.5, label = 'L1 & H1')
        

def Generate_file(IndexStart, IndexEnd):
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
    
    for i in tqdm(range(IndexStart, IndexEnd)):
                
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
        
         #验证一下welch方法得到的结果
        # pp.loglog(psd_h1_estimated.sample_frequencies, psd_h1_estimated ** 0.5, label='Estimated')
        # pylab.loglog(psd_l1_h1.sample_frequencies, psd_l1_h1**0.5, label = 'L1 & H1') 
        
        # try:
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
        start_time_signal = int(SampleParam[9][i]) - (hp_orign.duration - hp_20.duration) #保证20Hz的开始时间是在start time上
        start_time = int(SampleParam[9][i])
        
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
        
        #得到再逆傅立叶变换回去，为的是跟后面透镜的操作一致，也为了验证一下后面的透镜操作对不对。
        signal_h1_re = np.real(scipy.fft.ifft(tilde_signal_h1) / signal_h1.delta_t)
        signal_l1_re = np.real(scipy.fft.ifft(tilde_signal_l1) / signal_l1.delta_t)
        signal_v1_re = np.real(scipy.fft.ifft(tilde_signal_v1) / signal_v1.delta_t)
        
        
        signal_h1_func = interp1d(signal_h1.sample_times, signal_h1_re)
        signal_l1_func = interp1d(signal_l1.sample_times, signal_l1_re)
        signal_v1_func = interp1d(signal_v1.sample_times, signal_v1_re)
            
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
        
        plt.figure(i)
        plt.plot(signal_h1.sample_times, signal_h1.data, 'r', label = 'H1')
        plt.plot(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end], signal_h1_func(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end]), 'r--')
        plt.plot(signal_l1.sample_times, signal_l1.data, 'b', label = 'L1')
        plt.plot(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end], signal_l1_func(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end]), 'b--')
        plt.plot(signal_v1.sample_times, signal_v1.data, 'k', label = 'V1')
        plt.plot(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end], signal_v1_func(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end]), 'k--')
        plt.xlim(signal_h1.sample_times[len(signal_h1.sample_times)*1//2],signal_h1.sample_times[len(signal_h1.sample_times) - 1] + 0.1)
        plt.legend()
        plt.xlabel('t[s]')
        plt.ylabel('strain')
        plt.savefig('./IntermediaPlot_Gen_inter/waveform_' + str(i) + '.png', dpi=450)
        plt.close()
    
        #===============================
        """matched filter SNR"""
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
    
        
    
        frame.write_frame('./Sim_GW_Data/H1/H-H1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_h1.start_time.gpsSeconds) + '-4096_' + str(i) + '.gwf','H1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_h1)
        frame.write_frame('./Sim_GW_Data/L1/L-L1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_l1.start_time.gpsSeconds) + '-4096_' + str(i) + '.gwf','L1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_l1)
        frame.write_frame('./Sim_GW_Data/V1/V-V1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_v1.start_time.gpsSeconds) + '-4096_' + str(i) + '.gwf','V1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_v1)
        
            
    
    return H1_snr_set, L1_snr_set, V1_snr_set, trigger_time_H1, trigger_time_L1, trigger_time_V1, gw_duration, gw_duration_20, Error, dt_h1_from_earth_center, dt_l1_from_earth_center, dt_v1_from_earth_center
# pylab.plot(hoft_l1_h1.sample_times, hoft_l1_h1,'b')
# pp.plot(signal_h1.sample_times, signal_h1)
# pp.ylabel('Strain')
# pp.xlabel('Time (s)')
# pp.legend()
# pp.show()

# Generate_file(0,1)
if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = 267
    res_z = [0 for i in range(threadcount)]
    length = int(1068/4) #len(SampleParam[0])
    
     
    percount = length//threadcount

    print("percount = ", percount)
   
    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
        if i < threadcount - 1:
            res = pool.apply_async(func=Generate_file, args=(percount * i, percount * (i + 1),))
            res_z[i] = res
            # print("0",i)
        else:
            # print("1",i)
            res = pool.apply_async(func=Generate_file, args=(percount * i, length,)) 
            res_z[i] = res

    pool.close()
    pool.join() 
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
        
       
    
np.savetxt('./Sim_GW_Data/H1/H1_snr.csv', Result1, delimiter=',')
np.savetxt('./Sim_GW_Data/L1/L1_snr.csv', Result2, delimiter=',')
np.savetxt('./Sim_GW_Data/V1/V1_snr.csv', Result3, delimiter=',')
np.savetxt('./Sim_GW_Data/trigger_time_H1.csv', Result4, delimiter=',')
np.savetxt('./Sim_GW_Data/trigger_time_L1.csv', Result5, delimiter=',')
np.savetxt('./Sim_GW_Data/trigger_time_V1.csv', Result6, delimiter=',')
np.savetxt('./Sim_GW_Data/gw_duration.csv', Result7, delimiter=',')
np.savetxt('./Sim_GW_Data/gw_duration_20.csv', Result8, delimiter=',')
np.savetxt('./Sim_GW_Data/Error_index.csv', Result9, delimiter=',')
np.savetxt('./Sim_GW_Data/dt_h1_from_earth_center.csv', Result10, delimiter=',')
np.savetxt('./Sim_GW_Data/dt_l1_from_earth_center.csv', Result11, delimiter=',')      
np.savetxt('./Sim_GW_Data/dt_v1_from_earth_center.csv', Result12, delimiter=',')      