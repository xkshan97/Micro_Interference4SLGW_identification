'''
生成所有能被看到的强透镜事例的，不区分2像、3像和4像。
'''
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
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
import struct
import sys
sys.path.append(r"../Paper4_MicroLens_Modify/")
from GetMicroDiffrac import DiffracMicro
import matplotlib.pyplot as plt
import gc

# Generate a PSD using an analytic expression for 
# the full design Advanced LIGO noise curve
f_lower = 10
duration = 4096
sample_rate = 4096
sample_RATE = '4KHZ'
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 0.5)
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
lens_z_set = np.loadtxt('./SampleResult/lens_z.csv', delimiter=',')


SNR_network = np.loadtxt('./Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')
Index_SNR_gtr_12 = np.where(SNR_network>12)[0]


four_image = np.loadtxt('./Lensed_SampleResult/four_image_index.csv', delimiter=',') 
four_kappa = np.loadtxt('./Lensed_SampleResult/four_image_kappa.csv', delimiter=',')
four_gamma = np.loadtxt('./Lensed_SampleResult/four_image_gamma.csv', delimiter=',')
four_mag = np.loadtxt('./Lensed_SampleResult/four_image_mag.csv', delimiter=',') 
four_td = np.loadtxt('./Lensed_SampleResult/four_image_td.csv', delimiter=',')

num_min = 258
num_sad = 258
f1=open('../Paper4_MicroLens_Modify/ResultMinimum/TimeLength_quadruple.bin',"rb")
TimeLengthFileMinimum=struct.unpack("l"*num_min, f1.read(8*num_min))
f1.close()
TimeLengthFileMinimum = np.array(TimeLengthFileMinimum)

f2=open('../Paper4_MicroLens_Modify/ResultSaddle/TimeLength_quadruple.bin',"rb")
TimeLengthFileSaddle=struct.unpack("l"*num_sad, f2.read(8*num_sad))
f2.close()
TimeLengthFileSaddle = np.array(TimeLengthFileSaddle)

f3=open('../Paper4_MicroLens_Modify/ResultSaddle/X10_quadruple.bin',"rb")
X10File=struct.unpack("d"*num_sad, f3.read(8*num_sad))
f3.close()
X10File = np.array(X10File)

f4=open('../Paper4_MicroLens_Modify/ResultSaddle/X20_quadruple.bin',"rb")
X20File=struct.unpack("d"*num_sad, f4.read(8*num_sad))
f4.close()
X20File = np.array(X20File)



'''这里做一个查找表出来，给定一个透镜事例的imageindex，能找到在TimeLengthFile中的什么位置'''
LookUpTableMinimum = []
LookUpTableSaddle = []
for i_ in range(len(AcceptLensIndex)):
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = int(imagenum[i_])
    save_index = 0
    IndexSumimage = int(np.sum(imagenum[0:i_]))
    for j_ in range(imag_num):
        
        if SNR_network[IndexSumimage] > 12 and i_ in four_image:
            
            # print(SNR_network[IndexSumimage])
            sqrt_mag_abs = np.sqrt(np.abs(magnification[IndexSumimage]))
            if 1 - kappa[IndexSumimage] - gamma[IndexSumimage] > 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] > 0:    
                """================================="""
                Lens_Type = "Minimum"
                LookUpTableMinimum.append(IndexSumimage)
            elif (1 - kappa[IndexSumimage] - gamma[IndexSumimage]) * (1 - kappa[IndexSumimage] + gamma[IndexSumimage]) < 0:
                Lens_Type = "Saddle"
                LookUpTableSaddle.append(IndexSumimage)
        IndexSumimage += 1
                
LookUpTableMinimum = np.array(LookUpTableMinimum)                
LookUpTableSaddle = np.array(LookUpTableSaddle)   
                    

def GenWaveFile(start_index, end_index): 
    inject_time = []
    
    trigger_time_H1 = []
    trigger_time_L1 = []
    trigger_time_V1 = []
  
    gw_duration = [] #10Hz的持续时间
    gw_duration_20 = [] #20Hz的持续时间
    dt_h1_from_earth_center = []
    dt_l1_from_earth_center = []
    dt_v1_from_earth_center = []
    
    for i_ in tqdm(range(start_index, end_index)):
        i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
        imag_num = int(imagenum[i_])
        save_index = 0
        IndexSumimage = int(np.sum(imagenum[0:i_]))
        for j_ in range(imag_num):
            
            if SNR_network[IndexSumimage] > 12 and i_ in four_image:
                
                td_second = timedelay[IndexSumimage] * 24 * 3600
                SampleParam9_td = SampleParam[9][i] + td_second
            
                # print(SNR_network[IndexSumimage])
                sqrt_mag_abs = np.sqrt(np.abs(magnification[IndexSumimage]))
                if 1 - kappa[IndexSumimage] - gamma[IndexSumimage] > 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] > 0:    
                    """================================="""
                    Lens_Type = "Minimum"
                    num_minimum = np.where(LookUpTableMinimum == IndexSumimage)[0][0]
                elif (1 - kappa[IndexSumimage] - gamma[IndexSumimage]) * (1 - kappa[IndexSumimage] + gamma[IndexSumimage]) < 0:
                    Lens_Type = "Saddle"
                    num_saddle = np.where(LookUpTableSaddle == IndexSumimage)[0][0]
                
                if kappa[IndexSumimage] not in four_kappa or gamma[IndexSumimage] not in four_gamma or td_second not in four_td:
                    print('wrong !!!!!!!!!!!!!!!!!!!')
                    
                
                """透镜的放大"""
                if Lens_Type == "Minimum": 
                    M_L = 0.38
                    z_L = lens_z_set[i_]
                    fileArea = "../Paper4_MicroLens_Modify/ResultMinimum/Total" + str(i) + "_" + str(save_index) + "ReadAreaMinimum_quadruple.bin"
                    fileTime = "../Paper4_MicroLens_Modify/ResultMinimum/Total" + str(i) + "_" + str(save_index) + "ReadTimeMinimum_quadruple.bin"
                    binlength = TimeLengthFileMinimum[num_minimum]
                    num_minimum += 1
                    x_step = 0.005
                    Ff_abs, ThetaF = DiffracMicro().Minimum(M_L, z_L, kappa[IndexSumimage], gamma[IndexSumimage], fileArea, fileTime, binlength, x_step, lens_freq) 
                    
                    
                elif Lens_Type == "Saddle":
                    M_L = 0.38
                    z_L = lens_z_set[i_]
                    fileArea = "../Paper4_MicroLens_Modify/ResultSaddle/Total" + str(i) + "_" + str(save_index) + "ReadAreaSaddle_quadruple.bin"
                    fileTime = "../Paper4_MicroLens_Modify/ResultSaddle/Total" + str(i) + "_" + str(save_index) + "ReadTimeSaddle_quadruple.bin"
                    binlength = TimeLengthFileSaddle[num_saddle]
                    x_step = 0.005
                    x10 = X10File[num_saddle]
                    x20 = X20File[num_saddle]
                    num_saddle += 1
                    Ff_abs, ThetaF = DiffracMicro().Saddle(M_L, z_L, kappa[IndexSumimage], gamma[IndexSumimage], fileArea, fileTime, binlength, x_step, x10, x20, lens_freq)
                    ThetaF[np.where(lens_freq==0)[0][0]] = 0
                    
                
                # plt.figure(IndexSumimage)
                # plt.semilogx(lens_freq, Ff_abs)
                # plt.semilogx([30,1000],[sqrt_mag_abs, sqrt_mag_abs])   

                Ff_inter = interp1d(lens_freq, Ff_abs)
                ThetaF_inter = interp1d(lens_freq, ThetaF)


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
            
                # plt.plot(hoft_l1_h1.sample_times, hoft_l1_h1)
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
                inject_time.append(int(hoft_h1.start_time.gpsSeconds))

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

                
                np.savetxt("../Paper4_MicroLens_Modify/Ffabs" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) + "_quadruple.csv", Ff_inter(freq_h1), delimiter=',')
                np.savetxt("../Paper4_MicroLens_Modify/ThetaF" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) + "_quadruple.csv", ThetaF_inter(freq_h1), delimiter=',')
                Lens_Type += "_" + 'with_micro'
                #画一下加入透镜之前和之后在频域上的不同
                plt.loglog(freq_h1, np.abs(tilde_signal_h1))
                plt.loglog(freq_h1, np.abs(tilde_signal_h1_lensed)/sqrt_mag_abs, label='lensed')
                plt.legend()
                plt.title('SNR = ' + str(round(SNR_network[IndexSumimage],2)))
                plt.xlim(1, 1000)
                plt.savefig('./IntermediaPlot_Gen_inter_Micro_quadruple/Tilde_hf' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                plt.close()
                
                
                #得到透镜以后的时域
                signal_h1_lensed = np.real(scipy.fft.ifft(tilde_signal_h1_lensed) / signal_h1.delta_t)
                signal_l1_lensed = np.real(scipy.fft.ifft(tilde_signal_l1_lensed) / signal_l1.delta_t)
                signal_v1_lensed = np.real(scipy.fft.ifft(tilde_signal_v1_lensed) / signal_v1.delta_t)
                
                #画一下加入透镜之前和之后时域上的不同
                # plt.plot(signal_h1.sample_times, signal_h1)
                # plt.plot(signal_h1.sample_times, signal_h1_lensed, '--', label='lensed')
                # plt.legend()
                # plt.xlim(end_time_h1+hp_h1.duration*3/4, end_time_h1+hp_h1.duration)
                
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
                
                # plt.figure(i)
                plt.plot(signal_h1.sample_times, signal_h1_lensed, 'r', label = 'H1')
                plt.plot(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end], signal_h1_func(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end]), 'r--')
                plt.plot(signal_l1.sample_times, signal_l1_lensed, 'b', label = 'L1')
                plt.plot(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end], signal_l1_func(signal_with_noise_l1.sample_times[signal_with_noise_l1_index_start:signal_with_noise_l1_index_end]), 'b--')
                plt.plot(signal_v1.sample_times, signal_v1_lensed, 'k', label = 'V1')
                plt.plot(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end], signal_v1_func(signal_with_noise_v1.sample_times[signal_with_noise_v1_index_start:signal_with_noise_v1_index_end]), 'k--')
                plt.xlim(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end][len(signal_with_noise_h1.sample_times[signal_with_noise_h1_index_start:signal_with_noise_h1_index_end])*3//4], signal_with_noise_h1.sample_times[signal_with_noise_h1_index_end] + 0.05)
                plt.legend()
                plt.title('SNR = ' + str(SNR_network[IndexSumimage]))
                plt.xlabel('t[s]')
                plt.ylabel('strain')
                plt.savefig('./IntermediaPlot_Gen_inter_Micro_quadruple/waveform_' + str(i) + '_' +str(save_index) +'.png', dpi=450)
                plt.close()
                
                
                frame.write_frame('./Total_Sim_GW_Data_Lensed_' + Lens_Type + '_quadruple/H1/H-H1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_h1.start_time.gpsSeconds) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf','H1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_h1)
                frame.write_frame('./Total_Sim_GW_Data_Lensed_' + Lens_Type + '_quadruple/L1/L-L1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_l1.start_time.gpsSeconds) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf','L1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_l1)
                frame.write_frame('./Total_Sim_GW_Data_Lensed_' + Lens_Type + '_quadruple/V1/V-V1_GWOSC_' + sample_RATE + '_R1-' + str(hoft_v1.start_time.gpsSeconds) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf','V1:GWOSC_' + sample_RATE + '_R1_STRAIN',signal_with_noise_v1)
            
                
                save_index += 1
                del Ff_abs, ThetaF, hoft_h1, hoft_l1, hoft_v1, psd_h1_estimated, psd_l1_estimated, psd_v1_estimated, hp_orign, hc_orign, hp_20, hc_20, hp_h1, hc_h1, hp_l1, hc_l1, hp_v1, hc_v1, signal_h1, signal_l1, signal_v1, tilde_signal_h1, tilde_signal_l1, tilde_signal_v1, freq_h1, freq_l1, freq_v1, tilde_signal_h1_lensed, tilde_signal_l1_lensed, tilde_signal_v1_lensed, signal_h1_lensed, signal_l1_lensed, signal_v1_lensed, signal_with_noise_h1, signal_with_noise_l1, signal_with_noise_v1
                gc.collect()
            IndexSumimage += 1
            
    return inject_time, trigger_time_H1, trigger_time_L1, trigger_time_V1, gw_duration, gw_duration_20, dt_h1_from_earth_center, dt_l1_from_earth_center, dt_v1_from_earth_center
# SNR_network_2 = np.sqrt(np.array(H1_snr_set)**2 + np.array(L1_snr_set)**2)

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
            res = pool.apply_async(func=GenWaveFile, args=(percount * i, percount * (i + 1),))
            res_z[i] = res
            # print("0",i)
        else:
            # print("1",i)
            res = pool.apply_async(func=GenWaveFile, args=(percount * i, length,)) 
            res_z[i] = res
            
    pool.close()
    pool.join() 
    Result1, Result2, Result3, Result4, Result5, Result6, Result7, Result8, Result9 = list(res_z[0].get()[0]), list(res_z[0].get()[1]), list(res_z[0].get()[2]), list(res_z[0].get()[3]), list(res_z[0].get()[4]), list(res_z[0].get()[5]), list(res_z[0].get()[6]), list(res_z[0].get()[7]), list(res_z[0].get()[8])
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

np.savetxt('./Lensed_SampleResult_Micro_quadruple/Total_inject_time.csv', Result1)
np.savetxt('./Lensed_SampleResult_Micro_quadruple/trigger_time_H1.csv', Result2, delimiter=',')
np.savetxt('./Lensed_SampleResult_Micro_quadruple/trigger_time_L1.csv', Result3, delimiter=',')
np.savetxt('./Lensed_SampleResult_Micro_quadruple/trigger_time_V1.csv', Result4, delimiter=',')
np.savetxt('./Lensed_SampleResult_Micro_quadruple/gw_duration.csv', Result5, delimiter=',')
np.savetxt('./Lensed_SampleResult_Micro_quadruple/gw_duration_20.csv', Result6, delimiter=',')
np.savetxt('./Lensed_SampleResult_Micro_quadruple/dt_h1_from_earth_center.csv', Result7, delimiter=',')
np.savetxt('./Lensed_SampleResult_Micro_quadruple/dt_l1_from_earth_center.csv', Result8, delimiter=',')      
np.savetxt('./Lensed_SampleResult_Micro_quadruple/dt_v1_from_earth_center.csv', Result9, delimiter=',')      
      

# pylab.plot(hoft_l1_h1.sample_times, hoft_l1_h1,'b')
# plt.plot(signal_h1.sample_times, signal_h1)
# plt.ylabel('Strain')
# plt.xlabel('Time (s)')
# plt.legend()
# plt.show()
#画图检验
'''
SNR_selected = np.loadtxt('../Paper5_cWB_Modify/Total_inject_time_and_snr_and_match_and_outdir_name/SNR_selected' + f_star + '.csv', delimiter=',')
IndexSumimage = 0
for i_ in range(len(AcceptLensIndex)):
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = int(imagenum[i_])
    save_index = 0
    for j_ in range(imag_num):
        if SNR_network[IndexSumimage] > 12:
            if SNR_network[IndexSumimage] in SNR_selected:
                sqrt_mag_abs = np.sqrt(np.abs(magnification[IndexSumimage]))
                if 1 - kappa[IndexSumimage] - gamma[IndexSumimage] > 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] > 0:    
                    """================================="""
                    Lens_Type = "Minimum"
                elif (1 - kappa[IndexSumimage] - gamma[IndexSumimage]) * (1 - kappa[IndexSumimage] + gamma[IndexSumimage]) < 0:
                    Lens_Type = "Saddle"
                elif 1 - kappa[IndexSumimage] - gamma[IndexSumimage] < 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] < 0:
                    Lens_Type = "Maximum"
                    
                time = np.loadtxt("../Paper4_MicroLens_Modify/timenew" + Lens_Type + "/timenew_Total" + str(i) + "_" + str(save_index) +f_star + ".csv",  delimiter=',')
                Ft = np.loadtxt("../Paper4_MicroLens_Modify/Ftnew" + Lens_Type + "/Ftnew_Total" + str(i) + "_" + str(save_index) + f_star +".csv", delimiter=',')
                
                Ffabs = np.loadtxt("../Paper4_MicroLens_Modify/Ffabs" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) + f_star +  ".csv",  delimiter=',')
                ThetaF = np.loadtxt("../Paper4_MicroLens_Modify/ThetaF" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) + f_star +  ".csv", delimiter=',')
                plt.figure(IndexSumimage)
                plt.plot(Ffabs[0:len(Ffabs)//2])
                # # # plt.plot(Ffabs)
                plt.plot([0,2048],[sqrt_mag_abs, sqrt_mag_abs])
                # plt.semilogx(time, Ft, label=str(i))
                # plt.legend()
                # plt.plot(ThetaF[1::])
            save_index += 1
        
        IndexSumimage += 1

    
'''

