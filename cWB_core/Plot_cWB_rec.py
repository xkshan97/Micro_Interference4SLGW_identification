'''
读取Paper5_150914like中的文件（这个文件夹下有很多算自相关的程序），与cWB X Bilby相关的
只有Generate_150913_like.py
ReadAndEstimated_lens.py
cal_match_lensed.py
'''
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import scipy
from scipy.interpolate import interp1d
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
from astropy.cosmology import Planck18
from tqdm import tqdm
from  multiprocessing import Process,Pool
import time
import pycbc.psd


"读取cWB重构结果"
with open('../Paper5_150914like/cWB_rec/GW150914/L1H1_1126259458.438_1126259458.438/H1_wf_strain.dat') as f3:
    data = f3.readlines()
h1t_lensed_rec_t = []
h1t_lensed_rec = []
for i in range(len(data)):
    h1t_lensed_rec_t.append(float(data[i].split()[0]))
    h1t_lensed_rec.append(float(data[i].split()[1]))
h1f_lensed_rec = scipy.fft.fft(h1t_lensed_rec) * np.diff(h1t_lensed_rec_t)[0]
h1f_lensed_rec_f = scipy.fft.fftfreq(len(h1t_lensed_rec), np.diff(h1t_lensed_rec_t)[0])
func_lensed_rec = interp1d(h1f_lensed_rec_f, abs(h1f_lensed_rec))

index_yes = np.where((h1f_lensed_rec_f>40)&(h1f_lensed_rec_f<300))

freq_plot_rude = h1f_lensed_rec_f[index_yes]
freq_plot = np.arange(freq_plot_rude[0], freq_plot_rude[-1], 10**(-3))
'''以上：cWB'''

H1_inject_lensed = np.loadtxt('../Paper5_150914like/GW150914/H1/2image21626_00.090signal.csv', delimiter=',')
H1_freq = np.loadtxt('../Paper5_150914like/GW150914/H1/2image21626_00.090signal_freq.csv', delimiter=',')
func_inject = interp1d(H1_freq, H1_inject_lensed)


'''以上：输入的透镜波形'''    

Ffabs_freq = np.loadtxt('../Paper5_150914like/GW150914/freq/2image21626_00.090.csv', delimiter=',')
Ffabs = np.loadtxt('../Paper5_150914like/GW150914/Ffabs/2image21626_00.090.csv', delimiter=',')
func_Ffabs = interp1d(Ffabs_freq, Ffabs)
'''微透镜的Ffabs'''

"读取Bilby重构结果"
f_lower = 30
duration = 4096 #"""注意这里改成了4096，paper4_copy中是64."""
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)
func_psd = interp1d(psd_l1_h1.sample_frequencies, psd_l1_h1.data)

file1_lensed = '../Paper5_150914like/outdir_Lensed/outdir21626/H1_L1_result.json'
    
with open(file1_lensed) as f1:
    data_lensed = json.load(f1)

parameters = list(data_lensed['posterior']['content'].keys())
apx = 'IMRPhenomPv2'
SNR_tmp = []
match_micro = []


j = np.where(data_lensed['posterior']['content']['log_likelihood'] == np.max(data_lensed['posterior']['content']['log_likelihood']))[0][0]
q = data_lensed['posterior']['content']['mass_ratio'][j]
Mz_o = data_lensed['posterior']['content']['chirp_mass'][j]
luminosity_dis = data_lensed['posterior']['content']['luminosity_distance'][j]
dec = data_lensed['posterior']['content']['dec'][j]
ra = data_lensed['posterior']['content']['ra'][j]
theta_jn = data_lensed['posterior']['content']['theta_jn'][j]
psi = data_lensed['posterior']['content']['psi'][j]
phase = data_lensed['posterior']['content']['phase'][j]
a_1 = data_lensed['posterior']['content']['a_1'][j]
a_2 = data_lensed['posterior']['content']['a_2'][j]
spin_1x = data_lensed['posterior']['content']['spin_1x'][j]
spin_1y = data_lensed['posterior']['content']['spin_1y'][j]
spin_1z = data_lensed['posterior']['content']['spin_1z'][j]
spin_2x = data_lensed['posterior']['content']['spin_2x'][j]
spin_2y = data_lensed['posterior']['content']['spin_2y'][j]
spin_2z = data_lensed['posterior']['content']['spin_2z'][j]

tilt_1 = data_lensed['posterior']['content']['tilt_1'][j]
tilt_2 = data_lensed['posterior']['content']['tilt_2'][j]
phi_12 = data_lensed['posterior']['content']['phi_12'][j]
phi_jl = data_lensed['posterior']['content']['phi_jl'][j]
# time_jitter = data_lensed['posterior']['content']['time_jitter'][j]
geocent_time = data_lensed['posterior']['content']['geocent_time'][j]
m1_o = Mz_o * (1 + q)**(1/5) / q ** (3/5)
m2_o = q * m1_o
SNR_tmp.append(np.sqrt(data_lensed['posterior']['content']['H1_matched_filter_snr'][j]['real']**2 + data_lensed['posterior']['content']['L1_matched_filter_snr'][j]['real']**2))


hp, hc = get_td_waveform(approximant=apx,
                        mass1=m1_o, mass2=m2_o,
                        spin1x = spin_1x, spin1y = spin_1y,
                        spin1z = spin_1z, spin2x = spin_2x,
                        spin2y = spin_2y, spin2z = spin_2z,
                        inclination = theta_jn,
                        distance = luminosity_dis,
                        delta_t=1.0/sample_rate, f_lower=10, f_ref = 50, coa_phase = phase)


det_h1 = Detector('H1')



hp.start_time += geocent_time
hc.start_time += geocent_time
# 从地心到探测器的时间
dt_h1_re = det_h1.time_delay_from_earth_center(ra, dec, hp.start_time.gpsSeconds + hp.start_time.gpsNanoSeconds / 10 ** 9)
            
hp.start_time += dt_h1_re
hc.start_time += dt_h1_re

hp_h1, hc_h1 = hp.copy(), hc.copy()

# We get back the fp and fc antenna pattern weights.
fp_h1_re, fc_h1_re = det_h1.antenna_pattern(ra, dec, psi, geocent_time)

signal_h1 = fp_h1_re * hp_h1 + fc_h1_re * hc_h1

#得到频域的信号

tilde_signal_h1 = scipy.fft.fft(signal_h1.data) * signal_h1.delta_t

freq_h1 = scipy.fft.fftfreq(len(signal_h1), signal_h1.delta_t)


Bilby_func_h1 = interp1d(freq_h1, abs(tilde_signal_h1))


'''以上：Bilby重构的'''
plt.figure(i)
plt.title(i)
plt.semilogy(freq_plot, abs(func_lensed_rec(freq_plot)), label = 'cWB')
plt.semilogy(freq_plot, func_inject(freq_plot), label = 'inject')
plt.semilogy(freq_plot, Bilby_func_h1(freq_plot), label = 'Bilby')
# plt.xlim(40,250)
plt.ylim(10**(-25.5),10**(-23.5))
# plt.xlim(40,300)
plt.grid()

'''画一下Ffabs和cWB/Bilby'''
plt.plot(freq_plot, func_Ffabs(freq_plot)/Ffabs[0])
plt.plot(freq_plot, abs(func_lensed_rec(freq_plot))/Bilby_func_h1(freq_plot))
plt.ylim(0,2.5)

np.savetxt('Result_file4plot/freq_plot.csv', freq_plot, delimiter=',')
np.savetxt('Result_file4plot/cWB_hf.csv', abs(func_lensed_rec(freq_plot)), delimiter=',')
np.savetxt('Result_file4plot/inject_hf.csv', func_inject(freq_plot), delimiter=',')
np.savetxt('Result_file4plot/Bilby_hf.csv', abs(Bilby_func_h1(freq_plot)), delimiter=',')
np.savetxt('Result_file4plot/Ffabs.csv', func_Ffabs(freq_plot), delimiter=',')
