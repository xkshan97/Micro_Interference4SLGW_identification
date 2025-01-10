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
import pycbc.noise

apx = 'IMRPhenomPv2'

"""match function"""
def Match(strain1, strain2):
    cor_match = 4 * np.sum(strain1 * np.conjugate(strain2))
    match1 = 4 * np.sum(strain1 * np.conjugate(strain1))
    match2 = 4 * np.sum(strain2 * np.conjugate(strain2))
    return np.abs(cor_match / np.sqrt(match1 * match2))

f_lower = 30
duration = 4096 #"""注意这里改成了4096，paper4_copy中是64."""
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_ligo = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)
func_psd = interp1d(psd_ligo.sample_frequencies, psd_ligo.data)
htilde_h1 = pycbc.noise.frequency_noise_from_psd(psd_ligo, seed=2000)
hoft_h1 = htilde_h1.to_timeseries()


inject_time_set = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/inject_time_set.csv', delimiter=',')
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name_complex/max_match_det_lensed.csv', 'r') as f1:
    reference_det = f1.readlines()
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/with_micro.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据
freq_max = np.loadtxt('./cWB_micro/Lensed_freq.csv', delimiter=',', encoding='UTF-8-sig')
#上面是从图上用眼睛读出来的数据
match_noise = []
match_no_noise = []
for i in tqdm(range(len(inject_time_set))):

    match_micro = []
    reference_det_tmp = reference_det[i].split('\n')[0]
    try:
        with open(file_name[i].split('\n')[0]) as f1:
            data_lensed = json.load(f1)
        
        parameters = list(data_lensed['posterior']['content'].keys())
    except FileNotFoundError:
        
        continue
    
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
    iota = data_lensed['posterior']['content']['iota'][j]
    m1_o = Mz_o * (1 + q)**(1/5) / q ** (3/5)
    m2_o = q * m1_o
    
    hp, hc = get_td_waveform(approximant=apx,
                            mass1=m1_o, mass2=m2_o,
                            spin1x = spin_1x, spin1y = spin_1y,
                            spin1z = spin_1z, spin2x = spin_2x,
                            spin2y = spin_2y, spin2z = spin_2z,
                            inclination = iota,
                            distance = luminosity_dis,
                            delta_t=1.0/sample_rate, f_lower=10, f_ref = 50, coa_phase = phase)
    det_h1 = Detector(reference_det_tmp)
    
    
    
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
    
    #without noise
    tilde_signal_h1 = scipy.fft.fft(signal_h1.data) * np.diff(signal_h1.sample_times)[0]

    freq_h1 = scipy.fft.fftfreq(len(signal_h1.data), np.diff(signal_h1.sample_times)[0])



    Bilby_func_h1 = interp1d(freq_h1, tilde_signal_h1)
    

    #with noise
    tilde_signal_h1_noise = scipy.fft.fft(signal_h1.data + hoft_h1.data[0:len(signal_h1.data)]) * np.diff(signal_h1.sample_times)[0]
    Bilby_func_h1_noise = interp1d(freq_h1, tilde_signal_h1_noise)
    
    
    #频率范围
    index_yes = np.where((htilde_h1.sample_frequencies>f_lower)&(htilde_h1.sample_frequencies<freq_max[i]))

    freq_match = htilde_h1.sample_frequencies[index_yes]
    
    inject_freq = np.loadtxt('./Micro_signal_Total/freq_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    inject_signal = np.loadtxt('./Micro_signal_Total/lens_signal_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    inject_time = np.loadtxt('./Micro_signal_Total/time_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    time_inject_signal = np.loadtxt('./Micro_signal_Total/time_lens_signal_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    
    func_inject = interp1d(inject_freq, inject_signal)
    
    match_tmp_h1_noise = Match((func_inject(freq_match))/func_psd(freq_match)**(1/2), (np.abs(Bilby_func_h1_noise(freq_match)))/func_psd(freq_match)**(1/2))
    match_tmp_h1 = Match((func_inject(freq_match))/func_psd(freq_match)**(1/2), (np.abs(Bilby_func_h1(freq_match)))/func_psd(freq_match)**(1/2))
    
    match_noise.append(match_tmp_h1_noise)
    match_no_noise.append(match_tmp_h1)
    
    
    
    plt.semilogy(freq_match, np.abs(Bilby_func_h1_noise(freq_match)), label='noisy')
    plt.semilogy(freq_h1, np.abs(tilde_signal_h1), label='Bilby')
    plt.semilogy(freq_match, np.abs(Bilby_func_h1(freq_match)), label='Bilby inter')
    plt.semilogy(inject_freq, inject_signal, '--',label='injected')
    # plt.semilogy(freq_match, func_inject(freq_match), '--',label='injected')
    # plt.plot(freq_match, np.abs(htilde_h1.data[index_yes]))
    # plt.semilogy(freq_match, np.abs(htilde_h1.data[index_yes]+Bilby_func_h1(freq_match)))
    # plt.plot(freq_match, func_psd(freq_match)**0.5)
    plt.xlim(10, 200)
    plt.ylim(10**(-26), 10**(-22))
    plt.legend()
    plt.grid()
    plt.savefig('IntermediaPlot/Micro_cWB_inject/noisy_freq_' + file_name[i].split('/')[3][6::] + '.png', dpi=450)
    plt.close()
    
np.savetxt('./match_noisy_inject_VS_Bilby/micro_match_noise.csv', match_noise, delimiter=',')
np.savetxt('./match_noisy_inject_VS_Bilby/micro_match_no_noise.csv', match_no_noise, delimiter=',')
    
    
