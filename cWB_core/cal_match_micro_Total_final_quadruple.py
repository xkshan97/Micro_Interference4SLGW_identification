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


apx = 'IMRPhenomPv2'

"""match function"""
def Match(strain1, strain2):
    cor_match = 4 * np.sum(strain1 * strain2)
    match1 = 4 * np.sum(strain1 * strain1)
    match2 = 4 * np.sum(strain2 * strain2)
    return cor_match / np.sqrt(match1 * match2)

f_lower = 30
duration = 4096 #"""注意这里改成了4096，paper4_copy中是64."""
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_ligo = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)
func_psd = interp1d(psd_ligo.sample_frequencies, psd_ligo.data)

inject_time_set = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_quadruple/inject_time_set.csv', delimiter=',')
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name_quadruple/max_match_det_lensed.csv', 'r') as f1:
    reference_det = f1.readlines()
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name_quadruple/with_micro.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据
freq_max = np.loadtxt('./cWB_micro_quadruple/quad_freq.csv', delimiter=',', encoding='UTF-8-sig')
#上面是从图上用眼睛读出来的数据
match_filter_SNR = []
max_match_det_set = []
for i in tqdm(range(len(inject_time_set))):
    
    match_micro = []
    reference_det_tmp = reference_det[i].split('\n')[0]
    try:
        with open(file_name[i].split('\n')[0]) as f1:
            data_lensed = json.load(f1)
        
        parameters = list(data_lensed['posterior']['content'].keys())
    except FileNotFoundError:
        
        continue
        
    
    
    
    

    
    try:
        "读取cWB重构结果"
        cWB_path = './cWB_micro_quadruple/data/' + str(int(inject_time_set[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/' + reference_det_tmp + '_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
    except FileNotFoundError:
        
        continue
    
        
    h1t_rec_t = []
    h1t_rec = []
    for k in range(len(data)):
        h1t_rec_t.append(float(data[k].split()[0]))
        h1t_rec.append(float(data[k].split()[1]))
    h1f_rec = scipy.fft.fft(h1t_rec) * np.diff(h1t_rec_t)[0]
    h1f_rec_f = scipy.fft.fftfreq(len(h1t_rec), np.diff(h1t_rec_t)[0])
    func_rec = interp1d(h1f_rec_f, abs(h1f_rec))
    
    index_yes = np.where((h1f_rec_f>f_lower)&(h1f_rec_f<freq_max[i]))

    freq_match_rude = h1f_rec_f[index_yes]
    freq_match = np.arange(freq_match_rude[0], freq_match_rude[-1], 10**(-3))
    

    
    start_index = 0
    end_index = 100
    for j in range(start_index, end_index):
        try:
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
            
            hp, hc = get_td_waveform(approximant=apx,
                                    mass1=m1_o, mass2=m2_o,
                                    spin1x = spin_1x, spin1y = spin_1y,
                                    spin1z = spin_1z, spin2x = spin_2x,
                                    spin2y = spin_2y, spin2z = spin_2z,
                                    inclination = theta_jn,
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


            #得到频域的信号

            tilde_signal_h1 = scipy.fft.fft(signal_h1.data) * signal_h1.delta_t
            
            freq_h1 = scipy.fft.fftfreq(len(signal_h1), signal_h1.delta_t)
            
            
            Bilby_func_h1 = interp1d(freq_h1, abs(tilde_signal_h1))
            
            
            
            match_tmp_h1 = Match(abs(func_rec(freq_match))/func_psd(freq_match)**(1/2), abs(Bilby_func_h1(freq_match))/func_psd(freq_match)**(1/2))
            
            
                
                
                
            match_micro.append(match_tmp_h1)
        except ValueError:
            match_micro.append(match_tmp_h1)
    
    
    
    np.savetxt('./match_result_Total_quadruple/match_micro_'+file_name[i].split('/')[3][6::] + '.csv', match_micro, delimiter=',')
    
    
