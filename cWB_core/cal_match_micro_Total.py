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

inject_time_set = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/inject_time_set.csv', delimiter=',')
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/reference_det.csv', 'r') as f1:
    reference_det = f1.readlines()
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/with_micro.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据
freq_max = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/ringdown_freq.csv', delimiter=',', encoding='UTF-8-sig') + 100
freq_max = np.loadtxt('./cWB_micro/Lensed_freq.csv', delimiter=',', encoding='UTF-8-sig')
#上面是从图上用眼睛读出来的数据
match_filter_SNR = []
max_match_det_set = []
for i in tqdm(range(len(inject_time_set))):
    max_match = 0
    # SNR_tmp = []
    match_micro = []
    try:
        with open(file_name[i].split('\n')[0]) as f1:
            data_lensed = json.load(f1)
        
        parameters = list(data_lensed['posterior']['content'].keys())
    except FileNotFoundError:
        print(file_name[i].split('\n')[0])
        match_filter_SNR.append(0)
        max_match_det_set.append('H1')
        continue
    except KeyError:
        print(file_name[i].split('\n')[0])
        match_filter_SNR.append(0)
        continue
        
    
    
    
    

    
    try:
        "读取cWB重构结果"
        cWB_path = './cWB_micro/data/' + str(int(inject_time_set[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/H1_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
    except FileNotFoundError:
        print(cWB_path)
        match_filter_SNR.append(0)
        max_match_det_set.append('H1')
        continue
    
    for reference_det_tmp in ["H1", "L1", "V1"]:
        with open(file_name[i].split('\n')[0]) as f1:
            data_lensed = json.load(f1)
        
        parameters = list(data_lensed['posterior']['content'].keys())
        
        "读取cWB重构结果"
        cWB_path = './cWB_micro/data/' + str(int(inject_time_set[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/' + reference_det_tmp + '_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
        
        h1t_rec_t = []
        h1t_rec = []
        for k in range(len(data)):
            h1t_rec_t.append(float(data[k].split()[0]))
            h1t_rec.append(float(data[k].split()[1]))
        

    
        start_index = 0
        end_index = 1000
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
        SNR_tmp = np.sqrt(data_lensed['posterior']['content']['H1_matched_filter_snr'][j]['real']**2 + data_lensed['posterior']['content']['L1_matched_filter_snr'][j]['real']**2 + + data_lensed['posterior']['content']['V1_matched_filter_snr'][j]['real']**2)
    
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
        #先找到共同的时间切片
        max_t = np.max([h1t_rec_t[0], signal_h1.sample_times[0]])
        min_t = np.min([h1t_rec_t[-1], signal_h1.sample_times[-1]])
        index_cwb = np.where((h1t_rec_t>max_t)&(h1t_rec_t<min_t))
        index_bilby = np.where((signal_h1.sample_times>max_t)&(signal_h1.sample_times<min_t))
        #cwb
        h1t_rec_t = np.array(h1t_rec_t)[index_cwb]
        h1t_rec = np.array(h1t_rec)[index_cwb]
        #bilby
        signal_h1_sample_times = signal_h1.sample_times[index_bilby]
        signal_h1 = signal_h1[index_bilby]
        
        #cwb
        h1f_rec = scipy.fft.fft(h1t_rec) * np.diff(h1t_rec_t)[0]
        h1f_rec_f = scipy.fft.fftfreq(len(h1t_rec), np.diff(h1t_rec_t)[0])
        func_rec = interp1d(h1f_rec_f, h1f_rec)
        #bilby
        
        tilde_signal_h1 = scipy.fft.fft(signal_h1) * np.diff(signal_h1_sample_times)[0]

        freq_h1 = scipy.fft.fftfreq(len(signal_h1), np.diff(signal_h1_sample_times)[0])



        Bilby_func_h1 = interp1d(freq_h1, tilde_signal_h1)
        
        #频率范围
        index_yes = np.where((h1f_rec_f>f_lower)&(h1f_rec_f<freq_max[i]))

        freq_match_rude = h1f_rec_f[index_yes]
        freq_match = np.arange(freq_match_rude[0], freq_match_rude[-1], 10**(-3))
        
        
        match_tmp_h1 = Match((func_rec(freq_match))/func_psd(freq_match)**(1/2), (Bilby_func_h1(freq_match))/func_psd(freq_match)**(1/2))
        if match_tmp_h1 > max_match:
            max_match = match_tmp_h1
            max_match_det = reference_det_tmp
            func_rec_max = func_rec
            Bilby_func_max = Bilby_func_h1
        
             
            
            
    match_micro.append(max_match)
    max_match_det_set.append(max_match_det)
    plt.figure(i)
    plt.title('reference_det = ' + max_match_det + ' SNR = ' + str(SNR_tmp) + 'match = ' + str(round(max_match,3)))
    plt.semilogy(freq_match, abs(func_rec_max(freq_match)))
    plt.semilogy(freq_match, abs(Bilby_func_max(freq_match)), '--')
    plt.plot([freq_max[i], freq_max[i]], [10**(-30), 10**(-22)])
    plt.xlim(10, 1024)
    plt.ylim(10**(-30),10**(-22))
    plt.grid()
    plt.savefig("./IntermediaPlot_match_lensed_complex/abs_" + file_name[i].split('/')[3][6::] + "_" + str(inject_time_set[i]) + ".png", dpi=450)
    plt.close()
    plt.figure(i)
    plt.title('reference_det = ' + max_match_det + ' SNR = ' + str(SNR_tmp) + 'match = ' + str(round(max_match,3)))
    plt.semilogy(freq_match, np.unwrap(np.angle(func_rec_max(freq_match))))
    plt.semilogy(freq_match, np.unwrap(np.angle(Bilby_func_max(freq_match))), '--')
    # plt.plot([freq_max[i], freq_max[i]], [10**(-30), 10**(-22)])
    # plt.xlim(10, 1024)
    # plt.ylim(10**(-30),10**(-22))
    plt.grid()
    plt.savefig("./IntermediaPlot_match_lensed_complex/phase_" + file_name[i].split('/')[3][6::] + "_" + str(inject_time_set[i]) + ".png", dpi=450)
    plt.close()
    
    match_filter_SNR.append(SNR_tmp)
    
    
    
    # np.savetxt('./match_result_Total/match_micro_'+file_name[i].split('/')[3][6::] + '.csv', match_micro, delimiter=',')
    np.savetxt('./Bilby_rec_micro_Total_complex/rec_'+file_name[i].split('/')[3][6::] + '.csv', tilde_signal_h1, delimiter=',')
    np.savetxt('./Bilby_rec_micro_Total_complex/freq_'+file_name[i].split('/')[3][6::] + '.csv', freq_h1, delimiter=',')

np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_complex/match_filter_SNR.csv', match_filter_SNR, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_complex/max_match_det_lensed.csv', max_match_det_set, delimiter=',', fmt='%s')


    
