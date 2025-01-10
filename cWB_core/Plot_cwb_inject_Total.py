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

inject_time_set = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/inject_time_set.csv', delimiter=',')
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/with_micro.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据

#上面是从图上用眼睛读出来的数据
match_filter_SNR = []
for i in tqdm(range(len(inject_time_set))):

    try:
        "读取cWB重构结果"
        cWB_path = './cWB_micro/data/' + str(int(inject_time_set[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/H1_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
    except FileNotFoundError:
        print(cWB_path)
        match_filter_SNR.append(0)
        continue
    
    h1t_rec_t = []
    h1t_rec = []
    for k in range(len(data)):
        h1t_rec_t.append(float(data[k].split()[0]))
        h1t_rec.append(float(data[k].split()[1]))
    h1f_rec = scipy.fft.fft(h1t_rec) * np.diff(h1t_rec_t)[0]
    h1f_rec_f = scipy.fft.fftfreq(len(h1t_rec), np.diff(h1t_rec_t)[0])
    func_rec = interp1d(h1f_rec_f, abs(h1f_rec))
    
    inject_freq = np.loadtxt('./Micro_signal_Total/freq_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    inject_signal = np.loadtxt('./Micro_signal_Total/lens_signal_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    inject_time = np.loadtxt('./Micro_signal_Total/time_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
    time_inject_signal = np.loadtxt('./Micro_signal_Total/time_lens_signal_' + file_name[i].split('/')[3][6::] + '.csv', delimiter=',')
     
    
    
    plt.semilogy(h1f_rec_f, np.abs(h1f_rec), label='cWB')
    plt.semilogy(inject_freq, inject_signal, '--',label='injected')
    plt.xlim(10, 200)
    plt.ylim(10**(-26), 10**(-23))
    plt.legend()
    plt.grid()
    plt.savefig('IntermediaPlot/Micro_cWB_inject/freq_' + file_name[i].split('/')[3][6::] + '.png', dpi=450)
    plt.close()
    
    
    plt.plot(h1t_rec_t, h1t_rec, label='cWB')
    plt.plot(inject_time, time_inject_signal, '--',label='injected')
    plt.xlim(inject_time[len(inject_time)*9//10], inject_time[-1])
    plt.legend()
    plt.grid()
    plt.savefig('IntermediaPlot/Micro_cWB_inject/time_' + file_name[i].split('/')[3][6::] + '.png', dpi=450)
    plt.close()
    
    
