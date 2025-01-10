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

SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')
H1_snr = np.loadtxt('../Paper4_CE_Modify/Sim_GW_Data_big_sample/H1/H1_snr.csv', delimiter=',')
L1_snr = np.loadtxt('../Paper4_CE_Modify/Sim_GW_Data_big_sample/L1/L1_snr.csv', delimiter=',')
V1_snr = np.loadtxt('../Paper4_CE_Modify/Sim_GW_Data_big_sample/V1/V1_snr.csv', delimiter=',')

SNR_network = np.sqrt(H1_snr**2 + L1_snr**2 + V1_snr**2)
"""match function"""
def Match(strain1, strain2):
    cor_match = 4 * np.sum(strain1 * np.conjugate(strain2))
    match1 = 4 * np.sum(strain1 * np.conjugate(strain1))
    match2 = 4 * np.sum(strain2 * np.conjugate(strain2))
    return np.abs(cor_match / np.sqrt(match1 * match2))

f_lower = 40
duration = 4096 #"""注意这里改成了4096，paper4_copy中是64."""
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_ligo = pycbc.psd.CosmicExplorerP1600143(fsamples, df, 20)
func_psd = interp1d(psd_ligo.sample_frequencies, psd_ligo.data)

inject_time_set = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/Data_start_time_set_unlens.csv', delimiter=',')
with open('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/reference_det.csv', 'r') as f1:
    reference_det = f1.readlines()
with open('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/unlens.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据
freq_max = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/merger_freq_unlens.csv', delimiter=',', encoding='UTF-8-sig')
# freq_max = np.loadtxt('./cWB_unlens/Unlens_freq.csv', delimiter=',', encoding='UTF-8-sig')
freq_max[freq_max<128] = 128
freq_max[freq_max>512] = 512
name_event = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/excel_name_match_unlens.csv', dtype=str)

#上面是从图上用眼睛读出来的数据
match_filter_SNR = []
max_match_det_set = []
max_match_set = []
noise_mean = []
noise_std = []
for i in tqdm(range(len(inject_time_set))):
    
    max_match = 0
    # SNR_tmp = []
    match_unlens = []
    reference_det_tmp = reference_det[i][0:2]
    
    try:
        "读取cWB重构结果"
        path_tmp = '../Paper5_PycWB/Unlens_big_sample/' + name_event[i] +'/trigger'
        file_tmp = os.listdir(path_tmp)
        PycWB_res = np.loadtxt(path_tmp + '/' + file_tmp[0] + '/reconstructed_waveform_' + reference_det_tmp + '.txt')
        
    except Exception as e:
        print(path_tmp)
        match_filter_SNR.append(0)
        max_match_set.append(0)
        max_match_det_set.append('H1')
        noise_mean.append(10**5)
        noise_std.append(10**5)
        continue
    
    name_i = int(name_event[i].split('_')[-1])
    for reference_det_tmp in ["H1", "L1", "V1"]:
        
        f_lower = 40
        "读取cWB重构结果"
        path_tmp = '../Paper5_PycWB/Unlens_big_sample/' + name_event[i] +'/trigger'
        file_tmp = os.listdir(path_tmp)
        
        # SNR_tmp = SNR_network[name_i]
        
        
        
        luminosity_dis = Planck15.luminosity_distance(SampleParam[0][name_i]).value
        
        hp_orign, hc_orign = get_td_waveform(approximant=apx,
                    mass1=SampleParam[1][name_i] * (1 + SampleParam[0][name_i]), mass2=SampleParam[2][name_i] * (1 + SampleParam[0][name_i]),
                    spin1z = SampleParam[3][name_i], spin2z = SampleParam[4][name_i],
                    inclination = SampleParam[5][name_i],
                    distance = luminosity_dis,
                    delta_t=1.0/sample_rate, f_lower=10, f_ref = 50, coa_phase=0)

        hp_20, hc_20 = get_td_waveform(approximant=apx,
                            mass1=SampleParam[1][name_i] * (1 + SampleParam[0][name_i]), mass2=SampleParam[2][name_i] * (1 + SampleParam[0][name_i]),
                            spin1z = SampleParam[3][name_i], spin2z = SampleParam[4][name_i],
                            inclination = SampleParam[5][name_i],
                            distance = luminosity_dis,
                            delta_t=1.0/sample_rate, f_lower=20, f_ref = 50, coa_phase=0)
        

            
        #把时间调成0试试，这样做意味着警报时间和开始时间对齐了

        # merger_index = np.where(hp.sample_times <= 0)[0][-1]
        hp_orign.start_time = 0
        hc_orign.start_time = 0
        
        

        det_h1 = Detector(reference_det_tmp)
        


        #不同探测器的波形
        hp_h1, hc_h1 = hp_orign.copy(), hc_orign.copy()
    
        
        #得到不同探测器的时间
        start_time_signal = int(SampleParam[9][name_i]) - (hp_orign.duration - hp_20.duration) #保证30Hz的开始时间是在start time上
        start_time = int(SampleParam[9][name_i])
        
        declination = SampleParam[8][name_i]
        right_ascension = SampleParam[7][name_i]
        polarization = SampleParam[6][name_i]
        
        #校准不同探测的波形时间
        hp_h1.start_time += start_time_signal
        hc_h1.start_time += start_time_signal
        
        
        
        # We get back the fp and fc antenna pattern weights.
        fp_h1, fc_h1 = det_h1.antenna_pattern(right_ascension, declination, polarization, start_time)
        
        # 从地心到探测器的时间
        dt_h1 = det_h1.time_delay_from_earth_center(right_ascension, declination, start_time)
        
    
        
        #校准不同探测的波形时间
        hp_h1.start_time += dt_h1 
        hc_h1.start_time += dt_h1
        

        ## Apply the factors to get the detector frame strain
        signal_h1_inj = fp_h1 * hp_h1 + fc_h1 * hc_h1
        
        
        tilde_signal_h1 = scipy.fft.fft(signal_h1_inj.data) * signal_h1_inj.delta_t

        freq_h1 = scipy.fft.fftfreq(len(signal_h1_inj), signal_h1_inj.delta_t)
        

        Bilby_func_h1 = interp1d(freq_h1, np.abs(tilde_signal_h1))
        
        for file_i in range(len(file_tmp)):
            try:
                PycWB_res = np.loadtxt(path_tmp + '/' + file_tmp[file_i] + '/reconstructed_waveform_' + reference_det_tmp + '.txt')
                
                
                with open(path_tmp + '/' + file_tmp[file_i] + '/event.json') as event_josn:
                    read_res = event_josn.readlines()
                
                SNR_tmp = np.sqrt(np.sum([float(snr_) for snr_ in read_res[0].split("snr")[1].split('[')[1].split(']')[0].split(',')]))
                
                h1t_rec_t = PycWB_res[:,0]
                h1t_rec = PycWB_res[:,1]
                
                #cwb
                h1f_rec = scipy.fft.fft(h1t_rec) * np.diff(h1t_rec_t)[0]
                h1f_rec_f = scipy.fft.fftfreq(len(h1t_rec), np.diff(h1t_rec_t)[0])
                func_rec = interp1d(h1f_rec_f, np.abs(h1f_rec))


            
            
            #频率范围
            
                freq_max_bilby = freq_h1[freq_h1>f_lower][np.where(np.abs(tilde_signal_h1[freq_h1>f_lower])/np.abs(tilde_signal_h1[freq_h1>f_lower])[0] < 0.3)[0][0]]
                freq_max_cwb = h1f_rec_f[h1f_rec_f>f_lower][np.where(np.abs(h1f_rec[h1f_rec_f>f_lower])/np.abs(h1f_rec[h1f_rec_f>f_lower])[0] < 0.3)[0][0]]
                freq_max_ = np.min([freq_max_bilby, freq_max_cwb, 256])
                freq_max_ = np.max([freq_max_, 64])
            
                freq_min_cwb = h1f_rec_f[np.where(np.abs(h1f_rec) == np.max(np.abs(h1f_rec)))][0]
                if freq_min_cwb > freq_max_:
                    f_lower = freq_max_ - 10
                else:
                    f_lower = freq_min_cwb
            
                index_yes = np.where((h1f_rec_f>f_lower)&(h1f_rec_f<freq_max_))

                freq_match_rude = h1f_rec_f[index_yes]
                freq_match = np.arange(freq_match_rude[0], freq_match_rude[-1], 10**(-1))
                
            
                match_tmp_h1 = Match(np.abs(func_rec(freq_match))/func_psd(freq_match)**(1/2), np.abs(Bilby_func_h1(freq_match))/func_psd(freq_match)**(1/2))
                if match_tmp_h1 > max_match:
                    max_match = match_tmp_h1
                    max_match_det = reference_det_tmp
                    func_rec_max = func_rec
                    Bilby_func_max = Bilby_func_h1
                    max_snr = SNR_tmp
                    
            except Exception as e:
                print(e)
        
    match_unlens.append(max_match)
    max_match_set.append(max_match)
    max_match_det_set.append(max_match_det)
    noise_mean.append(np.mean(np.abs(func_rec(freq_match)) - np.abs(Bilby_func_h1(freq_match))))
    noise_std.append(np.std(np.abs(func_rec(freq_match)) - np.abs(Bilby_func_h1(freq_match))))
    if max_snr > 200:
        plt.figure(i)
        plt.title('reference_det = ' + max_match_det + ' SNR = ' + str(SNR_tmp) + 'match = ' + str(round(max_match,3)))
        plt.semilogy(freq_match, abs(func_rec_max(freq_match)))
        plt.semilogy(freq_match, abs(Bilby_func_max(freq_match)), '--')
        plt.plot([freq_max[i], freq_max[i]], [10**(-30), 10**(-22)])
        plt.xlim(10, 1024)
        # plt.ylim(10**(-25),10**(-22))
        plt.grid()
        plt.savefig("./IntermediaPlot_match_unlens_big_sample/abs_" + file_name[i].split('/')[3][6::] + "_" + str(inject_time_set[i]) + ".png", dpi=450)
        plt.close()
        plt.figure(i)
        plt.title('reference_det = ' + max_match_det + ' SNR = ' + str(SNR_tmp) + 'match = ' + str(round(max_match,3)))
        plt.semilogy(freq_match, np.unwrap(np.angle(func_rec_max(freq_match))))
        plt.semilogy(freq_match, np.unwrap(np.angle(Bilby_func_max(freq_match))), '--')
        # plt.plot([freq_max[i], freq_max[i]], [10**(-30), 10**(-22)])
        # plt.xlim(10, 1024)
        # plt.ylim(10**(-30),10**(-22))
        plt.grid()
        plt.savefig("./IntermediaPlot_match_unlens_big_sample/phase_" + file_name[i].split('/')[3][6::] + "_" + str(inject_time_set[i]) + ".png", dpi=450)
        plt.close()
    match_filter_SNR.append(max_snr)
    
    
    # np.savetxt('./match_result/match_unlens_'+file_name[i].split('/')[3][6::]+'.csv', match_unlens, delimiter=',')
  
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/match_filter_SNR_unlens.csv', match_filter_SNR, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/max_match_det_unlens.csv', max_match_det_set, delimiter=',', fmt='%s')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/max_match_set_unlens.csv', max_match_set, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/noise_mean.csv', noise_mean, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/noise_std.csv', noise_std, delimiter=',')

    




match_filter_SNR = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/match_filter_SNR_unlens.csv', delimiter=',')
max_match_set = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/max_match_set_unlens.csv', delimiter=',')
chirp_mass = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/chirp_mass_unlens.csv')
index_mass_gtr_20 = np.where(chirp_mass > 20)

noise_mean = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/noise_mean.csv', delimiter=',')
noise_std = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/noise_std.csv', delimiter=',')
noise_mean[noise_mean==10**5] = 0
noise_std[noise_std==10**5] = 0
plt.scatter(match_filter_SNR, np.abs(noise_mean), color='r')
plt.scatter(match_filter_SNR, np.abs(noise_std), color='b')
plt.semilogy()
# plt.ylim(0.94, 1)
# plt.xlim(40, 500)
plt.savefig("./IntermediaPlot_match_unlens_big_sample/noise_level.png", dpi=450)
plt.close()

plt.scatter(match_filter_SNR[index_mass_gtr_20], max_match_set[index_mass_gtr_20])
plt.ylim(0.94, 1)
plt.xlim(40, 500)
plt.savefig("./IntermediaPlot_match_unlens_big_sample/SNR_VS_match.png", dpi=450)