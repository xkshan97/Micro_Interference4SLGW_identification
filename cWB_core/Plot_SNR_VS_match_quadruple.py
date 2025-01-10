'''
计算所有微透镜事例（4像）的理论match，信噪比，以及chirp mass
'''
import numpy as np
import matplotlib.pyplot as plt
import shutil
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
from astropy.cosmology import Planck15
from scipy.interpolate import interp1d
import scipy
import pycbc.psd
import os
from tqdm import tqdm
import corner


Gravity_con = 6.67 * 10**(-11)
light_speed = 3 * 10**(8)
sun_mass = 2 * 10**(30)
apx = 'IMRPhenomPv2'
f_lower = 10
duration = 4096
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)
psd_inter = interp1d(psd_l1_h1.sample_frequencies, psd_l1_h1)


SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')
AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')


magnification = np.loadtxt('../Paper4_CE_Modify/SampleResult/magnification.csv', delimiter=',')
timedelay = np.loadtxt('../Paper4_CE_Modify/SampleResult/timedelayDays.csv', delimiter=',')
kappa = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('../Paper4_CE_Modify/SampleResult/gamma.csv', delimiter=',')
imagenum = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
lens_z_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')
Total_inject_time = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult_Micro/Total_inject_time.csv', delimiter=',')


#四像的信息
four_image = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/four_image_index.csv', delimiter=',') 
four_kappa = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/four_image_kappa.csv', delimiter=',')
four_gamma = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/four_image_gamma.csv', delimiter=',')
four_mag = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/four_image_mag.csv', delimiter=',') 
four_td = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/four_image_td.csv', delimiter=',')


#only macro
SNR_network = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')
# imagenum = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
# kappa = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa.csv', delimiter=',')
# gamma = np.loadtxt('../Paper4_CE_Modify/SampleResult/gamma.csv', delimiter=',')
#micro SNR
# SNR_network_micro = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/micro_net_snr.csv',delimiter=',')




#有微透镜的
match = []
Mz = []
SNR_gtr12_total = []
redshift_sampled = []
#把上面挑出来的事例，保存成cWB能用的文件名
# freq_max = [120, 160, 210, 90, 240, 160, 140, 170, 160, 190, 580, 110, 80, 50, 540, 110, 490, 390, 290, 170, 310, 50, 590, 240, 140, 410, 160, 60, 120, 360, 290, 160, 160, 70]
for i_ in tqdm(range(len(AcceptLensIndex))):
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = int(imagenum[i_])
    save_index = 0
    IndexSumimage = int(np.sum(imagenum[0:i_]))
    for j_ in range(imag_num):
        
        if SNR_network[IndexSumimage] > 12 and i_ in four_image:
            
            SNR_gtr12_total.append(SNR_network[IndexSumimage])
            
            td_second = timedelay[IndexSumimage] * 24 * 3600
            SampleParam9_td = SampleParam[9][i] + td_second
            
            sqrt_mag_abs = np.sqrt(np.abs(magnification[IndexSumimage]))
            if 1 - kappa[IndexSumimage] - gamma[IndexSumimage] > 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] > 0:    
                """================================="""
                Lens_Type = "Minimum"
            elif (1 - kappa[IndexSumimage] - gamma[IndexSumimage]) * (1 - kappa[IndexSumimage] + gamma[IndexSumimage]) < 0:
                Lens_Type = "Saddle"
            
        
            Ff_inter = np.loadtxt("../Paper4_MicroLens_Modify/Ffabs" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) +  "_quadruple.csv", delimiter=',')
            ThetaF_inter = np.loadtxt("../Paper4_MicroLens_Modify/ThetaF" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) +  "_quadruple.csv", delimiter=',')
            Ff_inter = Ff_inter / sqrt_mag_abs
            #顺便计算一下带微透镜波形和不带微透镜波形的match程度
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
                    
            Mz_tmp = (SampleParam[0][i] + 1) * (SampleParam[1][i] * SampleParam[2][i]) ** (3/5) / (SampleParam[1][i] + SampleParam[2][i]) ** (1/5)
            Mz.append(Mz_tmp)
            redshift_sampled.append(SampleParam[0][i])
            #innermost stable circular orbit, ISCO
            innermost_sco = light_speed ** 3 / 6 ** (3/2) / np.pi / Gravity_con / (sun_mass * (SampleParam[1][i] + SampleParam[2][i])* (1 + SampleParam[0][i])) 
            
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

            # hp, hc = get_td_waveform(approximant=apx,
            #     mass1=36., mass2=29.,
            #     distance=2000., inclination=0.4, polarization=2.659,
            #     ra=1.375, dec=-1.2108,delta_t=1.0/sample_rate, f_lower=f_lower)

            det_h1 = Detector('H1')

            #不同探测器的波形
            hp_h1, hc_h1 = hp_orign.copy(), hc_orign.copy()

            #得到不同探测器的时间
            start_time_signal = int(SampleParam9_td) - (hp_orign.duration - hp_20.duration) #保证20Hz的开始时间是在start time上
            start_time = int(SampleParam9_td)
            
            declination = SampleParam[8][i]
            right_ascension = SampleParam[7][i]
            polarization = SampleParam[6][i]

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
            signal_h1 = fp_h1 * hp_h1 + fc_h1 * hc_h1
            #得到频域的信号
            tilde_signal_h1 = scipy.fft.fft(signal_h1.data) * signal_h1.delta_t
            
            freq_h1 = scipy.fft.fftfreq(len(signal_h1), signal_h1.delta_t)
            
            # plt.figure(i_)
            # plt.semilogy(freq_h1, np.abs(tilde_signal_h1))
            # plt.title(i_)
            # plt.xlim(40,1000)
            # plt.grid()
            # plt.xticks(np.arange(40,1000,100))
            #透镜波形
            tilde_signal_h1_lensed = tilde_signal_h1 * Ff_inter * np.exp(complex(0,1) * ThetaF_inter)
            tilde_signal_h1_lensed_inter = interp1d(freq_h1, np.abs(tilde_signal_h1_lensed))
            tilde_signal_h1_inter = interp1d(freq_h1, np.abs(tilde_signal_h1))
            
            index_yes = np.where((freq_h1>10)&(freq_h1<1024))
            freq_match_rude = freq_h1[index_yes]
            freq_match = np.arange(freq_match_rude[0], freq_match_rude[-1], 10**(-4))
            
            match_tmp = np.sum(np.abs(tilde_signal_h1_inter(freq_match)) * np.abs(tilde_signal_h1_lensed_inter(freq_match)) / psd_inter(freq_match)) \
                / np.sqrt(np.sum(np.abs(tilde_signal_h1_inter(freq_match))**2 / psd_inter(freq_match)) * np.sum(np.abs(tilde_signal_h1_lensed_inter(freq_match))**2 / psd_inter(freq_match)))
            match.append(match_tmp)
            
            save_index += 1
        IndexSumimage += 1
    # plt.figure(i_)
    # plt.loglog(freq_h1[index_yes], np.abs(tilde_signal_h1[index_yes]))
    # plt.loglog(freq_h1[index_yes], np.abs(tilde_signal_h1_lensed[index_yes]), '--')
    # plt.title(str(i_)+'_'+str(round(match_tmp,3)))
        

# plt.plot(match)
np.savetxt('./Theory_match_total_quadruple/match.csv', match, delimiter=',')
np.savetxt('./Theory_match_total_quadruple/chirp_mass.csv', Mz, delimiter=',')
np.savetxt('./Theory_match_total_quadruple/SNR_total.csv', SNR_gtr12_total, delimiter=',')
np.savetxt('./Theory_match_total_quadruple/redshift_sampled.csv', redshift_sampled, delimiter=',')




#图
Mz_cut = 20
match = np.loadtxt('./Theory_match_total_quadruple/match.csv', delimiter=',')
Mz = np.loadtxt('./Theory_match_total_quadruple/chirp_mass.csv', delimiter=',')
SNR_Lens_gtr_12_sum = np.loadtxt('./Theory_match_total_quadruple/SNR_total.csv', delimiter=',')

samples = []
for _ in range(len(SNR_Lens_gtr_12_sum)):
    
    # m1 = chirp_mass[i]*(1 + mass_ratio[i])**(1/5)/mass_ratio[i]**(3/5)
    # m2 = mass_ratio[i]*m1
    if Mz[_] > Mz_cut:
        samples.append([SNR_Lens_gtr_12_sum[_], match[_]]) 

samples = np.array(samples)
figure = corner.corner(samples, labels=["SNR", "match"],
        show_titles=True, title_kwargs={"fontsize": 12},smooth=False, title_fmt='.4f')
plt.savefig('./Plot_result_quadruple/corner_SNR_VS_match.png', dpi=450)
frac_1 = np.sum((match<0.995)&(SNR_Lens_gtr_12_sum>50)&(SNR_Lens_gtr_12_sum<200)&(Mz>Mz_cut)) / len(SNR_Lens_gtr_12_sum[np.where((Mz>Mz_cut))])
frac_2 = np.sum((match<0.9999)&(SNR_Lens_gtr_12_sum>150)&(Mz>Mz_cut)) / len(SNR_Lens_gtr_12_sum[np.where((Mz>Mz_cut))])
frac_sum = frac_1 + frac_2
print(frac_sum)

