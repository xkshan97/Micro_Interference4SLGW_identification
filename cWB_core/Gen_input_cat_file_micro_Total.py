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


sample_RATE = '4KHZ'

G = 6.67 * 10**(-11)
c = 3 * 10**(8)
M_sun = 2 * 10**(30)
apx = 'IMRPhenomPv2'
f_lower = 10
duration = 4096
sample_rate = 4096
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

# #only macro
SNR_network = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')

H1_snr = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/H1_snr.csv', delimiter=',')
L1_snr = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/L1_snr.csv', delimiter=',')
V1_snr = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/V1_snr.csv', delimiter=',')
# Index_SNR_gtr_12 = np.where(SNR_network>12)[0]
# SNR_select_index = Index_SNR_gtr_12[0:50]
# SNR_selected = np.array(SNR_network)[SNR_select_index]
# 


#挑选出来要用来被bilby和cWB估计的事例。 
match_total = np.loadtxt('./Theory_match_total/match.csv', delimiter=',')
Mz_total = np.loadtxt('./Theory_match_total/chirp_mass.csv', delimiter=',')
SNR_Lens_gtr_12_sum = np.loadtxt('./Theory_match_total/SNR_total.csv', delimiter=',')
selected_index1 = np.where((match_total<0.995)&(SNR_Lens_gtr_12_sum>50)&(SNR_Lens_gtr_12_sum<150)&(Mz_total>20))
selected_index2 = np.where((match_total<0.999)&(SNR_Lens_gtr_12_sum>150)&(Mz_total>20)) 
selected_index = np.append(selected_index1, selected_index2)
# selected_index = selected_index1

match_selected = match_total[selected_index]
SNR_selected_tmp = SNR_Lens_gtr_12_sum[selected_index]
cummulative_num = []
for i in range(len(imagenum)):
    cummulative_num.append(np.sum(imagenum[0:i+1]))
cummulative_num = np.array(cummulative_num)
IndexInAccept = []
IndexInTotalSNR = []
for i in range(len(SNR_selected_tmp)):
    index_tmp_snr = np.where(SNR_network == SNR_selected_tmp[i])[0][0]
    IndexInTotalSNR.append(index_tmp_snr) #这里的Total SNR是SNR_network里的
    index_tmp_accept = np.where(cummulative_num > index_tmp_snr)[0][0]
    if index_tmp_accept not in IndexInAccept:
        IndexInAccept.append(index_tmp_accept)



a0 = 2.9740 * 10**(-1)
b0 = 4.4810 * 10**(-2)
c0 = 9.5560 * 10**(-2)
a1 = 5.9411 * 10**(-1)
b1 = 8.9794 * 10**(-2)
c1 = 1.9111 * 10**(-1)
def merger_freq(m1, m2, z):
    eta = m2 / m1
    M = (1 + z) * (m1 + m2)
    return c**3*(a0 * eta ** 2 + b0 * eta + c0)/np.pi / M / G / M_sun
def ringdown_freq(m1, m2, z):
    eta = m2 / m1
    M = (1 + z) * (m1 + m2)
    return c**3*(a1 * eta ** 2 + b1 * eta + c1)/np.pi / M / G / M_sun
   

SNR_selected = []
inject_time_set = []
excel_name_match = [] #用来手动挑选运行哪个。
match = []
Mz = []
mass_ratio = []
merger_freq_set = []
ringdown_freq_set = []
reference_det = []
#把上面挑出来的事例，保存成cWB能用的文件名
# freq_max = [120, 160, 210, 90, 240, 160, 140, 170, 160, 190, 580, 110, 80, 50, 540, 110, 490, 390, 290, 170, 310, 50, 590, 240, 140, 410, 160, 60, 120, 360, 290, 160, 160, 70]
freq_index = 0
with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/with_micro.txt', 'w') as f0:
    for i_ in IndexInAccept:
        i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
        imag_num = int(imagenum[i_])
        save_index = 0
        IndexSumimage = int(np.sum(imagenum[0:i_]))
        Index_inject = np.sum(SNR_network[0:IndexSumimage] >12)
        for j_ in range(imag_num):
            if SNR_network[IndexSumimage] > 12:
                if SNR_network[IndexSumimage] in SNR_selected_tmp:
                    SNR_selected.append(SNR_network[IndexSumimage])
                    if H1_snr[IndexSumimage] >= L1_snr[IndexSumimage] and H1_snr[IndexSumimage] >= V1_snr[IndexSumimage]:
                        reference_det.append('H1')
                    elif L1_snr[IndexSumimage] > H1_snr[IndexSumimage] and L1_snr[IndexSumimage] >= V1_snr[IndexSumimage]:
                        reference_det.append('L1')
                    elif V1_snr[IndexSumimage] > H1_snr[IndexSumimage] and V1_snr[IndexSumimage] > L1_snr[IndexSumimage]:
                        reference_det.append('V1')
                    
                    td_second = timedelay[IndexSumimage] * 24 * 3600
                    SampleParam9_td = SampleParam[9][i] + td_second
                    
                    sqrt_mag_abs = np.sqrt(np.abs(magnification[IndexSumimage]))
                    
                    if 1 - kappa[IndexSumimage] - gamma[IndexSumimage] > 0 and 1 - kappa[IndexSumimage] + gamma[IndexSumimage] > 0:    
                        """================================="""
                        Lens_Type = "Minimum"
                    elif (1 - kappa[IndexSumimage] - gamma[IndexSumimage]) * (1 - kappa[IndexSumimage] + gamma[IndexSumimage]) < 0:
                        Lens_Type = "Saddle"
                    
                    
                    print(Lens_Type)
                    
                    
                    Ff_inter = np.loadtxt("../Paper4_MicroLens_Modify/Ffabs" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) +  ".csv", delimiter=',')
                    ThetaF_inter = np.loadtxt("../Paper4_MicroLens_Modify/ThetaF" + Lens_Type + "/Total" + str(i) + "_" + str(save_index) +  ".csv", delimiter=',')
                    
                    Ff_inter = Ff_inter / sqrt_mag_abs
                    
                    
                    #input
                    with open('./input_micro_Total/H1_' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.in', 'w') as f:
                        f.write('/home/cWB_docker/CE_Micro/catalog/GWTC-1-confident/GW150914/FRAMES/' + Lens_Type + '/H1/H-H1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096'+ '_' + str(i) + '_' + str(save_index) + '.gwf')
                    with open('./input_micro_Total/L1_' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.in', 'w') as f:
                        f.write('/home/cWB_docker/CE_Micro/catalog/GWTC-1-confident/GW150914/FRAMES/' + Lens_Type + '/L1/L-L1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096'+ '_' + str(i) + '_' + str(save_index) + '.gwf')
                    with open('./input_micro_Total/V1_' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.in', 'w') as f:
                        f.write('/home/cWB_docker/CE_Micro/catalog/GWTC-1-confident/GW150914/FRAMES/' + Lens_Type + '/V1/V-V1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096'+ '_' + str(i) + '_' + str(save_index) + '.gwf')
                    
                    # #copy file to FRAMES
                    shutil.copy('../Paper4_CE_Modify/Total_Sim_GW_Data_Lensed_' + Lens_Type + '_with_micro' + '/H1/H-H1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf', './FRAMES_micro_Total/' + Lens_Type + '/H1/H-H1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf')
                    shutil.copy('../Paper4_CE_Modify/Total_Sim_GW_Data_Lensed_' + Lens_Type + '_with_micro' + '/L1/L-L1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf', './FRAMES_micro_Total/' + Lens_Type + '/L1/L-L1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf')
                    shutil.copy('../Paper4_CE_Modify/Total_Sim_GW_Data_Lensed_' + Lens_Type + '_with_micro' + '/V1/V-V1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf', './FRAMES_micro_Total/' + Lens_Type + '/V1/V-V1_GWOSC_4KHZ_R1-' + str(int(Total_inject_time[Index_inject])) + '-4096_' + str(i) + '_' + str(save_index) + '.gwf')
                    
                    #DQ
                    start_gps = int(Total_inject_time[Index_inject] + 1000)
                    end_gps = int(start_gps + 3000)
                    with open('./DQ_micro_Total/H1_' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.txt', 'w') as f:
                        f.write("%s %s\n"%(start_gps, end_gps))
                    with open('./DQ_micro_Total/L1_' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.txt', 'w') as f:
                        f.write("%s %s\n"%(start_gps, end_gps))
                    with open('./DQ_micro_Total/V1_' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.txt', 'w') as f:
                        f.write("%s %s\n"%(start_gps, end_gps))
                    
                    excel_name_match.append(str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index))
                    inject_time_set.append(int(Total_inject_time[Index_inject]))
                    
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
                            
                    
                    merger_index = np.where(hp_orign.sample_times <= 0)[0][-1]
                    if(merger_index == np.where(hc_orign.sample_times <= 0)[0][-1]):
                        pass
                    else:
                        print("Error !!!!!!!!!!!!!!, GW index i = ", i)
                    
                    Mz_tmp = (SampleParam[0][i] + 1) * (SampleParam[1][i] * SampleParam[2][i]) ** (3/5) / (SampleParam[1][i] + SampleParam[2][i]) ** (1/5)
                    mass_ratio_tmp = SampleParam[2][i] / SampleParam[1][i]
                    Mz.append(Mz_tmp)
                    mass_ratio.append(mass_ratio_tmp)
                    
                    # hp = hp.time_slice(hp.end_time - hp.delta_t * hp.sample_rate * 6,hp.end_time)
                    # hc = hc.time_slice(hc.end_time - hc.delta_t * hc.sample_rate * 6,hc.end_time)
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
                            
                    
                    #透镜波形
                    tilde_signal_h1_lensed = tilde_signal_h1 * Ff_inter * np.exp(complex(0,1) * ThetaF_inter)
                    tilde_signal_h1_lensed_inter = interp1d(freq_h1, np.abs(tilde_signal_h1_lensed))
                    tilde_signal_h1_inter = interp1d(freq_h1, np.abs(tilde_signal_h1))
                    signal_h1_lensed = np.real(scipy.fft.ifft(tilde_signal_h1_lensed * sqrt_mag_abs) / signal_h1.delta_t)
                    np.savetxt('./Micro_signal_Total/freq_'+str(int(i))+'_'+str(int(save_index))+'.csv', freq_h1, delimiter=',') 
                    np.savetxt('./Micro_signal_Total/lens_signal_'+str(int(i))+'_'+str(int(save_index))+'.csv', abs(tilde_signal_h1_lensed * sqrt_mag_abs), delimiter=',') 
                    np.savetxt('./Micro_signal_Total/time_'+str(int(i))+'_'+str(int(save_index))+'.csv', signal_h1.sample_times, delimiter=',') 
                    np.savetxt('./Micro_signal_Total/time_lens_signal_'+str(int(i))+'_'+str(int(save_index))+'.csv', signal_h1_lensed, delimiter=',') 
                    
                     
                    index_yes = np.where((freq_h1>10)&(freq_h1<1024))
                    freq_index += 1
                    freq_match_rude = freq_h1[index_yes]
                    freq_match = np.arange(freq_match_rude[0], freq_match_rude[-1], 10**(-4))
                    
                    match_tmp = np.sum(np.abs(tilde_signal_h1_inter(freq_match)) * np.abs(tilde_signal_h1_lensed_inter(freq_match)) / psd_inter(freq_match)) \
                        / np.sqrt(np.sum(np.abs(tilde_signal_h1_inter(freq_match))**2 / psd_inter(freq_match)) * np.sum(np.abs(tilde_signal_h1_lensed_inter(freq_match))**2 / psd_inter(freq_match)))
                    match.append(match_tmp)
                    
                    plt.figure(i+save_index)
                    plt.semilogy(freq_h1[index_yes], np.abs(tilde_signal_h1[index_yes]))
                    plt.semilogy(freq_h1[index_yes], np.abs(tilde_signal_h1_lensed[index_yes]), '--')
                    plt.title(str(i)+'_'+str(round(match_tmp,3))+'_'+str(round(SNR_selected[-1])))
                    # plt.plot([merger_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]), merger_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i])], [10**(-30), 10**(-22)], 'r', label='Merger freq')
                    # plt.plot([ringdown_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]), ringdown_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i])], [10**(-30), 10**(-22)], 'k', label='ringdown freq')
                    plt.plot([128,128], [10**(-30), 10**(-22)], 'r')
                    plt.plot([256,256], [10**(-30), 10**(-22)], 'b')
                    plt.plot([1024, 1024], [10**(-30), 10**(-22)], 'y')
                    plt.plot([512, 512], [10**(-30), 10**(-22)], 'k', '--')
                    plt.grid()
                    plt.ylim(10**(-27), 10**(-23))
                    plt.xticks(np.arange(10,1200,100))
                    plt.savefig('./IntermediaPlot/Micro_selected/' + str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '.png', dpi=450)
                    plt.close()
                    # plt.savefig('./FreqDomain_figure_micro_fstar/'+str(int(Total_inject_time[Index_inject])) + '_' + str(i) + '_' + str(save_index) + '_' + f_star+'.png', dpi=450)
                    # plt.close()
                    
                    
                    # plt.figure(i)
                    # plt.semilogx(freq_h1, Ff_inter)
                    # plt.title(str(i)+'_'+str(round(match_tmp,3)))
                    
                    f0.write('../Paper4_CE_Modify/outdir_Lensed_'+Lens_Type+'_with_micro/outdir'+str(int(i))+'_'+str(int(save_index))+'/H1_L1_V1_result.json'+'\n')
                    
                    
                    
                    merger_freq_set.append(merger_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]))
                    ringdown_freq_set.append(ringdown_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]))
                    
                save_index += 1
                Index_inject += 1
        
                
            IndexSumimage += 1
# plt.plot(match)
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/match.csv', match, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/excel_name_match.csv', excel_name_match, delimiter=',', fmt='%s')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/SNR_selected.csv', SNR_selected, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/IndexInAccept.csv', IndexInAccept, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/inject_time_set.csv', inject_time_set, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/chirp_mass.csv', Mz, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/mass_ratio.csv', mass_ratio, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/merger_freq.csv', merger_freq_set, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/ringdown_freq.csv', ringdown_freq_set, delimiter=',')
np.savetxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/reference_det.csv', reference_det, delimiter=',', fmt='%s')