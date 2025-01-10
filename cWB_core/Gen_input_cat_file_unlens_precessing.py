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

sample_RATE = '4KHZ'

G = 6.674 * 10**(-11)
c = 3 * 10**(8)
M_sun = 2 * 10**(30)
SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')
#所有采样点的chirp mass
Mz_total = (SampleParam[0] + 1) * (SampleParam[1] * SampleParam[2]) ** (3/5) / (SampleParam[1] + SampleParam[2]) ** (1/5) 
mass_ratio_total = SampleParam[2] / SampleParam[1]

H1_snr = np.loadtxt('../Paper4_CE_Modify/Sim_GW_Data/H1/H1_snr.csv', delimiter=',')
L1_snr = np.loadtxt('../Paper4_CE_Modify/Sim_GW_Data/L1/L1_snr.csv', delimiter=',')
V1_snr = np.loadtxt('../Paper4_CE_Modify/Sim_GW_Data/V1/V1_snr.csv', delimiter=',')


SNR_network = np.sqrt(H1_snr**2 + L1_snr**2 + V1_snr**2)
index_add = list(np.where(SNR_network>300)[0][0:10]) #SNR大于300
# Index_SNR_network_G_12 = np.append(np.arange(0,35), np.where((SNR_network >= 12)&(Mz_total>20))[0][23::])
Index_SNR_network_G_12 = np.where((SNR_network >= 12)&(SNR_network <= 300))[0]
Index100_Index_SNR_network_G_12 = np.append(Index_SNR_network_G_12[0:200], index_add)

SNR_selected_unlens = SNR_network[Index100_Index_SNR_network_G_12]
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/SNR_select_index_unlens.csv', Index100_Index_SNR_network_G_12, delimiter=',')
Mz = Mz_total[Index100_Index_SNR_network_G_12]
mass_ratio = mass_ratio_total[Index100_Index_SNR_network_G_12]

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

#unlens
excel_name_match_unlens = [] #保存文件，然后适合在excel中手动挑选。
#把上面挑出来的事例，保存成cWB能用的文件名
Data_start_time_set_unlens = []
merger_freq_set = []
ringdown_freq_set = []
reference_det = []
with open('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/unlens.txt', 'w') as f0:
    for index_i in tqdm(range(len(Index100_Index_SNR_network_G_12))):
        i = Index100_Index_SNR_network_G_12[index_i]
        if H1_snr[i] >= L1_snr[i] and H1_snr[i] >= V1_snr[i]:
            reference_det.append('H1')
        elif L1_snr[i] > H1_snr[i] and L1_snr[i] >= V1_snr[i]:
            reference_det.append('L1')
        elif V1_snr[i] > H1_snr[i] and V1_snr[i] > L1_snr[i]:
            reference_det.append('V1')
            
        
        
        #找到文件所对应的inject time
        file_name_list_tmp = os.listdir('../Paper4_CE_Modify/Sim_GW_Data_precessing/H1')
        for pp in range(len(file_name_list_tmp)):
            if '-4096_' + str(i) + '.gwf' in file_name_list_tmp[pp]:
                file_name_tmp = file_name_list_tmp[pp]
        Data_start_time = file_name_tmp.split('-')[2]
        
        with open('./input_unlens_precessing/H1_' + Data_start_time + '_' + str(i) + '.in', 'w') as f:
            f.write('/disk1/home/shanxk/work/Paper4_CE_Modify/Sim_GW_Data_precessing/H1/H-H1_GWOSC_' + sample_RATE + '_R1-' + Data_start_time + '-4096'+ '_' + str(i) + '.gwf')
        with open('./input_unlens_precessing/L1_' + Data_start_time + '_' + str(i) + '.in', 'w') as f:
            f.write('/disk1/home/shanxk/work/Paper4_CE_Modify/Sim_GW_Data_precessing/L1/L-L1_GWOSC_' + sample_RATE + '_R1-' + Data_start_time + '-4096'+ '_' + str(i) + '.gwf')
        with open('./input_unlens_precessing/V1_' + Data_start_time + '_' + str(i) + '.in', 'w') as f:
            f.write('/disk1/home/shanxk/work/Paper4_CE_Modify/Sim_GW_Data_precessing/V1/V-V1_GWOSC_' + sample_RATE + '_R1-' + Data_start_time + '-4096'+ '_' + str(i) + '.gwf')
        
        
            #画图看一下需要分析的频率范围
        luminosity_dis = Planck15.luminosity_distance(SampleParam[0][i]).value
        # NOTE: Inclination runs from 0 to pi, with poles at 0 and pi
        #       coa_phase runs from 0 to 2 pi.
        hp, hc = get_fd_waveform(approximant='IMRPhenomPv2',
                        mass1=SampleParam[1][i] * (1 + SampleParam[0][i]), mass2=SampleParam[2][i] * (1 + SampleParam[0][i]),
                        spin1z = SampleParam[3][i], spin2z = SampleParam[4][i],
                        inclination = SampleParam[5][i], polarization = SampleParam[6][i],
                        right_ascension = SampleParam[7][i], declination = SampleParam[8][i],
                        distance = luminosity_dis,
                        delta_f=1/4, f_lower=10)
        
        
        plt.semilogy(hp.sample_frequencies, abs(hp))
        # plt.plot([merger_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]), merger_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i])], [10**(-30), 10**(-22)], 'r', label='Merger freq')
        # plt.plot([ringdown_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]), ringdown_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i])], [10**(-30), 10**(-22)], 'k', label='ringdown freq')
        plt.plot([128,128], [10**(-30), 10**(-22)], 'r')
        plt.plot([256,256], [10**(-30), 10**(-22)], 'b')
        plt.plot([1024, 1024], [10**(-30), 10**(-22)], 'y')
        plt.plot([512, 512], [10**(-30), 10**(-22)], 'k', '--')
        plt.title('i = ' + str(i) + ' Mz = ' + str(round(Mz_total[i],2)) + ' q = ' + str(round(mass_ratio_total[i],2)))
        plt.grid()
        plt.xticks(np.arange(10,1200,100))
        plt.xlim(0,1200)
        plt.savefig('./IntermediaPlot/Unlens_selected_precessing/' + Data_start_time + '_' + str(i) + '.png', dpi=450)
        plt.close()
        # plt.savefig('./FreqDomain_figure/'+Data_start_time + '_' + str(i)+'.png', dpi=450)
        # plt.close()
        #DQ
        start_gps = int(int(Data_start_time) + 1000)
        end_gps = int(start_gps + 3000)
        with open('./DQ_unlens_precessing/H1_' + Data_start_time + '_' + str(i) + '.txt', 'w') as f:
            f.write("%s %s\n"%(start_gps, end_gps))
        with open('./DQ_unlens_precessing/L1_' + Data_start_time + '_' + str(i) + '.txt', 'w') as f:
            f.write("%s %s\n"%(start_gps, end_gps))
        with open('./DQ_unlens_precessing/V1_' + Data_start_time + '_' + str(i) + '.txt', 'w') as f:
            f.write("%s %s\n"%(start_gps, end_gps))
        
        excel_name_match_unlens.append(Data_start_time + '_' + str(i))
        
        Data_start_time_set_unlens.append(int(Data_start_time))

        merger_freq_set.append(merger_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]))
        ringdown_freq_set.append(ringdown_freq(SampleParam[1][i], SampleParam[2][i], SampleParam[0][i]))
        
        f0.write('../Paper4_CE_Modify/outdir_unlensed_precessing/outdir'+str(i) + '/H1_L1_V1_result.json'+'\n')


        
   
     
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/Data_start_time_set_unlens.csv', Data_start_time_set_unlens, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/SNR_selected_unlens.csv', SNR_selected_unlens, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/excel_name_match_unlens.csv', excel_name_match_unlens, delimiter=',', fmt='%s')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/chirp_mass_unlens.csv', Mz, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/mass_ratio_unlens.csv', mass_ratio, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/merger_freq_unlens.csv', merger_freq_set, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/ringdown_freq_unlens.csv', ringdown_freq_set, delimiter=',')
np.savetxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/reference_det.csv', reference_det, delimiter=',', fmt='%s')
  
    
    
    
