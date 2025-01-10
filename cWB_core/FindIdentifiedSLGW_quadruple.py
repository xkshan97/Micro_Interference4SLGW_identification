import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import json
import os

Mz_cut = 20
low_SNR_envelope = np.loadtxt('./Result_file4plot_quadruple/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
up_SNR_envelope = np.loadtxt('./Result_file4plot_quadruple/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')

low_envelope = np.loadtxt('./Result_file4plot_quadruple/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
up_envelope = np.loadtxt('./Result_file4plot_quadruple/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')

SNR_tmp_lensed = np.loadtxt('./Result_file4plot_quadruple/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
match_mean_lensed = np.loadtxt('./Result_file4plot_quadruple/match_mean_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
theory_match = np.loadtxt('./Result_file4plot_quadruple/theory_match' + str(int(Mz_cut)) + '.csv', delimiter=',')
match_std_lensed = np.loadtxt('./Result_file4plot_quadruple/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')

Data_start_time_set_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_quadruple/inject_time_set.csv', delimiter=',')
SNR_selected_micro_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_quadruple/match_filter_SNR.csv', delimiter=',')


with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name_quadruple/with_micro.txt', 'r') as f:
    file_name = f.readlines() 
    
lower_envelope_func = interp1d(low_SNR_envelope, low_envelope)
SNR_test = np.arange(low_SNR_envelope[0], low_SNR_envelope[-1], 0.1)

plt.plot(low_SNR_envelope, low_envelope)
plt.plot(SNR_test, lower_envelope_func(SNR_test), '--')

#找到在背景之下的事例
Selected_file_name = []
read_index = 0
for i in range(len(file_name)):
    try:
        with open(file_name[i].split('\n')[0]) as f1:
            data_lens = json.load(f1)
        
    except FileNotFoundError:
        print(file_name[i].split('\n')[0])
        continue
    except KeyError:
        print(file_name[i].split('\n')[0])
        continue
        
    try:
        "读取cWB重构结果"
        cWB_path = './cWB_micro_quadruple/data/' + str(int(Data_start_time_set_Total[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/H1_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
    except FileNotFoundError:
        print(cWB_path)
        continue
    try:
        tmp = np.loadtxt('./match_result_Total_quadruple/match_micro_'+file_name[i].split('/')[3][6::]+'.csv', delimiter=',')
    except FileNotFoundError:
        print('./match_result_Total_quadruple/match_micro_'+file_name[i].split('/')[3][6::]+'.csv')
        continue
    
    try:
        Tmp_low_envelope = lower_envelope_func(SNR_tmp_lensed[read_index])
        if Tmp_low_envelope > match_mean_lensed[read_index] + 1.5*match_std_lensed[read_index]:
            Selected_file_name.append(file_name[i].split('\n')[0])
        else:
            pass
        read_index += 1
    except ValueError:
        read_index += 1
        continue
        
#找到背景之下事例的index，去除多像造成的重复
Selected_index = []
Selected_file_name_not_duplication = []
for i in range(len(Selected_file_name)):
    tmp_index = Selected_file_name[i].split('/')[3][6::].split('_')[0]
    if tmp_index not in Selected_index:
        Selected_index.append(tmp_index)
        Selected_file_name_not_duplication.append(Selected_file_name[i])
        
np.savetxt('./PairResult_quadruple/identified_index.csv', Selected_index, delimiter=',', fmt='%s')
np.savetxt('./PairResult_quadruple/identified_file_name.csv', Selected_file_name, delimiter=',', fmt='%s')
np.savetxt('./PairResult_quadruple/identified_file_name_not_duplication.csv', Selected_file_name_not_duplication, delimiter=',', fmt='%s')