import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

Mz_cut = 20
low_SNR_envelope = np.loadtxt('./Result_file4plot/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
up_SNR_envelope = np.loadtxt('./Result_file4plot/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')

low_envelope = np.loadtxt('./Result_file4plot/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
up_envelope = np.loadtxt('./Result_file4plot/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')

SNR_tmp_lensed = np.loadtxt('./Result_file4plot/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
match_mean_lensed = np.loadtxt('./Result_file4plot/match_mean_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
theory_match = np.loadtxt('./Result_file4plot/theory_match' + str(int(Mz_cut)) + '.csv', delimiter=',')
match_std_lensed = np.loadtxt('./Result_file4plot/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')

Data_start_time_set_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/inject_time_set.csv', delimiter=',')
SNR_selected_micro_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/match_filter_SNR.csv', delimiter=',')


with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/with_micro.txt', 'r') as f:
    file_name = f.readlines() 
    
lower_envelope_func = interp1d(low_SNR_envelope, low_envelope)
SNR_test = np.arange(low_SNR_envelope[0], low_SNR_envelope[-1], 0.1)

plt.plot(low_SNR_envelope, low_envelope)
plt.plot(SNR_test, lower_envelope_func(SNR_test), '--')

#找到在背景之下的事例
Selected_file_name = []
for i in range(len(file_name)):
    try:
        Tmp_low_envelope = lower_envelope_func(SNR_tmp_lensed[i])
        if Tmp_low_envelope > match_mean_lensed[i] + 1.5*match_std_lensed[i]:
            Selected_file_name.append(file_name[i].split('\n')[0])
        else:
            pass
    except ValueError:
        continue
        
#找到背景之下事例的index，去除多像造成的重复
Selected_index = []
Selected_file_name_not_duplication = []
for i in range(len(Selected_file_name)):
    tmp_index = Selected_file_name[i].split('/')[3][6::].split('_')[0]
    if tmp_index not in Selected_index:
        Selected_index.append(tmp_index)
        Selected_file_name_not_duplication.append(Selected_file_name[i])
        
np.savetxt('./PairResult/identified_index.csv', Selected_index, delimiter=',', fmt='%s')
np.savetxt('./PairResult/identified_file_name.csv', Selected_file_name, delimiter=',', fmt='%s')
np.savetxt('./PairResult/identified_file_name_not_duplication.csv', Selected_file_name_not_duplication, delimiter=',', fmt='%s')