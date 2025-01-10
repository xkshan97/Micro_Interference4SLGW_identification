import numpy as np
import matplotlib.pyplot as plt
import os
import json

Mz_cut = 20

SNR_selected_unlens = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/match_filter_SNR_unlens.csv', delimiter=',')

SNR_selected_micro_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_precessing/match_filter_SNR.csv', delimiter=',')

Data_start_time_set_unlens = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/Data_start_time_set_unlens.csv', delimiter=',')
Data_start_time_set_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_precessing/inject_time_set.csv', delimiter=',')

Mz_unlens = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/chirp_mass_unlens.csv', delimiter=',')
Mz_micro_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_precessing/chirp_mass.csv', delimiter=',')
name_event = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/excel_name_match_unlens.csv', dtype=str)
name_event_micro = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_precessing/excel_name_match.csv', dtype=str)

#unlens
match_sum_unlens = []
match_mean_unlens = []
match_std_unlens = []
match_best_unlens = [] #不是一个好的统计量
SNR_tmp_unlens = [] #在测试的时候有的事例没有运行完，所以把运行完的事例的SNR单独存起来。
max_match_set = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/max_match_set_unlens.csv', delimiter=',')

#unlens total
with open('./Data_start_time_and_snr_and_match_and_outdir_name_precessing/unlens.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据
for i in range(len(Data_start_time_set_unlens)):

    try:
        "读取cWB重构结果"
        path_tmp = '../Paper5_PycWB/Unlens_precessing/' + name_event[i] +'/trigger'
        file_tmp = os.listdir(path_tmp)
        PycWB_res = np.loadtxt(path_tmp + '/' + file_tmp[0] + '/reconstructed_waveform_' + "H1" + '.txt')
        
    except Exception as e:
        print(path_tmp)
        continue
    
    
    match_sum_unlens.append(max_match_set[i])
    
    if Mz_unlens[i] > Mz_cut:
        SNR_tmp_unlens.append(SNR_selected_unlens[i])
        match_mean_unlens.append(np.mean(match_sum_unlens[-1]))
        match_std_unlens.append(np.std(match_sum_unlens[-1]))
        match_best_unlens.append(np.max(match_sum_unlens[-1]))
    



match_sum_lensed = []
match_mean_lensed = []
match_std_lensed = []
match_best_lensed = [] #不是一个好的统计量
SNR_tmp_lensed = [] #有的事例没有用cWB找到
theory_match = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name_precessing/match.csv', delimiter=',')
theory_match_tmp = []


with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name_precessing/with_micro.txt', 'r') as f:
    file_name = f.readlines() 
for i in range(len(Data_start_time_set_Total)):
    try:
        with open(file_name[i].split('\n')[0]) as f1:
            data_lensed = json.load(f1)
        
        parameters = list(data_lensed['posterior']['content'].keys())
    
        "读取cWB重构结果"
        path_tmp = '../Paper5_PycWB/Micro_precessing/' + name_event_micro[i] +'/trigger'
        file_tmp = os.listdir(path_tmp)
        PycWB_res = np.loadtxt(path_tmp + '/' + file_tmp[0] + '/reconstructed_waveform_' + 'H1' + '.txt')
        with open(path_tmp + '/' + file_tmp[0] + '/event.json') as event_josn:
            read_res = event_josn.readlines()
        SNR_tmp = np.sqrt(np.sum([float(snr_) for snr_ in read_res[0].split("snr")[1].split('[')[1].split(']')[0].split(',')]))
        
    except FileNotFoundError:
        continue
        
    try:
        match_sum_lensed.append(np.loadtxt('./match_result_Total_precessing/match_micro_'+file_name[i].split('/')[3][6::]+'.csv', delimiter=','))
    except FileNotFoundError:
        print('./match_result_Total_precessing/match_micro_'+file_name[i].split('/')[3][6::]+'.csv')
        continue
    
    if Mz_micro_Total[i] > Mz_cut:
        
        
        SNR_tmp_lensed.append(SNR_selected_micro_Total[i]) #注意这里用的是micro，因为我发现与macro差别还挺大的，因为有的微透镜会压低很多。
        theory_match_tmp.append(theory_match[i])
        match_mean_lensed.append(np.mean(match_sum_lensed[-1]))
        match_std_lensed.append(np.std(match_sum_lensed[-1]))
        match_best_lensed.append(np.max(match_sum_lensed[-1]))



#unlens envelop
up_sum = np.array(match_mean_unlens) + 1.5*np.array(match_std_unlens)
low_sum = np.array(match_mean_unlens) - 1.5*np.array(match_std_unlens)
Z1 = zip(SNR_tmp_unlens, up_sum)
Z1 = sorted(Z1)
Z2 = zip(SNR_tmp_unlens, low_sum)
Z2 = sorted(Z2)
Z3 = zip(SNR_tmp_unlens, match_mean_unlens)
Z3 = sorted(Z3)
Z4 = zip(SNR_tmp_unlens, match_std_unlens)
Z4 = sorted(Z4)
SNR_tmp_unlens, up_sum = zip(*Z1)
SNR_tmp_unlens, low_sum = zip(*Z2)
SNR_tmp_unlens, match_mean_unlens = zip(*Z3)
SNR_tmp_unlens, match_std_unlens = zip(*Z4)
up_SNR_envelope = [SNR_tmp_unlens[0]]
low_SNR_envelope = [SNR_tmp_unlens[-1]]
up_envelope = [up_sum[0]]
u_i = 0
low_envelope = [low_sum[-1]]
for i in range(1, len(SNR_tmp_unlens)):
    if up_sum[i]<up_envelope[u_i]:
        pass
    else:
        up_envelope.append(up_sum[i])
        up_SNR_envelope.append(SNR_tmp_unlens[i])
        u_i += 1
    if low_sum[len(SNR_tmp_unlens) - i - 1] > low_envelope[0]:
        pass
    else:
        low_envelope = [low_sum[len(SNR_tmp_unlens) - i - 1]] + low_envelope
        low_SNR_envelope = [SNR_tmp_unlens[len(SNR_tmp_unlens) - i - 1]] + low_SNR_envelope


#plot unlens total
# plt.scatter(SNR_tmp_unlens, match_mean_unlens, marker='s')
# plt.plot(up_SNR_envelope, up_envelope)
# plt.plot(low_SNR_envelope, low_envelope)
plt.fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='grey', label='chirp mass >' + str(int(Mz_cut)))
# plt.ylim(0.935,1)
plt.xlim(20, 700)
# plt.errorbar(SNR_tmp_unlens, match_mean_unlens, yerr = np.array(match_std_unlens), elinewidth=2, ecolor='red', capsize=5, capthick=1, linestyle='none')



plt.scatter(SNR_tmp_lensed, match_mean_lensed, c=theory_match, marker='s', cmap='viridis_r', vmin = 0.9, vmax = 1)
plt.errorbar(SNR_tmp_lensed, match_mean_lensed, yerr = np.array(match_std_lensed), elinewidth=2, ecolor='black', capsize=5, capthick=1, linestyle='none')

cb = plt.colorbar()
cb.set_label('Theory Match')
plt.xlabel('SNR')
plt.ylabel('match')
plt.grid()
# plt.ylim(0.95,1)
plt.xlim(30,600)
# plt.loglog()
plt.legend()
plt.savefig('./Plot_result_precessing/match' + str(int(Mz_cut)) + '.png', dpi=450)

np.savetxt('./Result_file4plot_precessing/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', low_SNR_envelope, delimiter=',')
np.savetxt('./Result_file4plot_precessing/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', up_SNR_envelope, delimiter=',')

np.savetxt('./Result_file4plot_precessing/low_envelope' + str(int(Mz_cut)) + '.csv', low_envelope, delimiter=',')
np.savetxt('./Result_file4plot_precessing/up_envelope' + str(int(Mz_cut)) + '.csv', up_envelope, delimiter=',')

np.savetxt('./Result_file4plot_precessing/SNR_lensed' + str(int(Mz_cut)) + '.csv', SNR_tmp_lensed, delimiter=',')
np.savetxt('./Result_file4plot_precessing/match_mean_lensed' + str(int(Mz_cut)) + '.csv', match_mean_lensed, delimiter=',')
np.savetxt('./Result_file4plot_precessing/theory_match' + str(int(Mz_cut)) + '.csv', theory_match, delimiter=',')
np.savetxt('./Result_file4plot_precessing/match_std_lensed' + str(int(Mz_cut)) + '.csv', match_std_lensed, delimiter=',')



# plt.fill_between(SNR_tmp_lensed003, np.array(match_mean_lensed003) + np.array(match_std_lensed003), np.array(match_mean_lensed003) - np.array(match_std_lensed003))

# plt.fill_between(SNR_tmp_lensed006, np.array(match_mean_lensed006) + np.array(match_std_lensed006), np.array(match_mean_lensed006) - np.array(match_std_lensed006))

# plt.fill_between(SNR_tmp_lensed009, np.array(match_mean_lensed009) + np.array(match_std_lensed009), np.array(match_mean_lensed009) - np.array(match_std_lensed009))


