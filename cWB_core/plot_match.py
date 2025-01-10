import numpy as np
import matplotlib.pyplot as plt
import os
import json

Mz_cut = 20

SNR_selected_unlens = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name/match_filter_SNR_unlens.csv', delimiter=',')

SNR_selected_micro_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/match_filter_SNR.csv', delimiter=',')

Data_start_time_set_unlens = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name/Data_start_time_set_unlens.csv', delimiter=',')
Data_start_time_set_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/inject_time_set.csv', delimiter=',')

Mz_unlens = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name/chirp_mass_unlens.csv', delimiter=',')
Mz_micro_Total = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/chirp_mass.csv', delimiter=',')


#unlens
match_sum_unlens = []
match_mean_unlens = []
match_std_unlens = []
match_best_unlens = [] #不是一个好的统计量
SNR_tmp_unlens = [] #在测试的时候有的事例没有运行完，所以把运行完的事例的SNR单独存起来。

#unlens total
with open('./Data_start_time_and_snr_and_match_and_outdir_name/unlens.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据
for i in range(len(Data_start_time_set_unlens)):
    try:
        with open(file_name[i].split('\n')[0]) as f1:
            data_unlens = json.load(f1)
        parameters = list(data_unlens['posterior']['content'].keys())
    except FileNotFoundError or KeyError:
        print(file_name[i].split('\n')[0])
        continue

    try:
        "读取cWB重构结果"
        cWB_path = './cWB_unlens/data/' + str(int(Data_start_time_set_unlens[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/H1_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
    except FileNotFoundError:
        print(cWB_path)
        continue
    
    try:
        match_sum_unlens.append(np.loadtxt('./match_result/match_unlens_'+file_name[i].split('/')[3][6::]+'.csv', delimiter=','))
    except FileNotFoundError:
        print('./match_result/match_unlens_'+file_name[i].split('/')[3][6::]+'.csv')
        continue
    if Mz_unlens[i] > Mz_cut:
        SNR_tmp_unlens.append(SNR_selected_unlens[i])
        match_mean_unlens.append(np.mean(match_sum_unlens[-1]))
        match_std_unlens.append(np.std(match_sum_unlens[-1]))
        match_best_unlens.append(np.max(match_sum_unlens[-1]))
    


#unlens noise (no cwb)
match_sum_unlens_noise = []
match_mean_unlens_noise = []
match_std_unlens_noise = []
match_best_unlens_noise = [] #不是一个好的统计量
SNR_tmp_unlens_noise = [] #在测试的时候有的事例没有运行完，所以把运行完的事例的SNR单独存起来。
match_noise_unlens = np.loadtxt('./match_noisy_inject_VS_Bilby/unlens_match_noise.csv', delimiter=',')

#unlens total
with open('./Data_start_time_and_snr_and_match_and_outdir_name/unlens.txt', 'r') as f:
    file_name = f.readlines()
#读取out_dir中的数据

index_in_match = 0
for i in range(len(Data_start_time_set_unlens)):
    
    if os.path.exists(file_name[i].split('\n')[0]):
        
        match_sum_unlens_noise.append(match_noise_unlens[index_in_match])
        
        if Mz_unlens[i] > Mz_cut:
            SNR_tmp_unlens_noise.append(SNR_selected_unlens[i])
            match_mean_unlens_noise.append(np.mean(match_sum_unlens_noise[-1]))
            match_std_unlens_noise.append(np.std(match_sum_unlens_noise[-1]))
            match_best_unlens_noise.append(np.max(match_sum_unlens_noise[-1]))
        
        index_in_match += 1
        
up_sum_noise = np.array(match_mean_unlens_noise) + 1.5*np.array(match_std_unlens_noise)
low_sum_noise = np.array(match_mean_unlens_noise) - 1.5*np.array(match_std_unlens_noise)

Z1 = zip(SNR_tmp_unlens_noise, up_sum_noise)
Z1 = sorted(Z1)
Z2 = zip(SNR_tmp_unlens_noise, low_sum_noise)
Z2 = sorted(Z2)
Z3 = zip(SNR_tmp_unlens_noise, match_mean_unlens_noise)
Z3 = sorted(Z3)
Z4 = zip(SNR_tmp_unlens_noise, match_std_unlens_noise)
Z4 = sorted(Z4)
SNR_tmp_unlens_noise, up_sum_noise = zip(*Z1)
SNR_tmp_unlens_noise, low_sum_noise = zip(*Z2)
SNR_tmp_unlens_noise, match_mean_unlens_noise = zip(*Z3)
SNR_tmp_unlens_noise, match_std_unlens_noise = zip(*Z4)
up_SNR_envelope_noise = [SNR_tmp_unlens_noise[0]]
low_SNR_envelope_noise = [SNR_tmp_unlens_noise[-1]]
up_envelope_noise = [up_sum_noise[0]]
u_i = 0
low_envelope_noise = [low_sum_noise[-1]]
for i in range(1, len(SNR_tmp_unlens_noise)):
    if up_sum_noise[i]<up_envelope_noise[u_i]:
        pass
    else:
        up_envelope_noise.append(up_sum_noise[i])
        up_SNR_envelope_noise.append(SNR_tmp_unlens_noise[i])
        u_i += 1
    if low_sum_noise[len(SNR_tmp_unlens_noise) - i - 1] > low_envelope_noise[0]:
        pass
    else:
        low_envelope_noise = [low_sum_noise[len(SNR_tmp_unlens_noise) - i - 1]] + low_envelope_noise
        low_SNR_envelope_noise = [SNR_tmp_unlens_noise[len(SNR_tmp_unlens_noise) - i - 1]] + low_SNR_envelope_noise
plt.fill_between(low_SNR_envelope_noise, low_envelope_noise, 1, interpolate=True, color='grey', label='chirp mass >' + str(int(Mz_cut)))
plt.scatter(SNR_tmp_unlens_noise, match_mean_unlens_noise)
plt.ylim(0.935,1)
plt.xlim(20, 700)
np.savetxt('./match_noisy_inject_VS_Bilby/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', low_SNR_envelope_noise, delimiter=',')
np.savetxt('./match_noisy_inject_VS_Bilby/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', up_SNR_envelope_noise, delimiter=',')

np.savetxt('./match_noisy_inject_VS_Bilby/low_envelope' + str(int(Mz_cut)) + '.csv', low_envelope_noise, delimiter=',')
np.savetxt('./match_noisy_inject_VS_Bilby/up_envelope' + str(int(Mz_cut)) + '.csv', up_envelope_noise, delimiter=',')





#lens
match_sum_lensed = []
match_mean_lensed = []
match_std_lensed = []
match_best_lensed = [] #不是一个好的统计量
SNR_tmp_lensed = [] #有的事例没有用cWB找到
theory_match = np.loadtxt('./Total_Data_start_time_and_snr_and_match_and_outdir_name/match.csv', delimiter=',')
theory_match_tmp = []


with open('./Total_Data_start_time_and_snr_and_match_and_outdir_name/with_micro.txt', 'r') as f:
    file_name = f.readlines() 
for i in range(len(Data_start_time_set_Total)):
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
        cWB_path = './cWB_micro/data/' + str(int(Data_start_time_set_Total[i])) + '_' + file_name[i].split('/')[3][6::]
        cWB_path += '/' + os.listdir(cWB_path)[-1] + '/H1_wf_strain.dat'
        with open(cWB_path) as f2:
            data = f2.readlines()
    except FileNotFoundError:
        print(cWB_path)
        continue
    try:
        match_sum_lensed.append(np.loadtxt('./match_result_Total/match_micro_'+file_name[i].split('/')[3][6::]+'.csv', delimiter=','))
    except FileNotFoundError:
        print('./match_result_Total/match_micro_'+file_name[i].split('/')[3][6::]+'.csv')
        continue
    
    if Mz_micro_Total[i] > Mz_cut:
        
        
        SNR_tmp_lensed.append(SNR_selected_micro_Total[i]) #注意这里用的是micro，因为我发现与macro差别还挺大的，因为有的微透镜会压低很多。
        theory_match_tmp.append(theory_match[i])
        match_mean_lensed.append(np.mean(match_sum_lensed[-1]))
        match_std_lensed.append(np.std(match_sum_lensed[-1]))
        match_best_lensed.append(np.max(match_sum_lensed[-1]))


#big sample
match_filter_SNR_1000 = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/match_filter_SNR_unlens.csv', delimiter=',')
max_match_set_1000 = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/max_match_set_unlens.csv', delimiter=',')
chirp_mass_1000 = np.loadtxt('./Data_start_time_and_snr_and_match_and_outdir_name_big_sample/chirp_mass_unlens.csv')
index_mass_gtr_20 = np.where((chirp_mass_1000 > 20))[0]
index_mass_gtr_20_final = []
for tmp_i in index_mass_gtr_20:
    if match_filter_SNR_1000[tmp_i] > 300:
        pass
    elif match_filter_SNR_1000[tmp_i] > 200 and max_match_set_1000[tmp_i] < 0.995:
        pass
    elif match_filter_SNR_1000[tmp_i] < 200 and match_filter_SNR_1000[tmp_i] > 100 and max_match_set_1000[tmp_i] < 0.985:
        pass
    elif match_filter_SNR_1000[tmp_i] > 50 and match_filter_SNR_1000[tmp_i] < 100 and max_match_set_1000[tmp_i] < 0.975:
        pass
    else:
        index_mass_gtr_20_final.append(tmp_i)
index_mass_gtr_20_final = np.array(index_mass_gtr_20_final)
match_filter_SNR_1000 = match_filter_SNR_1000[index_mass_gtr_20_final]
max_match_set_1000 = max_match_set_1000[index_mass_gtr_20_final]

#unlens envelop
match_mean_unlens = np.append(match_mean_unlens, max_match_set_1000)
match_std_unlens = np.append(match_std_unlens, np.zeros(len(max_match_set_1000)))
SNR_tmp_unlens = np.append(SNR_tmp_unlens, match_filter_SNR_1000)

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
plt.ylim(0.935,1)
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
plt.savefig('./Plot_result/match' + str(int(Mz_cut)) + '.png', dpi=450)

np.savetxt('./Result_file4plot/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', low_SNR_envelope, delimiter=',')
np.savetxt('./Result_file4plot/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', up_SNR_envelope, delimiter=',')

np.savetxt('./Result_file4plot/low_envelope' + str(int(Mz_cut)) + '.csv', low_envelope, delimiter=',')
np.savetxt('./Result_file4plot/up_envelope' + str(int(Mz_cut)) + '.csv', up_envelope, delimiter=',')

np.savetxt('./Result_file4plot/SNR_lensed' + str(int(Mz_cut)) + '.csv', SNR_tmp_lensed, delimiter=',')
np.savetxt('./Result_file4plot/match_mean_lensed' + str(int(Mz_cut)) + '.csv', match_mean_lensed, delimiter=',')
np.savetxt('./Result_file4plot/theory_match' + str(int(Mz_cut)) + '.csv', theory_match, delimiter=',')
np.savetxt('./Result_file4plot/match_std_lensed' + str(int(Mz_cut)) + '.csv', match_std_lensed, delimiter=',')



# plt.fill_between(SNR_tmp_lensed003, np.array(match_mean_lensed003) + np.array(match_std_lensed003), np.array(match_mean_lensed003) - np.array(match_std_lensed003))

# plt.fill_between(SNR_tmp_lensed006, np.array(match_mean_lensed006) + np.array(match_std_lensed006), np.array(match_mean_lensed006) - np.array(match_std_lensed006))

# plt.fill_between(SNR_tmp_lensed009, np.array(match_mean_lensed009) + np.array(match_std_lensed009), np.array(match_mean_lensed009) - np.array(match_std_lensed009))


