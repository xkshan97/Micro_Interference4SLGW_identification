#GW随机抽样了100个位置
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import corner
from getdist import plots, MCSamples
import getdist
import matplotlib.lines as mlines
import IPython
from  multiprocessing import Process,Pool
from tqdm import tqdm
import gc

det = 'JWST'
#引力波的时间延迟
index_GW = 120
timedelay_GW_set = np.loadtxt('./Random_pin_down_GW/t_days.csv', delimiter=',')
S_host_quadruple = np.loadtxt('./SampleResult_GWGalaxy/S_host_quadruple.csv', delimiter=',')
S_unhost_quadruple = np.loadtxt('./SampleResult_Galaxy/S_unhost_quadruple.csv', delimiter=',')


def GW_TD_consistence_GG(index_GW_random):
    delta_t_GW_true = timedelay_GW_set[index_GW_random][1::] - timedelay_GW_set[index_GW_random][0]



    #unhost时间延迟
    Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_25p5_and_4_images = np.loadtxt('SampleResult_Galaxy/Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images.csv', delimiter=',')
    factor_unhost_10 = []
    factor_unhost_20 = []
    index_in_S_unhost_quadruple = 0
    for index in Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_25p5_and_4_images:
        index = int(index)
        print(index_in_S_unhost_quadruple)
        try:
            mcmc_new_list_unhost = np.loadtxt('./' + det + '_unhost_26_little_sigma/' + str(index) + 'corner_Res_data.csv', delimiter=',')
            # if len(mcmc_new_list_unhost) < 5000:
            #     pass
            # else:
            #     mcmc_new_list_unhost = mcmc_new_list_unhost[0:5000]
            #labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
            t0_unhost = mcmc_new_list_unhost[:,5]
            t1_unhost = mcmc_new_list_unhost[:,6]
            t2_unhost = mcmc_new_list_unhost[:,7]
            t3_unhost = mcmc_new_list_unhost[:,8]
            SFR_unhost = mcmc_new_list_unhost[:,9]
            delta_t_10_unhost = t1_unhost - t0_unhost
            delta_t_20_unhost = t2_unhost - t0_unhost
            delta_t_30_unhost = t3_unhost - t0_unhost
            del mcmc_new_list_unhost, t0_unhost, t1_unhost, t2_unhost, t3_unhost
            gc.collect()

            
            #未加权的
            delta_t_ratio_10_20_unhost = delta_t_10_unhost / delta_t_20_unhost
            delta_t_ratio_10_30_unhost = delta_t_10_unhost / delta_t_30_unhost
            
            factor_tmp_1_index_unhost = np.where((delta_t_ratio_10_20_unhost <= delta_t_GW_true[0]/delta_t_GW_true[1] + 0.001)&((delta_t_ratio_10_20_unhost >= delta_t_GW_true[0]/delta_t_GW_true[1] - 0.001)))[0]
            factor_tmp_2_index_unhost = np.where((delta_t_ratio_10_30_unhost <= delta_t_GW_true[0]/delta_t_GW_true[2] + 0.001)&((delta_t_ratio_10_30_unhost >= delta_t_GW_true[0]/delta_t_GW_true[2] - 0.001)))[0]
            
             

            factor_tmp_1_unhost = len(factor_tmp_1_index_unhost) / len(delta_t_ratio_10_20_unhost) * S_unhost_quadruple[index_in_S_unhost_quadruple] * np.sum(SFR_unhost[factor_tmp_1_index_unhost]) / np.sum(SFR_unhost)
            factor_tmp_2_unhost = len(factor_tmp_2_index_unhost) / len(delta_t_ratio_10_30_unhost) * S_unhost_quadruple[index_in_S_unhost_quadruple] * np.sum(SFR_unhost[factor_tmp_2_index_unhost]) / np.sum(SFR_unhost)
            factor_unhost_10.append(factor_tmp_1_unhost)
            factor_unhost_20.append(factor_tmp_2_unhost)
            index_in_S_unhost_quadruple += 1
            del delta_t_10_unhost, delta_t_20_unhost, delta_t_30_unhost, delta_t_ratio_10_20_unhost, delta_t_ratio_10_30_unhost
            gc.collect()
            
        except FileNotFoundError:
            factor_unhost_10.append(0)
            factor_unhost_20.append(0)
            index_in_S_unhost_quadruple += 1
        

    factor_unhost_10 = np.array(factor_unhost_10) 
    factor_unhost_20 = np.array(factor_unhost_20)
    factor_unhost_sum = factor_unhost_10 * factor_unhost_20

    #host时间延迟
    factor_host_10 = []
    factor_host_20 = []
    mag_Res = np.loadtxt('./SampleResult_GWGalaxy/mag_Res_' + str(index_GW) + '.csv', delimiter=',')
    for index_tmp in tqdm(range(len(mag_Res))): #分别估计宿主星系星等最小和最大的情况
        try:
            mcmc_new_list_host = np.loadtxt('./' + det + '_host_pop_little_sigma/' + str(index_GW) + '_' + str(index_tmp) + 'corner_Res_data.csv', delimiter=',')
            # if len(mcmc_new_list_host) < 5000:
            #     pass
            # else:
            #     mcmc_new_list_host = mcmc_new_list_host[0:5000]
            #labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
            t0_host = mcmc_new_list_host[:,5]
            t1_host = mcmc_new_list_host[:,6]
            t2_host = mcmc_new_list_host[:,7]
            t3_host = mcmc_new_list_host[:,8]
            SFR_host = mcmc_new_list_host[:,9]
            delta_t_10_host = t1_host - t0_host
            delta_t_20_host = t2_host - t0_host
            delta_t_30_host = t3_host - t0_host
            del mcmc_new_list_host, t0_host, t1_host, t2_host, t3_host
            gc.collect()
            #未加权的
            delta_t_ratio_10_20_host = delta_t_10_host / delta_t_20_host
            delta_t_ratio_10_30_host = delta_t_10_host / delta_t_30_host
            
            factor_tmp_1_index_host = np.where((delta_t_ratio_10_20_host <= delta_t_GW_true[0]/delta_t_GW_true[1] + 0.001)&((delta_t_ratio_10_20_host >= delta_t_GW_true[0]/delta_t_GW_true[1] - 0.001)))[0]
            factor_tmp_2_index_host = np.where((delta_t_ratio_10_30_host <= delta_t_GW_true[0]/delta_t_GW_true[2] + 0.001)&((delta_t_ratio_10_30_host >= delta_t_GW_true[0]/delta_t_GW_true[2] - 0.001)))[0]
            
            
            factor_tmp_1_host = len(factor_tmp_1_index_host) / len(delta_t_ratio_10_20_host) * S_host_quadruple[index_tmp] * np.sum(SFR_host[factor_tmp_1_index_host]) / np.sum(SFR_host)
            factor_tmp_2_host = len(factor_tmp_2_index_host) / len(delta_t_ratio_10_30_host) * S_host_quadruple[index_tmp] * np.sum(SFR_host[factor_tmp_2_index_host]) / np.sum(SFR_host)
            factor_host_10.append(factor_tmp_1_host)
            factor_host_20.append(factor_tmp_2_host)
            del delta_t_10_host, delta_t_20_host, delta_t_30_host, delta_t_ratio_10_20_host, delta_t_ratio_10_30_host
            gc.collect()
        except FileNotFoundError:
            factor_host_10.append(0)
            factor_host_20.append(0)
        

    factor_host_10 = np.array(factor_host_10) 
    factor_host_20 = np.array(factor_host_20)
    factor_host_sum = factor_host_10 * factor_host_20
    JWST_factor_host_sum_weighted = []
    for tmp_i in range(len(mag_Res)):
        if mag_Res[tmp_i][0] < 26:
            weighted_factor = int(mag_Res[tmp_i][5] / np.min(mag_Res[:,5]))
            for tmp_j in range(weighted_factor):
                JWST_factor_host_sum_weighted.append(factor_host_sum[tmp_i])
                






    #贝叶斯因子画图
    plt.hist(factor_unhost_sum, bins=1000, histtype="step", cumulative=-1, label='unhost', density='True')
    plt.hist(factor_host_sum, bins=1000, histtype="step", cumulative=-1, label='host', density='True')
    plt.xlabel('Bayes factor')
    plt.ylabel('PDF')
    plt.legend()
    plt.grid()
    # plt.ylim(0, 100)
    plt.savefig('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_statistc_long_time_S_area.png', dpi=450)
    plt.close()
    plt.hist(factor_unhost_sum, bins=1000, histtype="step", cumulative=-1, label='unhost', density='True')
    plt.hist(JWST_factor_host_sum_weighted, bins=1000, histtype="step", cumulative=-1, label='host', density='True')
    plt.xlabel('Bayes factor')
    plt.ylabel('PDF')
    plt.legend()
    plt.grid()
    # plt.ylim(0, 100)
    plt.savefig('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_statistc_long_time_SFR_weighted.png', dpi=450)
    plt.close()
    np.savetxt('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_factor_unhost_long_time.csv', factor_unhost_sum, delimiter=',')
    np.savetxt('./IntermediaPlot/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time.csv', factor_host_sum, delimiter=',')
    np.savetxt('./IntermediaPlot/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time_SFR_weighted.csv', JWST_factor_host_sum_weighted, delimiter=',')


if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = len(timedelay_GW_set)
    length = threadcount#int(threadcount * 150) #
    # threadcount = 1
    
    
    percount = length//threadcount

    print("percount = ", percount)

    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
    
        res = pool.apply_async(func=GW_TD_consistence_GG, args=(i,))
        
    pool.close()
    pool.join()

color_set = ['indianred','royalblue','goldenrod', 'black', 'seagreen']
#所有100个采样点合起来的图
factor_unhost_sum = []
JWST_factor_host_sum_weighted = []
for index_GW_random in range(len(timedelay_GW_set)):
    factor_unhost_sum.extend(np.loadtxt('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_factor_unhost_long_time.csv', delimiter=','))
    JWST_factor_host_sum_weighted.extend(np.loadtxt('./IntermediaPlot/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time_SFR_weighted.csv', delimiter=','))
mean_unhost = np.percentile(factor_unhost_sum, 50)
percentile_unhost = np.percentile(factor_unhost_sum, [2.2, 97.8])
error_range_unhost = [[mean_unhost - percentile_unhost[0]], [percentile_unhost[1] - mean_unhost]]
mean_host = np.percentile(JWST_factor_host_sum_weighted, 50)
percentile_host = np.percentile(JWST_factor_host_sum_weighted, [2.2, 97.8])
error_range_host = [[mean_host - percentile_host[0]], [percentile_host[1] - mean_host]]
plt.errorbar(0, mean_host, yerr=error_range_host, fmt='o:', capsize=3, color=color_set[0], label='host')
plt.errorbar(1, mean_unhost, yerr=error_range_unhost, fmt='o:', capsize=3, color=color_set[1], label='unhost')
# plt.semilogy()
plt.xlim(-5, 5)
plt.xticks([])
plt.ylabel('Area[kpc$^2$]')
plt.grid()
plt.savefig('./IntermediaPlot/Area_sum_100_random_point.png', dpi=450)
plt.close()

#100个采样点分开画
plt.figure(figsize=(20,4))
unhost_lower = []
unhost_upper = []
host_lower = []
host_upper = []
for index_GW_random in range(len(timedelay_GW_set)):
    factor_unhost_sum = np.loadtxt('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_factor_unhost_long_time.csv', delimiter=',')[10:15]
    JWST_factor_host_sum_weighted = np.loadtxt('./IntermediaPlot/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time_SFR_weighted.csv', delimiter=',')
    mean_unhost = np.percentile(factor_unhost_sum, 50)
    percentile_unhost = np.percentile(factor_unhost_sum, [15.87, 84.13])
    error_range_unhost = [[mean_unhost - percentile_unhost[0]], [percentile_unhost[1] - mean_unhost]]
    mean_host = np.percentile(JWST_factor_host_sum_weighted, 50)
    percentile_host = np.percentile(JWST_factor_host_sum_weighted, [15.87, 84.13])
    error_range_host = [[mean_host - percentile_host[0]], [percentile_host[1] - mean_host]]
    unhost_lower.append(percentile_unhost[0])
    unhost_upper.append(percentile_unhost[1])
    host_lower.append(percentile_host[0])
    host_upper.append(percentile_host[1])

    plt.errorbar(index_GW_random, mean_host, yerr=error_range_host, fmt='o:', capsize=3, color=color_set[0])
    # plt.errorbar(index_GW_random, mean_unhost, yerr=error_range_unhost, fmt='o:', capsize=3, color=color_set[1], label='unhost')

plt.errorbar(index_GW_random, mean_host, yerr=error_range_host, fmt='o:', capsize=3, color=color_set[0], label='host')
plt.semilogy()
# plt.errorbar(index_GW_random, mean_unhost, yerr=error_range_unhost, fmt='o:', capsize=3, color=color_set[1], label='unhost')
plt.fill_between(range(len(timedelay_GW_set)), unhost_lower, unhost_upper, color = color_set[1], label='unhost')
# plt.fill_between(range(len(timedelay_GW_set)), unhost_lower, unhost_upper)
plt.xlim(-1, 101)
# plt.xticks([])
plt.ylabel('Area[kpc$^2$]')
plt.grid()
plt.legend()

plt.savefig('./IntermediaPlot/Area_100_random_point.png', dpi=450)
plt.close()


#计算误报率
FAP_peryear1sigma = []
FAP_peryear2sigma = []
FAP_peryear3sigma = []
year_obs = 11
for index_GW_random in range(len(timedelay_GW_set)):
    factor_unhost_sum = np.loadtxt('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_factor_unhost_long_time.csv', delimiter=',')
    JWST_factor_host_sum_weighted = np.loadtxt('./IntermediaPlot/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time_SFR_weighted.csv', delimiter=',')
    
    mean_host = np.percentile(JWST_factor_host_sum_weighted, 50)
    percentile_host = np.percentile(JWST_factor_host_sum_weighted, [15.87, 2.2, 0.13])
    lower_limit1 = percentile_host[0]
    lower_limit2 = percentile_host[1]
    lower_limit3 = percentile_host[2]
    FAP_peryear1sigma.append(np.sum(factor_unhost_sum > lower_limit1) / year_obs)
    FAP_peryear2sigma.append(np.sum(factor_unhost_sum > lower_limit2) / year_obs)
    FAP_peryear3sigma.append(np.sum(factor_unhost_sum > lower_limit3) / year_obs)

plt.plot(FAP_peryear1sigma, '--', label='1$\sigma$') 
plt.plot(FAP_peryear2sigma, '--', label='2$\sigma$') 
plt.plot(FAP_peryear3sigma, '--', label='3$\sigma$') 
plt.legend()
plt.xlim(-1, 101)
# plt.xticks([])
plt.ylabel('FAR per year')
plt.grid()

plt.savefig('./IntermediaPlot/FAR_per_year.png', dpi=450)
plt.close()

    
#找出来能mimic最小星等的事例
# index_mimic_min_range = np.where(factor_unhost_sum>factor_host_sum[0])[0]
# index_mimic_max_range = np.where(factor_unhost_sum>factor_host_sum[1])[0]

# mimic_event_index = 0
# corner_max_mimic = [[], [], [], [], []]
# for index_tmp in index_mimic_max_range:
#     index = Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_25p5_and_4_images[index_tmp]
#     mcmc_new_list_unhost = np.loadtxt('./Euclid_unhost_25p5/' + str(int(index)) + 'corner_Res_data.csv', delimiter=',')
#     #labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
#     t0_unhost = mcmc_new_list_unhost[:,5]
#     t1_unhost = mcmc_new_list_unhost[:,6]
#     t2_unhost = mcmc_new_list_unhost[:,7]
#     t3_unhost = mcmc_new_list_unhost[:,8]
#     delta_t_10_unhost = t1_unhost - t0_unhost
#     delta_t_20_unhost = t2_unhost - t0_unhost
#     delta_t_30_unhost = t3_unhost - t0_unhost
#     delta_t_ratio_10_20_unhost = delta_t_10_unhost / delta_t_20_unhost
#     delta_t_ratio_10_30_unhost = delta_t_10_unhost / delta_t_30_unhost
    
#     for j in range(len(delta_t_ratio_10_20_unhost)):
#         corner_max_mimic[mimic_event_index].append([delta_t_ratio_10_20_unhost[j], delta_t_ratio_10_30_unhost[j]])
    
#     mimic_event_index += 1



# mcmc_new_list_host = np.loadtxt('./Euclid_host_25p5/' + str(index_GW) + '_mag_mincorner_Res_data.csv', delimiter=',')
# #labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
# t0_host = mcmc_new_list_host[:,5]
# t1_host = mcmc_new_list_host[:,6]
# t2_host = mcmc_new_list_host[:,7]
# t3_host = mcmc_new_list_host[:,8]
# delta_t_10_host = t1_host - t0_host
# delta_t_20_host = t2_host - t0_host
# delta_t_30_host = t3_host - t0_host
# delta_t_ratio_10_20_host = delta_t_10_host / delta_t_20_host
# delta_t_ratio_10_30_host = delta_t_10_host / delta_t_30_host
# corner_min = []
# for j in range(len(delta_t_ratio_10_20_host)):
#     corner_min.append([delta_t_ratio_10_20_host[j], delta_t_ratio_10_30_host[j]])


# mcmc_new_list_host = np.loadtxt('./Euclid_host_25p5/' + str(index_GW) + '_mag_maxcorner_Res_data.csv', delimiter=',')
# #labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
# t0_host = mcmc_new_list_host[:,5]
# t1_host = mcmc_new_list_host[:,6]
# t2_host = mcmc_new_list_host[:,7]
# t3_host = mcmc_new_list_host[:,8]
# delta_t_10_host = t1_host - t0_host
# delta_t_20_host = t2_host - t0_host
# delta_t_30_host = t3_host - t0_host
# delta_t_ratio_10_20_host = delta_t_10_host / delta_t_20_host
# delta_t_ratio_10_30_host = delta_t_10_host / delta_t_30_host
# corner_max = []
# for j in range(len(delta_t_ratio_10_20_host)):
#     corner_max.append([delta_t_ratio_10_20_host[j], delta_t_ratio_10_30_host[j]])
    



# labels = ['\\frac{\Delta\;t_{10}}{\Delta\;t_{20}}', '\\frac{\Delta\;t_{10}}{\Delta\;t_{30}}']
# names = ['\\frac{\Delta\;t_{10}}{\Delta\;t_{20}}', '\\frac{\Delta\;t_{10}}{\Delta\;t_{30}}']

# samples1 = MCSamples(samples=np.array(corner_min), names = names, labels = labels, label='mag min')
# samples2= MCSamples(samples=np.array(corner_max),names = names, labels = labels,label='mag min')

# samples1_mimic = MCSamples(samples=np.array(corner_max_mimic[0]), names = names, labels = labels, label='mimic $1$')
# samples2_mimic= MCSamples(samples=np.array(corner_max_mimic[1]),names = names, labels = labels,label='mimic $2$')
# samples3_mimic = MCSamples(samples=np.array(corner_max_mimic[2]), names = names, labels = labels, label='mimic $3$')
# samples4_mimic= MCSamples(samples=np.array(corner_max_mimic[3]),names = names, labels = labels,label='mimic $4$')
# samples5_mimic = MCSamples(samples=np.array(corner_max_mimic[4]), names = names, labels = labels, label='mimic $5$')
# # samples31= MCSamples(samples=samples2,names = names1, labels = labels1,label='QSO[XUV]+QSO[AS]+BAO',settings={'smooth_scale_2D':0.7})

# g = plots.get_subplot_plotter()
# samples1.updateSettings({'contours': [0.6826, 0.9544]})
# samples2.updateSettings({'contours': [0.6826, 0.9544]})
# samples1_mimic.updateSettings({'contours': [0.6826, 0.9544]})
# samples2_mimic.updateSettings({'contours': [0.6826, 0.9544]})
# samples3_mimic.updateSettings({'contours': [0.6826, 0.9544]})
# samples4_mimic.updateSettings({'contours': [0.6826, 0.9544]})
# samples5_mimic.updateSettings({'contours': [0.6826, 0.9544]})
# # samples31.updateSettings({'contours': [0.6826, 0.9544]})

# g.settings.legend_fontsize = 10
# g.settings.axes_fontsize = 15
# g.settings.axes_labelsize = 15
# g.settings.axis_tick_x_rotation = 45
# g.triangle_plot([samples1_mimic, samples2_mimic, samples3_mimic, samples4_mimic, samples5_mimic, samples1, samples2], colors=['grey', 'grey', 'grey', 'grey', 'grey', color_set[0], color_set[1]],filled=True,legend_loc='upper right',
#                 line_args=[{'ls':'--','color':'grey'},{'ls':'--','color':'grey'},{'ls':'--','color':'grey'},{'ls':'--','color':'grey'},{'ls':'--','color':'grey'},{'color':color_set[0]},{'color':color_set[1]}], markers={'\\frac{\Delta\;t_{10}}{\Delta\;t_{20}}':delta_t_GW_true[0]/delta_t_GW_true[1], '\\frac{\Delta\;t_{10}}{\Delta\;t_{30}}':delta_t_GW_true[0]/delta_t_GW_true[2]}
#             ,marker_args={'color':'black','lw':20,'lw':1} ,contour_lws=[1,1,1,1,1,1.5,1.5], shaded=True)

# plt.savefig('./IntermediaPlot/Bayes_corner.png', dpi=450)