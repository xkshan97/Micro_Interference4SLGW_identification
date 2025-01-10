#GW在宿主星系中心
# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import corner
from getdist import plots, MCSamples
import getdist
import matplotlib.lines as mlines
import IPython
import gc
from tqdm import tqdm

SFR_weight_in_Galaxy = 'No' #如果是No的话，就把后面的SFR加权去掉，如果是''，则加上SFR加权
det = 'JWST'
#引力波的时间延迟
index_GW = 120
timedelay_GW_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/timedelayDays.csv', delimiter=',')
imagenum_GW_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
GW_index_start_in_td = int(np.sum(imagenum_GW_set[0:120]))
GW_index_end_in_td = GW_index_start_in_td + 4
td_GW_true = timedelay_GW_set[GW_index_start_in_td: GW_index_end_in_td]
delta_t_GW_true = td_GW_true[1:] - td_GW_true[0]
S_host_quadruple = np.loadtxt('./SampleResult_GWGalaxy/S_host_quadruple.csv', delimiter=',')
S_unhost_quadruple = np.loadtxt('./SampleResult_Galaxy/S_unhost_quadruple.csv', delimiter=',')




color_set = ['indianred','royalblue','goldenrod', 'black', 'seagreen']

#unhost时间延迟
Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_25p5_and_4_images = np.loadtxt('SampleResult_Galaxy/Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images.csv', delimiter=',')
factor_unhost_10 = []
factor_unhost_20 = []
index_in_S_unhost_quadruple = 0
for index in Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_25p5_and_4_images:
    index = int(index)
    print(index_in_S_unhost_quadruple)
    try:
        mcmc_new_list_unhost = np.loadtxt('./' + det + '_unhost_26/' + str(index) + 'corner_Res_data.csv', delimiter=',')
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
        
            

        factor_tmp_1_unhost = len(factor_tmp_1_index_unhost) / len(delta_t_ratio_10_20_unhost) * S_unhost_quadruple[index_in_S_unhost_quadruple] #* np.sum(SFR_unhost[factor_tmp_1_index_unhost]) / np.sum(SFR_unhost)
        factor_tmp_2_unhost = len(factor_tmp_2_index_unhost) / len(delta_t_ratio_10_30_unhost) * S_unhost_quadruple[index_in_S_unhost_quadruple] #* np.sum(SFR_unhost[factor_tmp_2_index_unhost]) / np.sum(SFR_unhost)
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
        mcmc_new_list_host = np.loadtxt('./' + det + '_host_pop/' + str(index_GW) + '_' + str(index_tmp) + 'corner_Res_data.csv', delimiter=',')
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
        
        
        factor_tmp_1_host = len(factor_tmp_1_index_host) / len(delta_t_ratio_10_20_host) * S_host_quadruple[index_tmp] #* np.sum(SFR_host[factor_tmp_1_index_host]) / np.sum(SFR_host)
        factor_tmp_2_host = len(factor_tmp_2_index_host) / len(delta_t_ratio_10_30_host) * S_host_quadruple[index_tmp] #* np.sum(SFR_host[factor_tmp_2_index_host]) / np.sum(SFR_host)
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
plt.hist(JWST_factor_host_sum_weighted, bins=1000, histtype="step", cumulative=-1, label='host', density='True')
plt.xlabel('Bayes factor')
plt.ylabel('PDF')
plt.legend()
plt.grid()
# plt.show()
# plt.ylim(0, 100)
plt.savefig('./IntermediaPlot/' + det + 'Bayes_statistc_long_time_Center' + SFR_weight_in_Galaxy + '.png', dpi=450)
plt.close()
np.savetxt('./IntermediaPlot/' + det + 'Bayes_factor_unhost_long_time_Center' + SFR_weight_in_Galaxy + '.csv', factor_unhost_sum, delimiter=',')
np.savetxt('./IntermediaPlot/' + det + 'Bayes_factor_host_long_time_Center' + SFR_weight_in_Galaxy + '.csv', factor_host_sum, delimiter=',')
np.savetxt('./IntermediaPlot/' + det + 'Bayes_factor_host_long_time_SFR_weighted_Center' + SFR_weight_in_Galaxy + '.csv', JWST_factor_host_sum_weighted, delimiter=',')






#计算误报率
FAP_peryear1sigma = []
FAP_peryear2sigma = []
FAP_peryear3sigma = []
year_obs = 11

factor_unhost_sum = np.loadtxt('./IntermediaPlot/' + det + 'Bayes_factor_unhost_long_time_Center' + SFR_weight_in_Galaxy + '.csv', delimiter=',')
JWST_factor_host_sum_weighted = np.loadtxt('./IntermediaPlot/' + det + 'Bayes_factor_host_long_time_SFR_weighted_Center' + SFR_weight_in_Galaxy + '.csv', delimiter=',')

mean_host = np.percentile(JWST_factor_host_sum_weighted, 50)
percentile_host = np.percentile(JWST_factor_host_sum_weighted, [15.87, 2.2, 0.13])
lower_limit1 = percentile_host[0]
lower_limit2 = percentile_host[1]
lower_limit3 = percentile_host[2]
FAP_peryear1sigma.append(np.sum(factor_unhost_sum > lower_limit1) / year_obs)
FAP_peryear2sigma.append(np.sum(factor_unhost_sum > lower_limit2) / year_obs)
FAP_peryear3sigma.append(np.sum(factor_unhost_sum > lower_limit3) / year_obs)

plt.scatter(0, FAP_peryear1sigma, label='1$\sigma$') 
plt.scatter(0, FAP_peryear2sigma, label='2$\sigma$') 
plt.scatter(0, FAP_peryear3sigma, label='3$\sigma$') 
plt.legend()
# plt.xlim(-1, 101)
# plt.xticks([])
plt.ylabel('FAR per year')
plt.grid()

plt.savefig('./IntermediaPlot/FAR_per_year_Center' + SFR_weight_in_Galaxy + '.png', dpi=450)
plt.close()