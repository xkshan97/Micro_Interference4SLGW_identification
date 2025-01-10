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

det = 'JWST'
#引力波的时间延迟
index_GW = 120
timedelay_GW_set = np.loadtxt('./Random_pin_down_GW/t_days.csv', delimiter=',')

for index_GW_random in range(len(timedelay_GW_set)):
    factor_unhost_sum = np.loadtxt('./IntermediaPlot/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_factor_unhost_long_time.csv', delimiter=',')
    factor_host_sum = np.loadtxt('./IntermediaPlot/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time.csv', delimiter=',')


if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = len(timedelay_GW_set)
    length = threadcount#int(threadcount * 150) #
    
    
    percount = length//threadcount

    print("percount = ", percount)

    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
    
        res = pool.apply_async(func=GW_TD_consistence_GG, args=(i,))
        
    pool.close()
    pool.join()
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