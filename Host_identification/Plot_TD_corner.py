# %matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import corner
from getdist import plots, MCSamples
import getdist
import matplotlib.lines as mlines
import IPython

#引力波的时间延迟
index_GW = 120
index_tmp = 25

timedelay_GW_set = np.loadtxt('./Random_pin_down_GW/t_days.csv', delimiter=',')



color_set = ['indianred','royalblue','goldenrod', 'black', 'seagreen']

# CSST_mcmc_new_list_host = np.loadtxt('./CSST_host_pop/' + str(index_GW) + '_' + str(index_tmp) + 'corner_Res_data.csv', delimiter=',')
JWST_mcmc_new_list_host = np.loadtxt('./JWST_host_pop/' + str(index_GW) + '_' + str(index_tmp) + 'corner_Res_data.csv', delimiter=',')
#labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
# CSST_t0_host = CSST_mcmc_new_list_host[:,5]
# CSST_t1_host = CSST_mcmc_new_list_host[:,6]
# CSST_t2_host = CSST_mcmc_new_list_host[:,7]
# CSST_t3_host = CSST_mcmc_new_list_host[:,8]

JWST_t0_host = JWST_mcmc_new_list_host[:,5]
JWST_t1_host = JWST_mcmc_new_list_host[:,6]
JWST_t2_host = JWST_mcmc_new_list_host[:,7]
JWST_t3_host = JWST_mcmc_new_list_host[:,8]
JWST_delta_t_10_host = JWST_t1_host - JWST_t0_host
JWST_delta_t_20_host = JWST_t2_host - JWST_t0_host
JWST_delta_t_30_host = JWST_t3_host - JWST_t0_host
JWST_delta_t_ratio_10_20_host = JWST_delta_t_10_host / JWST_delta_t_20_host
JWST_delta_t_ratio_10_30_host = JWST_delta_t_10_host / JWST_delta_t_30_host


GW_t0 = timedelay_GW_set[:,0]
GW_t1 = timedelay_GW_set[:,1]
GW_t2 = timedelay_GW_set[:,2]
GW_t3 = timedelay_GW_set[:,3]
GW_delta_t10 = GW_t1 - GW_t0
GW_delta_t20 = GW_t2 - GW_t0
GW_delta_t30 = GW_t3 - GW_t0
GW_delta_t_ratio_10_20 = GW_delta_t10 / GW_delta_t20
GW_delta_t_ratio_10_30 = GW_delta_t10 / GW_delta_t30


# CSST_corner= []
JWST_corner= []
GW_corner= []

# for i in range(len(CSST_t0_host)):
#     CSST_corner.append([CSST_t0_host[i], CSST_t1_host[i], CSST_t2_host[i], CSST_t3_host[i]])

# for i in range(len(JWST_delta_t_ratio_10_20_host)):
#     JWST_corner.append([JWST_delta_t_ratio_10_20_host[i], JWST_delta_t_ratio_10_30_host[i]])
    


for i in range(len(GW_delta_t10)):
    GW_corner.append([GW_delta_t10[i], GW_delta_t20[i], GW_delta_t30[i]])
    
for i in range(len(JWST_delta_t_10_host)):
    JWST_corner.append([JWST_delta_t_10_host[i], JWST_delta_t_20_host[i], JWST_delta_t_30_host[i]])
    

labels = ['t_{21}', 't_{31}', 't_{41}']
names = ['t_{21}', 't_{31}', 't_{41}']


# labels = ['\\frac{\Delta\;t_{21}}{\Delta\;t_{31}}', '\\frac{\Delta\;t_{21}}{\Delta\;t_{41}}']
# names = ['\\frac{\Delta\;t_{21}}{\Delta\;t_{31}}', '\\frac{\Delta\;t_{21}}{\Delta\;t_{41}}']

# samples1 = MCSamples(samples=np.array(CSST_corner), names = names, labels = labels, label='CSST')
# samples2 = MCSamples(samples=np.array(JWST_corner), names = names, labels = labels, label='JWST')


# g = plots.get_subplot_plotter()
# samples1.updateSettings({'contours': [0.6826, 0.9544]})
# samples2.updateSettings({'contours': [0.6826, 0.9544]})


# g.settings.legend_fontsize = 10
# g.settings.axes_fontsize = 15
# g.settings.axes_labelsize = 15
# g.settings.axis_tick_x_rotation = 45
# g.triangle_plot([samples1, samples2], colors=[color_set[0], color_set[1]],filled=True,legend_loc='upper right',
#                 line_args=[{'color':color_set[0]},{'color':color_set[1]}], markers={'t_0':td_GW_true[0], 't_1':td_GW_true[1], 't_2':td_GW_true[2], 't_3':td_GW_true[3]}
#             ,marker_args={'color':'black','lw':20,'lw':1} ,contour_lws=[1.5,1.5], shaded=True)

# plt.savefig('./IntermediaPlot/TD_corner.png', dpi=450)



samples1 = MCSamples(samples=np.array(GW_corner), names = names, labels = labels, label='GW')

samples2 = MCSamples(samples=np.array(JWST_corner), names = names, labels = labels, label='Host Galaxy')


g = plots.get_subplot_plotter()
samples1.updateSettings({'contours': [0.6826, 0.995]})
samples2.updateSettings({'contours': [0.6826, 0.995]})


g.settings.legend_fontsize = 10
g.settings.axes_fontsize = 15
g.settings.axes_labelsize = 15
g.settings.axis_tick_x_rotation = 45
g.triangle_plot([samples1, samples2], colors=[color_set[0], color_set[1]],filled=True,legend_loc='upper right',
                line_args=[{'color':color_set[0]}, {'color':color_set[1]}]
            ,marker_args={'color':'grey','lw':2,'lw':1} ,contour_lws=[1.5])

plt.savefig('./IntermediaPlot/TD_corner_' + str(index_tmp) + '_.png', dpi=450)