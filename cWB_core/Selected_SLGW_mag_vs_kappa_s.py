%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import corner
from getdist import plots, MCSamples
import getdist
import matplotlib.lines as mlines
import IPython



#NA revision1版本中能被认证出的事例58/510
#被match认证出的事例
with open('./PairResult/identified_file_name.csv', 'r') as f1:
    Selected_file_name = f1.readlines()

#NA revision2版本中能被认证出的事例91/516
#被match认证出的quadruple事例
with open('./PairResult_quadruple/identified_file_name.csv', 'r') as f1:
    Selected_file_name_quadruple = f1.readlines()


# #only macro
SNR_network = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')

four_image_index = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/four_image_index.csv', delimiter=',')

AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')

SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')

magnification = np.loadtxt('../Paper4_CE_Modify/SampleResult/magnification.csv', delimiter=',')

kappa = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('../Paper4_CE_Modify/SampleResult/gamma.csv', delimiter=',')
imagenum = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
kappa_s_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa_s.csv', delimiter=',')

color_set = ['indianred','royalblue','goldenrod', 'black', 'seagreen']
labels=['\log{\\rho}', '\log{|\mu|}', '\log{\kappa_{\\ast}}']
names = labels
corner_res = []
corner_res_total = []

#能在第一步认证的
for i_ in range(len(Selected_file_name)):
    tmp = Selected_file_name[i_].split('/')[3].split('_')
    i = int(tmp[0][6::])
    save_index = int(tmp[1])
    tmp_i_ = np.where(AcceptLensIndex == i)[0][0]
    if tmp_i_ >= 256:
        print("Error !!!!!!")
    else:
        sum_num = int(np.sum(imagenum[0:tmp_i_]))
        sum_num += save_index
        corner_res.append([np.log10(SNR_network[sum_num]), np.log10(abs(magnification[sum_num])), np.log10(kappa_s_set[sum_num])]) 

#第一步中所有的
for i_ in range(len(AcceptLensIndex[0:256])):
    i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
    imag_num = int(imagenum[i_])
    save_index = 0
    IndexSumimage = int(np.sum(imagenum[0:i_]))
    for j_ in range(imag_num):
        
        if SNR_network[IndexSumimage] > 12:
            
            corner_res_total.append([np.log10(SNR_network[IndexSumimage]), np.log10(abs(magnification[IndexSumimage])), np.log10(kappa_s_set[IndexSumimage])])
    
        IndexSumimage += 1
            

samples1 = MCSamples(samples=np.array(corner_res), names = names, labels = labels, label='Microlensing identifiable')

samples2 = MCSamples(samples=np.array(corner_res_total), names = names, labels = labels, label='Total events')


g = plots.get_subplot_plotter()
samples1.updateSettings({'contours': [0.6826, 0.99]})

# samples31.updateSettings({'contours': [0.6826, 0.9544]})

g.settings.legend_fontsize = 15
g.settings.axes_fontsize = 15
g.settings.axes_labelsize = 15
g.settings.axis_tick_x_rotation = 45
g.triangle_plot([samples1, samples2], colors=[color_set[0], color_set[1]],filled=False,
                line_args=[{'ls':'--','color':color_set[0]},{'ls':'--','color':color_set[1]}]
            ,marker_args={'color':'black','lw':2,'lw':1} ,contour_lws=[1.5])

plt.savefig('./Plot_statistic/statistic.png', dpi=450)