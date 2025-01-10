'''
copy的Paper4_copy_CE中的File4Plot_final_modify.py。
画一下Paper5_cWB_Modify中能挑出来事例的参数bias，表明只有ra和dec的bias小。
'''
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
from tqdm import tqdm
import scipy as sp
from sklearn.neighbors import KernelDensity
import time
from  multiprocessing import Process,Pool
import os
import json
from getdist import plots, MCSamples
import getdist

    

# AcceptLensIndex = np.loadtxt('./SampleResult/AcceptLensIndex.csv', delimiter=',')
# SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')

# with open('../Paper5_cWB_Modify/PairResult/identified_index.csv', 'r') as f0:
#     Selected_index = f0.readlines()
    
# #被match认证出的事例
# with open('../Paper5_cWB_Modify/PairResult/identified_file_name_not_duplication.csv', 'r') as f1:
#     Selected_file_name_not_duplication = f1.readlines()
    
  
# #with micro
# filesminimum = os.listdir('./outdir_Lensed_Minimum_with_micro')
# filesminimum_index = [] #储存不带save_index的index
# filessaddle = os.listdir('./outdir_Lensed_Saddle_with_micro') 
# filessaddle_index = []
# #获得参数估计文件中的文件index
# for i in range(len(filesminimum)):
#     tmp_index = filesminimum[i].split('r')[1]
#     filesminimum_index.append(tmp_index.split('_')[0])
    
# for i in range(len(filessaddle)):
#     tmp_index = filessaddle[i].split('r')[1]
#     filessaddle_index.append(tmp_index.split('_')[0])

# sum_num = 0
# Selected_file_sum = []
# for i in range(len(Selected_index)):
#     Selected_file_sum.append(Selected_file_name_not_duplication[i].split('\n')[0])
#     tmp_index = Selected_index[i].split('\n')[0]
#     index1 = np.where(np.array(filesminimum_index) == tmp_index)[0]
#     index2 = np.where(np.array(filessaddle_index) == tmp_index)[0]
#     if len(index1 != 0): 
#         for j in range(len(index1)):
#             file_tmp1 = './outdir_Lensed_Minimum_with_micro/' + filesminimum[index1[j]] + '/H1_L1_V1_result.json'
#             if file_tmp1 not in Selected_file_sum:
#                 Selected_file_sum.append(file_tmp1)
#             else:
#                 pass
#     if len(index2 != 0):
#         for j in range(len(index2)):
#             file_tmp2 = './outdir_Lensed_Saddle_with_micro/' + filessaddle[index2[j]] + '/H1_L1_V1_result.json'
#             if file_tmp2 not in Selected_file_sum:
#                 Selected_file_sum.append(file_tmp2)
#             else:
#                 pass 

Micro_file_min = os.listdir('outdir_Lensed_Minimum_with_micro')
for i in range(len(Micro_file_min)):
    Micro_file_min[i] = './outdir_Lensed_Minimum_with_micro/' + Micro_file_min[i]  + '/H1_L1_V1_result.json'
Micro_file_sad = os.listdir('outdir_Lensed_Saddle_with_micro')
for i in range(len(Micro_file_sad)):
    Micro_file_sad[i] = './outdir_Lensed_Saddle_with_micro/' + Micro_file_sad[i]  + '/H1_L1_V1_result.json'
    
Micro_file_micro = Micro_file_min + Micro_file_sad   
       
Micro_file_unlens = os.listdir('outdir_unlensed')
for i in range(len(Micro_file_unlens)):
    Micro_file_unlens[i] = './outdir_unlensed/' + Micro_file_unlens[i] + '/H1_L1_V1_result.json'
 
corner_micro_pro = []
corner_unlens_pro = []

for i in tqdm(range(len(Micro_file_micro))):
    file_micro_tmp = Micro_file_micro[i]
    
    with open(file_micro_tmp) as f_micro:
        data_micro = json.load(f_micro)
    j = np.where(data_micro['posterior']['content']['log_likelihood'] == np.max(data_micro['posterior']['content']['log_likelihood']))[0][0]
    spin_1x = data_micro['posterior']['content']['spin_1x'][j]
    spin_1y = data_micro['posterior']['content']['spin_1y'][j]
    spin_1z = data_micro['posterior']['content']['spin_1z'][j]
    spin_2x = data_micro['posterior']['content']['spin_2x'][j]
    spin_2y = data_micro['posterior']['content']['spin_2y'][j]
    spin_2z = data_micro['posterior']['content']['spin_2z'][j]
    
    corner_micro_pro.append([spin_1x, spin_1y, spin_2x, spin_2y])
    

for i in tqdm(range(len(Micro_file_unlens))):
    file_unlens_tmp = Micro_file_unlens[i]
    
    with open(file_unlens_tmp) as f_unlens:
        data_unlens = json.load(f_unlens)
    j = np.where(data_unlens['posterior']['content']['log_likelihood'] == np.max(data_unlens['posterior']['content']['log_likelihood']))[0][0]
    spin_1x = data_unlens['posterior']['content']['spin_1x'][j]
    spin_1y = data_unlens['posterior']['content']['spin_1y'][j]
    spin_1z = data_unlens['posterior']['content']['spin_1z'][j]
    spin_2x = data_unlens['posterior']['content']['spin_2x'][j]
    spin_2y = data_unlens['posterior']['content']['spin_2y'][j]
    spin_2z = data_unlens['posterior']['content']['spin_2z'][j]
    
    corner_unlens_pro.append([spin_1x, spin_1y, spin_2x, spin_2y])
    

labels = ['|s_{1x}|', '|s_{1y}|', '|s_{2x}|', '|s_{2y}|']
names = ['|s_{1x}|', '|s_{1y}|', '|s_{2x}|', '|s_{2y}|']


samples1 = MCSamples(samples=np.abs(np.array(corner_unlens_pro)), names = names, labels = labels, label='Unlens')
samples2 = MCSamples(samples=np.abs(np.array(corner_micro_pro)), names = names, labels = labels, label='Micro')


g = plots.get_subplot_plotter()
samples1.updateSettings({'contours': [0.6826, 0.995]})
samples2.updateSettings({'contours': [0.6826, 0.995]})
# g.add_legend(['Unlens ', 'Micro'])

g.settings.legend_fontsize = 20
g.settings.axes_fontsize = 15
g.settings.axes_labelsize = 15
g.settings.axis_tick_x_rotation = 45
g.triangle_plot([samples1, samples2], colors=['grey', 'indianred'],filled=True,legend_loc='upper right',
                line_args=[{'color':'grey'}, {'color':'indianred'}] ,contour_lws=[1.5])

plt.savefig('./IntermediaPlot/corner_processing.png', dpi=450)
