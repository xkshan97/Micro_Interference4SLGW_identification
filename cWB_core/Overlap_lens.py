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
import corner
from getdist import plots, MCSamples
import getdist
import matplotlib.lines as mlines
import IPython
%matplotlib inline




SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')

with open('./PairResult/identified_index.csv', 'r') as f0:
    Selected_index = f0.readlines()
    
#被match认证出的事例
with open('./PairResult/identified_file_name_not_duplication.csv', 'r') as f1:
    Selected_file_name_not_duplication = f1.readlines()

#with micro
filesminimum = os.listdir('../Paper4_CE_Modify/outdir_Lensed_Minimum_with_micro')
filesminimum_index = [] #储存不带save_index的index
filessaddle = os.listdir('../Paper4_CE_Modify/outdir_Lensed_Saddle_with_micro') 
filessaddle_index = []
#获得参数估计文件中的文件index
for i in range(len(filesminimum)):
    tmp_index = filesminimum[i].split('r')[1]
    filesminimum_index.append(tmp_index.split('_')[0])
    
for i in range(len(filessaddle)):
    tmp_index = filessaddle[i].split('r')[1]
    filessaddle_index.append(tmp_index.split('_')[0])

combine_num = 0
Selected_file_name_pair = []
for i in range(len(Selected_index)):
    Selected_file_name_pair.append([Selected_file_name_not_duplication[i].split('\n')[0]])
    tmp_index = Selected_index[i].split('\n')[0]
    index1 = np.where(np.array(filesminimum_index) == tmp_index)[0]
    index2 = np.where(np.array(filessaddle_index) == tmp_index)[0]
    if len(index1 != 0): 
        for j in range(len(index1)):
            file_tmp1 = '../Paper4_CE_Modify/outdir_Lensed_Minimum_with_micro/' + filesminimum[index1[j]] + '/H1_L1_V1_result.json'
            if file_tmp1 not in Selected_file_name_pair[i]:
                Selected_file_name_pair[i].append(file_tmp1)
            else:
                pass
    if len(index2 != 0):
        for j in range(len(index2)):
            file_tmp2 = '../Paper4_CE_Modify/outdir_Lensed_Saddle_with_micro/' + filessaddle[index2[j]] + '/H1_L1_V1_result.json'
            if file_tmp2 not in Selected_file_name_pair[i]:
                Selected_file_name_pair[i].append(file_tmp2)
            else:
                pass 
            
                
            


#2参数
labels = ['\mathrm{ra}', '\mathrm{\sin{(dec)}}']
names = ['\mathrm{ra}', '\mathrm{\sin{(dec)}}']
BBHParam = ['ra', 'dec']
character_Sky_RA_map_error = []
character_Sky_DEC_map_error = []

#画一下两参数的corner图，然后把数据读进Poster中。
for i in tqdm(range(len(Selected_file_name_pair))):
    Poster = [[]] * len(Selected_file_name_pair[i])
    samples_sum = [[]] * len(Selected_file_name_pair[i])
    index_in_SampleParam = int(Selected_file_name_pair[i][0].split('outdir')[-1].split('_')[0])
    
    for j in range(len(Selected_file_name_pair[i])):
        Poster[j] = [[]] * len(BBHParam)
        filej = Selected_file_name_pair[i][j]
        with open(filej) as f2:
            dataj = json.load(f2)
        nest_samplesj = dataj['posterior']['content']
        
        for k in range(len(Poster[j])):
            Poster[j][k] = nest_samplesj[BBHParam[k]]  
        if j == 0:
            character_Sky_RA_map_error.append(np.std(nest_samplesj[BBHParam[0]]))
            character_Sky_DEC_map_error.append(np.std(nest_samplesj[BBHParam[1]]))
        Poster[j] = np.array(Poster[j])
        Poster[j][1] = np.sin(Poster[j][1])    
        
        Plot_corner = []
        for _ in range(len(Poster[j][0])):
            Plot_corner.append([Poster[j][0][_], Poster[j][1][_]])
        Plot_corner = np.array(Plot_corner)
        
        

        samples_sum[j] = MCSamples(samples=np.array(Plot_corner), names = names, labels = labels, label = 'GW signal ' + str(j))
        
        
        samples_sum[j].updateSettings({'contours': [0.6826, 0.9544]})

    g = plots.get_subplot_plotter()
    g.settings.legend_fontsize = 10
    g.settings.axes_fontsize = 15
    g.settings.axes_labelsize = 15
    g.settings.axis_tick_x_rotation = 45
    g.triangle_plot(samples_sum, colors=['grey', 'grey', 'grey', 'grey'],filled=True,legend_loc='upper right',
            line_args=[{'color':'grey'},{'color':'grey'},{'color':'grey'},{'color':'grey'}], markers={'\mathrm{ra}':SampleParam[7][index_in_SampleParam], '\mathrm{\sin{(dec)}}':np.sin(SampleParam[8][index_in_SampleParam])}
        ,marker_args={'color':'black','lw':20,'lw':1} ,contour_lws=[1,1,1,1])
        # corner.corner(Plot_corner, fig = figure, labels=["ra", "sin(dec)"],
        #         show_titles=True, title_kwargs={"fontsize": 12},smooth=False, title_fmt='.4f')
        
        
        
    plt.savefig('./IntermediaPlot_corner_Pair/' + Selected_index[i].split('\n')[0] + '.png', dpi=450)
    plt.close()
    
    #以第一个像为基准，计算Overlap，因为第一个像是被认证出的。
    Bayes_factor = []
    for j in range(len(Selected_file_name_pair[i]) -1):
        minset, maxset = [], []
        for ii in range(len(Poster[0])):
            minset.append(np.max([Poster[0][ii].min(), Poster[j+1][ii].min()]))
            maxset.append(np.min([Poster[0][ii].max(), Poster[j+1][ii].max()]))
        ra_, dec_ = np.mgrid[minset[0]:maxset[0]:100j, minset[1]:maxset[1]:100j]

        positions = np.vstack([ra_.ravel(), dec_.ravel()])

        
        kernel1 = sp.stats.gaussian_kde(Poster[0])
        kernel2 = sp.stats.gaussian_kde(Poster[j+1])
        # kernelPrior = sp.stats.gaussian_kde(Prior1)

        
        z1 = kernel1(positions)
        z2 = kernel2(positions)

        Bayse_ = np.sum(np.array(z1)*np.array(z2)) * np.abs(ra_[1][0] - ra_[0][0]) * np.abs(dec_[0][1] - dec_[0][0]) \
             * 2 * 2 * np.pi
            
        Bayes_factor.append(Bayse_)
    np.savetxt('./Bayes_Lens/' + Selected_index[i].split('\n')[0] + '.csv', Bayes_factor, delimiter=',')
        
        
Bayes_factor_sum = []
for i in range(len(Selected_file_name_pair)):
    try:
        Bayes_factor = [float(np.loadtxt('./Bayes_Lens/' + Selected_index[i].split('\n')[0] + '.csv', delimiter=','))]
    except TypeError:
        Bayes_factor = np.loadtxt('./Bayes_Lens/' + Selected_index[i].split('\n')[0] + '.csv', delimiter=',')
    for j in range(len(Bayes_factor)):
        Bayes_factor_sum.append(Bayes_factor[j])
        
plt.hist(Bayes_factor_sum, bins=100)