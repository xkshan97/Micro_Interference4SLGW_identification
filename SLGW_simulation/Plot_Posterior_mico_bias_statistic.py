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

    

AcceptLensIndex = np.loadtxt('./SampleResult/AcceptLensIndex.csv', delimiter=',')
SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')

with open('../Paper5_cWB_Modify/PairResult/identified_index.csv', 'r') as f0:
    Selected_index = f0.readlines()
    
#被match认证出的事例
with open('../Paper5_cWB_Modify/PairResult/identified_file_name_not_duplication.csv', 'r') as f1:
    Selected_file_name_not_duplication = f1.readlines()
    
  
#with micro
filesminimum = os.listdir('./outdir_Lensed_Minimum_with_micro')
filesminimum_index = [] #储存不带save_index的index
filessaddle = os.listdir('./outdir_Lensed_Saddle_with_micro') 
filessaddle_index = []
#获得参数估计文件中的文件index
for i in range(len(filesminimum)):
    tmp_index = filesminimum[i].split('r')[1]
    filesminimum_index.append(tmp_index.split('_')[0])
    
for i in range(len(filessaddle)):
    tmp_index = filessaddle[i].split('r')[1]
    filessaddle_index.append(tmp_index.split('_')[0])

sum_num = 0
Selected_file_sum = []
for i in range(len(Selected_index)):
    Selected_file_sum.append(Selected_file_name_not_duplication[i].split('\n')[0])
    tmp_index = Selected_index[i].split('\n')[0]
    index1 = np.where(np.array(filesminimum_index) == tmp_index)[0]
    index2 = np.where(np.array(filessaddle_index) == tmp_index)[0]
    if len(index1 != 0): 
        for j in range(len(index1)):
            file_tmp1 = './outdir_Lensed_Minimum_with_micro/' + filesminimum[index1[j]] + '/H1_L1_V1_result.json'
            if file_tmp1 not in Selected_file_sum:
                Selected_file_sum.append(file_tmp1)
            else:
                pass
    if len(index2 != 0):
        for j in range(len(index2)):
            file_tmp2 = './outdir_Lensed_Saddle_with_micro/' + filessaddle[index2[j]] + '/H1_L1_V1_result.json'
            if file_tmp2 not in Selected_file_sum:
                Selected_file_sum.append(file_tmp2)
            else:
                pass 
            
       

BBHParam = ['mass_ratio', 'chirp_mass', 'ra', 'dec', 'theta_jn', 'psi', 'a_1', 'a_2']



def Cal_Rel_bias(inject_value, data_micro):
    
    Qrec = [[], [], [], [], [], [], [], []]
    nest_samples_micro = data_micro['posterior']['content']
    
    for ii in range(len(BBHParam)):
        inject_tmp = inject_value[BBHParam[ii]]
        mean_micro = np.mean(nest_samples_micro[BBHParam[ii]])
        std_micro = np.std(nest_samples_micro[BBHParam[ii]])
        tmp = np.abs(inject_tmp - mean_micro) / std_micro
        
        Qrec[ii]=tmp
    
    return Qrec



Qrec_res = []

for i in tqdm(range(len(Selected_file_sum))):
    image_index = int(Selected_file_sum[i].split('/')[-2].split('_')[0].split('r')[-1])
    save_index = Selected_file_sum[i].split('/')[-2].split('_')[-1]
    file_micro_tmp = Selected_file_sum[i]
    

    with open(file_micro_tmp) as f_micro:
        data_micro = json.load(f_micro)
        
    Mz = (SampleParam[0][image_index] + 1) * (SampleParam[1][image_index] * SampleParam[2][image_index]) ** (3/5) / (SampleParam[1][image_index] + SampleParam[2][image_index]) ** (1/5)
    injection_parameters = dict(
    mass_ratio = SampleParam[2][image_index]/SampleParam[1][image_index], chirp_mass = Mz,  a_1=SampleParam[3][image_index], 
    a_2=SampleParam[4][image_index],
    theta_jn=SampleParam[5][image_index], psi=SampleParam[6][image_index],
    ra=SampleParam[7][image_index], dec=SampleParam[8][image_index])
    
    Qrec_res.append(Cal_Rel_bias(injection_parameters, data_micro))



np.savetxt('Final_Bias_show/Bias_Paper5_Selected.csv', Qrec_res, delimiter=',')
