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

# H1_snr = np.loadtxt('/data/Shanxikai/Paper4/Sim_GW_Data/H1/H1_snr.csv', delimiter=',')
# L1_snr = np.loadtxt('/data/Shanxikai/Paper4/Sim_GW_Data/L1/L1_snr.csv', delimiter=',')
# V1_snr = np.loadtxt('/data/Shanxikai/Paper4/Sim_GW_Data/V1/V1_snr.csv', delimiter=',')
# K1_snr = np.loadtxt('/data/Shanxikai/Paper4/Sim_GW_Data/K1/K1_snr.csv', delimiter=',')
# SNR_network = np.sqrt(H1_snr**2 + L1_snr**2 + V1_snr**2)
# Index_SNR_network_G_12 = np.where(SNR_network >= 12)[0]

files = os.listdir('../Paper4_CE_Modify/outdir_unlensed')
files_minimum = os.listdir('../Paper4_CE_Modify/outdir_Lensed_Minimum_with_micro') 
files_saddle = os.listdir('../Paper4_CE_Modify/outdir_Lensed_Saddle_with_micro') 


# #4参数
# BBHParam = ['mass_ratio', 'chirp_mass', 'ra', 'dec']
#3参数
BBHParam = ['ra', 'dec']

files_combine = []
for i in range(len(files_minimum)):
    for j in range(len(files)):
        files_minimum_tmp = files_minimum[i]
        lens_index = files_minimum_tmp.split('_')[0][6::]
        files_tmp = files[j]
        event_index =  files_tmp[6::]
        if lens_index != event_index:
            files_combine.append(['../Paper4_CE_Modify/outdir_Lensed_Minimum_with_micro/' + files_minimum_tmp + '/H1_L1_V1_result.json', '../Paper4_CE_Modify/outdir_unlensed/' + files_tmp + '/H1_L1_V1_result.json'])

for i in range(len(files_saddle)):
    for j in range(len(files)):
        files_saddle_tmp = files_saddle[i]
        lens_index = files_saddle_tmp.split('_')[0][6::]
        files_tmp = files[j]
        event_index =  files_tmp[6::]
        if lens_index != event_index:
            files_combine.append(['../Paper4_CE_Modify/outdir_Lensed_Saddle_with_micro/' + files_saddle_tmp + '/H1_L1_V1_result.json', '../Paper4_CE_Modify/outdir_unlensed/' + files_tmp + '/H1_L1_V1_result.json'])
    
   

        
def OverLap(Index1, Index2):
    Bayse = []
    for i in tqdm(range(Index1, Index2)):
        Poster1 = [[]]*len(BBHParam)
        file1_ = files_combine[i][0]
        file2_ = files_combine[i][1]
        # file1_ = '../Paper4_CE_Modify/outdir_unlensed/' + file1 + '/H1_L1_V1_result.json'
        try:
            with open(file1_) as f1:
                data1 = json.load(f1)
        except FileNotFoundError:
            print(file1_)
            continue
        nest_samples1 = data1['posterior']['content']
        for ii in range(len(BBHParam)):
            Poster1[ii] = nest_samples1[BBHParam[ii]]  
        Poster1 = np.array(Poster1)
        Poster1[1] = np.sin(Poster1[1])
        #=================上面是事例1的
        Poster2 = [[]]*len(BBHParam)
        # file2_ = '../Paper4_CE_Modify/outdir_unlensed/' + file2 + '/H1_L1_V1_result.json'
        try:
            with open(file2_) as f2:
                data2 = json.load(f2)
        except FileNotFoundError:
            print(file2_)
            continue
            
        nest_samples2 = data2['posterior']['content']
        for jj in range(len(BBHParam)):
            Poster2[jj] = nest_samples2[BBHParam[jj]] 
        Poster2 = np.array(Poster2)
        Poster2[1] = np.sin(Poster2[1])

        minset, maxset = [], []
        for ii in range(len(Poster1)):
            minset.append(np.max([Poster1[ii].min(), Poster2[ii].min()]))
            maxset.append(np.min([Poster1[ii].max(), Poster2[ii].max()]))
            


        
        #3参数
        ra_, dec_ = np.mgrid[minset[0]:maxset[0]:100j, minset[1]:maxset[1]:100j]

        positions = np.vstack([ra_.ravel(), dec_.ravel()])




        kernel1 = sp.stats.gaussian_kde(Poster1)
        kernel2 = sp.stats.gaussian_kde(Poster2)
        # kernelPrior = sp.stats.gaussian_kde(Prior1)

        
        z1 = kernel1(positions)
        z2 = kernel2(positions)

        Bayse_ = np.sum(np.array(z1)*np.array(z2)) * np.abs(ra_[1][0] - ra_[0][0]) * np.abs(dec_[0][1] - dec_[0][0]) \
            * 2 * 2 * np.pi
            
        Bayse.append(Bayse_)

    return Bayse
    # return 0








if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = 362
    res_z = [0 for i in range(threadcount)]
    length = len(files_combine)
   
     
    percount = length//threadcount

    print("percount = ", percount)
   
    t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
        if i < threadcount - 1:
            res = pool.apply_async(func=OverLap, args=(percount * i, percount * (i + 1),))
            res_z[i] = res
            # print("0",i)
        else:
            # print("1",i)
            res = pool.apply_async(func=OverLap, args=(percount * i, length,)) 
            res_z[i] = res

    pool.close()
    pool.join() 
    Result1 = list(res_z[0].get())
    # Result1, Result2 = list(res_z[0].get()[0]), list(res_z[0].get()[1])
    for i in range(threadcount-1):
        Result1.extend(list(res_z[i+1].get()))
      
    
    print(time.time()-t0)
    
    np.savetxt('./Bayes_Unlens/overlap_LU.csv', Result1, delimiter=',')
  
