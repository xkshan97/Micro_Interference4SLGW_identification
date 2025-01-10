#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 21:29:33 2022

@author: astroshan
"""
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy as sp
import matplotlib.lines as mlines
import seaborn as sns
from tqdm import tqdm
import random
from astropy.modeling.models import Sersic2D
import struct
from scipy.interpolate import interp1d
import json
#plt.style.use('science')
color_set = ['indianred','royalblue','goldenrod', 'black', 'seagreen']


'''图1'''
Mz_cut_set = [0, 20]
fig, ax = plt.subplots(2, 1, figsize=(6.5, 5))
#gs1 = plt.GridSpec(2, 1)
#gs.update(wspace=0.4)
#plt.figure(figsize=(6.5,6.5))
for i in range(len(Mz_cut_set)):
    Mz_cut = Mz_cut_set[i]
    low_SNR_envelope = np.loadtxt('../../data/Result_file4plot/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    low_envelope = np.loadtxt('../../data/Result_file4plot/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_SNR_envelope = np.loadtxt('../../data/Result_file4plot/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_envelope = np.loadtxt('../../data/Result_file4plot/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
        
    SNR_tmp_lensed_sum = np.loadtxt('../../data/Result_file4plot/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_mean_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_mean_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    theory_match_sum_tmp = np.loadtxt('../../data/Result_file4plot/theory_match' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_std_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    
#    
#    SNR_tmp_lensed04 = np.loadtxt('../../data/Result_file4plot/SNR_lensed_04' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    match_mean_lensed04 = np.loadtxt('../../data/Result_file4plot/match_mean_lensed_04' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    theory_match04_tmp = np.loadtxt('../../data/Result_file4plot/theory_match_04' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    match_std_lensed04 = np.loadtxt('../../data/Result_file4plot/match_std_lensed_04' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    
#    
#    SNR_tmp_lensed06 = np.loadtxt('../../data/Result_file4plot/SNR_lensed_06' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    match_mean_lensed06 = np.loadtxt('../../data/Result_file4plot/match_mean_lensed_06' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    theory_match06_tmp = np.loadtxt('../../data/Result_file4plot/theory_match_06' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    match_std_lensed06 = np.loadtxt('../../data/Result_file4plot/match_std_lensed_06' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    
#    
#    SNR_tmp_lensed08 = np.loadtxt('../../data/Result_file4plot/SNR_lensed_08' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    match_mean_lensed08 = np.loadtxt('../../data/Result_file4plot/match_mean_lensed_08' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    theory_match08_tmp = np.loadtxt('../../data/Result_file4plot/theory_match_08' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    match_std_lensed08 = np.loadtxt('../../data/Result_file4plot/match_std_lensed_08' + str(int(Mz_cut)) + '.csv', delimiter=',')

    
    
#    ax1 = plt.subplot(gs1[i,0])
    if i == 0:
        ax[0].fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='lightgrey', label='Total Unlens')
         #0.090
#        ax0_04 = ax[0].scatter(SNR_tmp_lensed04, match_mean_lensed04, c=theory_match04_tmp, s = 20, marker='v', cmap='Spectral')
#        ax0_06 = ax[0].scatter(SNR_tmp_lensed06, match_mean_lensed06, c=theory_match06_tmp, s = 20, marker='p', cmap='Spectral')
#        ax0_08 = ax[0].scatter(SNR_tmp_lensed08, match_mean_lensed08, c=theory_match08_tmp, s = 20, marker='x', cmap='Spectral')
#        
        ax0 = ax[0].scatter(SNR_tmp_lensed_sum, match_mean_lensed_sum, c='#0C5DA5', s = 10, marker='p')
        ax[0].errorbar(SNR_tmp_lensed_sum, match_mean_lensed_sum, yerr = 1.5*np.array(match_std_lensed_sum), elinewidth=1, ecolor='black', capsize=3, capthick=1, linestyle='none')
        
    else:
        ax[1].fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='lightgrey', label='Chirp Mass $\geq$' + str(int(Mz_cut)) + '$\mathrm{M}_\odot$')
         #0.090
        ax1 = ax[1].scatter(SNR_tmp_lensed_sum, match_mean_lensed_sum, c='#0C5DA5', s= 10, marker='p')
        ax[1].errorbar(SNR_tmp_lensed_sum, match_mean_lensed_sum, yerr = 1.5*np.array(match_std_lensed_sum), elinewidth=1, ecolor='black', capsize=3, capthick=1, linestyle='none')
        
       
    # plt.errorbar(SNR_tmp_unlens, match_mean_unlens, yerr = np.array(match_std_unlens), elinewidth=2, ecolor='black', capsize=5, capthick=1, linestyle='none')
    
   
    
    if i == 1:
        ax[1].set_xlabel('SNR')
    if i == 0:
        plt.setp(ax[0].get_xticklabels(), visible=False)
    ax[i].set_ylabel('Match')
    ax[i].grid()
    ax[i].set_ylim(0.935,1)
    ax[i].set_xlim(35,500)
#    ax1.set_yticks([0.95, 0.96, 0.97, 0.98, 0.99])
    if i == 1:
        yticks = ax[1].yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
#    ax1.set_xscale('log')
#    ax1.set_yscale('log')
    ax[i].legend()
plt.subplots_adjust(hspace=.0)
#fig.colorbar(ax0, ax=ax[0], fraction=0.017, label='Theoretical Match')
#fig.colorbar(ax1, ax=ax[1], fraction=0.017, label='Theoretical Match')
plt.savefig('../../results/match.png', dpi=450)
plt.close()
#    ax1.set_xticks(np.arange(50, 350, 50))
#    ax1.set_yticks(np.arange(0.9,1,0.02))
    


'''图二'''
freq_plot = np.loadtxt('../../data/Result_file4plot/freq_plot.csv', delimiter=',')
cWB_rec = np.loadtxt('../../data/Result_file4plot/cWB_hf.csv', delimiter=',')
inject_hf = np.loadtxt('../../data/Result_file4plot/inject_hf.csv', delimiter=',')
Bilby_rec = np.loadtxt('../../data/Result_file4plot/Bilby_hf.csv', delimiter=',')
Ffabs = np.loadtxt('../../data/Result_file4plot/Ffabs.csv', delimiter=',')


plt.figure(figsize=(6.5,5))
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogy(freq_plot, cWB_rec, label = 'cWB')
main_ax.semilogy(freq_plot, inject_hf, 'k', label = 'Injected')
main_ax.semilogy(freq_plot, Bilby_rec, 'r', label = 'Bilby')
# plt.xlim(40,250)
main_ax.set_ylim(10**(-25.3),10**(-23.5))
#    ax[_].set_xlim(40,300)
main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('Strain/Hz')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.plot(freq_plot, Ffabs / 8.139, 'k', label = '$F(f)/\sqrt{\mu}$')
sub_ax.plot(freq_plot, cWB_rec/Bilby_rec, color='#0C5DA5', label = 'cWB/Bilby')
sub_ax.grid()
sub_ax.legend()
sub_ax.set_ylim(-0.1,2.5)
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('Fraction')
plt.savefig('../../results/rec.png', dpi=450)
plt.close()


'''图三，pair认证结果'''
with open('../../data/PairResult/identified_index.csv', 'r') as f0:
    Selected_index = f0.readlines()  

#画图，x axis是FAP per pair，纵坐标是index
FAP_per_pair_sum = []
trial_factor = []
fig = plt.figure(figsize=(5,15))
ax2 = fig.add_subplot(111)
ax1 = ax2.twiny()
for i in range(len(Selected_index)):
    try:
        FAP_per_pair = [float(np.loadtxt('../../data/Bayes_Lens/FAP_per_pair_' + Selected_index[i].split('\n')[0] + '.csv', delimiter=','))]
    except TypeError:
        FAP_per_pair = np.loadtxt('../../data/Bayes_Lens/FAP_per_pair_' + Selected_index[i].split('\n')[0] + '.csv', delimiter=',')
    for j in range(len(FAP_per_pair)):
        FAP_per_pair_sum.append(FAP_per_pair[j])
        trial_factor.append(1 - (1 - FAP_per_pair[j])**100000)
        ax1.scatter(FAP_per_pair[j], i, color='grey', alpha=0)
        if FAP_per_pair[j] * 100000 < 1:
            ax2.scatter(FAP_per_pair[j] * 100000, i, color='indianred', marker='*', s=100)
        else:
            ax2.scatter(FAP_per_pair[j] * 100000, i, color='grey', s=20)
            
ax1.semilogx()
ax2.semilogx()
ax1.set_xlim(10**(-11), 1)
ax2.set_xlim(10**(-11) * 100000, 1 * 100000)
ax2.fill_between([1, 100000], y1=-1, y2=len(Selected_index), color='lightgrey', alpha=0.5)
ax1.set_ylim(-0.5, len(Selected_index) - 0.5)
ax2.set_ylim(-0.5, len(Selected_index) - 0.5)
ax1.set_yticks(range(len(Selected_index)))
ax2.set_yticks(range(len(Selected_index)))
ax1.set_xlabel('FAP')
ax2.set_xlabel('Number of false pair per year')
ax2.set_ylabel('Event index')
ax2.set_xticks([0.1, 10**(0), 10**(1), 10**(2), 10**(3), 10**(4), 10**5])
ax2.set_xticklabels(['$\leq0.1$', '$1$', '$10$', '$100$', '$1000$', '$10^4$', '$10^5$']) # 设置刻度标签
ax1.set_xticks([10**(-6), 10**(-5)])
ax1.set_xticklabels(['$0.09$', '$0.63\leq$'], ha='left') # 设置刻度标签
ax2.grid()
plt.savefig('../../results/identified_pair.png', dpi=450)
plt.close()








'''图4，Host galaxy的星等分布'''
index_set_selected = [46, 120, 170]
title_index = [27, 35, 42]
zs_plot = [3.9, 1.6, 5.9]
fig, ax = plt.subplots(3, 1, figsize=(6.5, 10))
plt.subplots_adjust(hspace=0.25)
ax_index = 0
for index in index_set_selected:
    mag_range = np.loadtxt('../../data/SampleResult_GWGalaxy/mag_range_' + str(index) + '.csv', delimiter=',') 
    mag_weight_by_SFR_50 = np.loadtxt('../../data/SampleResult_GWGalaxy/mag_range_weight_by_SFR_50' + str(index) + '.csv', delimiter=',')
    
    JWST_mag_range = np.loadtxt('../../data/SampleResult_GWGalaxy/JWST_mag_range_' + str(index) + '.csv', delimiter=',') 
    JWST_mag_weight_by_SFR_50 = np.loadtxt('../../data/SampleResult_GWGalaxy/JWST_mag_range_weight_by_SFR_50' + str(index) + '.csv', delimiter=',')
    
    mag_let_26_index = np.where(mag_weight_by_SFR_50 < 26)[0]
    mag_gtr_26_index = np.where(mag_weight_by_SFR_50 >= 26)[0]
    
    mag_weight_by_SFR_50_let_26 = mag_weight_by_SFR_50[mag_let_26_index]
    mag_weight_by_SFR_50_gtr_26 = mag_weight_by_SFR_50[mag_gtr_26_index]
    
    JWST_mag_weight_by_SFR_50_let_26 = JWST_mag_weight_by_SFR_50[mag_let_26_index]
    JWST_mag_weight_by_SFR_50_gtr_26 = JWST_mag_weight_by_SFR_50[mag_gtr_26_index]
    
    ax[ax_index].hist(mag_weight_by_SFR_50, bins=15, density=True, color=color_set[0], linewidth=5, histtype = 'step')
    ax[ax_index].hist(JWST_mag_weight_by_SFR_50, bins=15, density=True, color=color_set[1], linewidth=5, histtype = 'step')
    
#    ax[ax_index].hist([mag_weight_by_SFR_50_let_26, mag_weight_by_SFR_50_gtr_26], density=True, bins=15, color=[color_set[0], 'rosybrown'], linewidth=5, histtype = 'step')
#    ax[ax_index].hist([JWST_mag_weight_by_SFR_50_let_26, JWST_mag_weight_by_SFR_50_gtr_26], density=True, bins=10, color=[color_set[1], 'lightsteelblue'], linewidth=5, histtype = 'step')
    ax[ax_index].hist([10000], bins=3, label='CSST r band', color=color_set[0], histtype = 'step')
    
    ax[ax_index].hist([10000], bins=3, label='JWST F$200$W band', color=color_set[1], histtype = 'step')
#    ax[ax_index].plot([26, 26], [0,3], '--', color='black', linewidth=2)
    ax[ax_index].fill_between([26, 37],[10, 10], color='lightgrey', alpha=0.5)
    ax[ax_index].set_xlim(21, 37)
    ax[ax_index].grid()
    ax[ax_index].set_title('ID - ' + str(title_index[ax_index]), fontsize=20)
    if ax_index == 0:
        ax[ax_index].text(32,0.5,'$z_s = $' + str(zs_plot[ax_index]),fontsize=20, fontweight='bold')
        ax[ax_index].set_ylim(0,0.7)
    elif ax_index == 1:
        ax[ax_index].text(32,0.6,'$z_s = $' + str(zs_plot[ax_index]),fontsize=20, fontweight='bold')
        ax[ax_index].set_ylim(0,0.8)
        
    else:
        ax[ax_index].text(22,0.4,'$z_s = $' + str(zs_plot[ax_index]),fontsize=20, fontweight='bold')
        ax[ax_index].set_ylim(0,0.5)
#    if ax_index != 2:
#        plt.setp(ax[ax_index].get_xticklabels(), visible=False)
#    if ax_index != 0:
#        yticks = ax[ax_index].yaxis.get_major_ticks()
#        yticks[-1].label1.set_visible(False)
    if ax_index == 2:
        ax[ax_index].set_xlabel('$M_\mathrm{AB}$')
    ax[ax_index].set_ylabel('PDF')
    ax_index += 1
ax[0].legend(bbox_to_anchor=(0.83, 1.3), ncol=2)
#plt.subplots_adjust(hspace=.0)
plt.savefig('../../results/Host_Galaxy_mag.png', dpi=450)
plt.close()

'''图5，TD贝叶斯因子'''
CSST_factor_unhost_sum = np.loadtxt('../../data/IntermediaPlot_TD_Bayes_factor/CSSTBayes_factor_unhost_long_time.csv', delimiter=',')
CSST_factor_host_sum = np.loadtxt('../../data/IntermediaPlot_TD_Bayes_factor/CSSTBayes_factor_host_long_time.csv', delimiter=',')

JWST_factor_unhost_sum = np.loadtxt('../../data/IntermediaPlot_TD_Bayes_factor/JWSTBayes_factor_unhost_long_time.csv', delimiter=',')
JWST_factor_host_sum = np.loadtxt('../../data/IntermediaPlot_TD_Bayes_factor/JWSTBayes_factor_host_long_time.csv', delimiter=',')
mag_Res = np.loadtxt('../../data/SampleResult_GWGalaxy/mag_Res_120.csv', delimiter=',')
CSST_factor_host_sum_weighted = []
JWST_factor_host_sum_weighted = []
for tmp_i in range(len(mag_Res)):
    if mag_Res[tmp_i][0] < 26:
        weighted_factor = int(mag_Res[tmp_i][5] / np.min(mag_Res[:,5]))
        for tmp_j in range(weighted_factor):
            CSST_factor_host_sum_weighted.append(CSST_factor_host_sum[tmp_i])
            JWST_factor_host_sum_weighted.append(JWST_factor_host_sum[tmp_i])
            
        
    
#贝叶斯因子画图
#plt.hist(CSST_factor_unhost_sum[0:10], bins=100, density=True, histtype="stepfilled", cumulative=-1, label='CSST unhost', color=color_set[0])
#
#plt.hist(CSST_factor_host_sum, bins=40, density=True, histtype="step", cumulative=-1, label='CSST host', color=color_set[0])
#plt.hist(JWST_factor_host_sum, bins=40, density=True, histtype="step", cumulative=-1, label='JWST host', color=color_set[1])

#plt.hist(CSST_factor_host_sum_weighted, bins=100, density=True, histtype="step", cumulative=-1, label='CSST host', color=color_set[0])
plt.hist(JWST_factor_host_sum_weighted, bins=1000, density=True, histtype="step", cumulative=-1, label='JWST host', linewidth=2, color=color_set[0])
plt.hist(JWST_factor_unhost_sum, bins=1000, density=True, histtype="step", cumulative=-1, label='JWST unhost', linewidth=2, color=color_set[1])

#plt.plot([np.max(CSST_factor_unhost_sum[0:10]), np.max(CSST_factor_unhost_sum[0:10])], [0,2], '--', color=color_set[0])
plt.plot([np.max(JWST_factor_unhost_sum), np.max(JWST_factor_unhost_sum)], [0,2], '--', linewidth=2, color=color_set[1])
#plt.plot([0, 10], [0.8693650793650793, 0.8693650793650793], 'k--')
plt.xlabel('Discriminator (a.u.)')
plt.ylabel('Reversed CDF')
plt.legend(loc='best')
plt.grid()
plt.semilogx()
plt.xlim(0.00001, 2)
plt.ylim(0,1.25)
#plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 0.86, 1.0, 1.2])
plt.savefig('../../results/TD_Bayes_factor.png', dpi=450)
plt.close()

JWST_Identification = np.sum(mag_Res[:,5]) * 0.9 / 458.19174
CSST_Identification = np.sum(mag_Res[:,5]) * 0.2 / 458.19174



'''图6，挑选出事例的Bias分布图'''
#修正以后的图，替换Q_rec

BBHParam = ['mass_ratio', 'chirp_mass', 'ra', 'dec', 'theta_jn', 'psi', 'a_1', 'a_2']
#PSI: Polarization angle of the source
#theta_jn: Zenith angle between the total angular momentum and the line of sight
#a_i: Dimensionless spin magnitude of the ith object
Qrec_res = np.loadtxt('../../data/Final_Bias_show/Bias_Paper5_Selected.csv', delimiter=',')



#CDF
fig, ax = plt.subplots(1, 4, figsize=(12, 3.5))


ax0 = ax[0].hist(Qrec_res[:,0], bins=100000, density=True, cumulative=True, histtype='step', color = 'grey')

ax[0].plot([1,1],[0,1], 'k--', linewidth = 1)
ax[0].plot([2,2],[0,1], 'k--', linewidth = 2)
ax[0].plot([3,3],[0,1], 'k--', linewidth = 3)
ax[0].set_xlim(0, 3.2)
ax[0].set_ylim(0, 1)
ax[0].text(0.2,0.8,'$q$',fontsize=20)
ax[0].set_ylabel('CDF',fontsize=20)
ax[0].set_xlabel('bias$_\mathrm{ml}\ [\sigma]$',fontsize=20)
ax[0].grid()
ax[0].legend(ncol = 3, bbox_to_anchor=(1.5, 1))

ax1 = ax[1].hist(Qrec_res[:,1], bins=100000, density=True, cumulative=True, histtype='step', color = 'grey')

ax[1].plot([1,1],[0,1], 'k--', linewidth = 1)
ax[1].plot([2,2],[0,1], 'k--', linewidth = 2)
ax[1].plot([3,3],[0,1], 'k--', linewidth = 3)
ax[1].set_xlim(0, 3.2)
ax[1].set_ylim(0, 1)
ax[1].text(0.2,0.8,'$\mathcal{M}$',fontsize=20)
ax[1].grid()
ax[1].set_xlabel('bias$_\mathrm{ml}\ [\sigma]$',fontsize=20)

ax2 = ax[2].hist(Qrec_res[:,2], bins=100000, density=True, cumulative=True, histtype='step', color = 'grey')

ax[2].plot([1,1],[0,1], 'k--', linewidth = 1)
ax[2].plot([2,2],[0,1], 'k--', linewidth = 2)
ax[2].plot([3,3],[0,1], 'k--', linewidth = 3)
ax[2].set_xlim(0, 3.2)
ax[2].set_ylim(0, 1)
ax[2].text(0.2,0.8,'ra',fontsize=20)
ax[2].grid()
ax[2].set_xlabel('bias$_\mathrm{ml}\ [\sigma]$',fontsize=20)

ax3 = ax[3].hist(Qrec_res[:,3], bins=100000, density=True, cumulative=True, histtype='step', color = 'grey')

ax[3].plot([1,1],[0,1], 'k--', linewidth = 1)
ax[3].plot([2,2],[0,1], 'k--', linewidth = 2)
ax[3].plot([3,3],[0,1], 'k--', linewidth = 3)
ax[3].set_xlim(0, 3.2)
ax[3].set_ylim(0, 1)
ax[3].text(0.2,0.8,'dec',fontsize=20)
ax[3].grid()
#ax[1,3].set_xlabel('$\\frac{\\bar{x}_\\mathrm{micro} - \\bar{x}_\\mathrm{macro}}{\\sigma(x_\\mathrm{macro})}$', fontsize=20, x=-1.45)
#ax[3].set_xlabel('$\\mathrm{BIAS}$', fontsize=20, x=-1.3)
ax[3].set_xlabel('bias$_\mathrm{ml}\ [\sigma]$',fontsize=20)


plt.savefig('../../results/Rel_bias_CDF.png', dpi=450, bbox_inches = 'tight')
plt.close()









'''NA revision2 图5补充版，如果说BBH不放在中心，计算TD贝叶斯因子'''

det = 'JWST'
#引力波的时间延迟
index_GW = 120
#根据星系内的SFR和星系之间的恒星形成率，将100个随机实现以及40个宿主星系类型合并。
total_SFR_all = np.loadtxt('../../data/NA_revision2_data/Random_pin_down_GW/SFR_in_Galaxy.csv', delimiter=',')
mag_Res = np.loadtxt('../../data/NA_revision2_data/SampleResult_GWGalaxy/mag_Res_' + str(index_GW) + '.csv', delimiter=',')
factor_unhost_sum_total = []
factor_host_sum_total = []
weighted_factor_sum_all = 0
weighted_factor_in_Galaxy_all = 0
for index_GW_random in range(len(total_SFR_all)):
    factor_unhost_sum = np.sqrt(np.loadtxt('../../data/NA_revision3_data/IntermediaPlot_100weightedBySFR/' + det + '_' + str(index_GW_random) + '_' + 'Bayes_factor_unhost_long_time.csv', delimiter=','))
    factor_host_sum = np.sqrt(np.loadtxt('../../data/NA_revision3_data/IntermediaPlot_100weightedBySFR/' + det +  '_' + str(index_GW_random) + '_' + 'Bayes_factor_host_long_time.csv', delimiter=','))
    for tmp_i in range(len(mag_Res)):
        if mag_Res[tmp_i][0] < 26:
            weighted_factor_sum = int(mag_Res[tmp_i][5] / np.min(mag_Res[:,5])) * 100000
            weighted_factor_sum_all += weighted_factor_sum
            weighted_factor_in_Galaxy = int(weighted_factor_sum * total_SFR_all[index_GW_random, tmp_i] / np.sum(total_SFR_all[:,tmp_i]))
            weighted_factor_in_Galaxy_all += weighted_factor_in_Galaxy
            if weighted_factor_in_Galaxy < 1:
                print(index_GW_random, tmp_i)
            for tmp_tmp in range(weighted_factor_in_Galaxy):
                factor_host_sum_total.append(factor_host_sum[tmp_i])
    for tmptmp_i in range(len(factor_unhost_sum)):
        if factor_unhost_sum[tmptmp_i] == 0:
            factor_unhost_sum[tmptmp_i] = 10**(-15)
    factor_unhost_sum_total.append(factor_unhost_sum)
      
factor_unhost_sum_total = np.array(factor_unhost_sum_total)
    


mean_unhost = []
error_bar_unhost1 = []
error_bar_unhost2 = []
mean_host = np.mean(np.log10(factor_host_sum_total)) #
percentile_tmp = np.percentile(np.log10(factor_host_sum_total), [0.15, 99.85])
error_bar_host2 = percentile_tmp[1] - mean_host
error_bar_host1 = mean_host - percentile_tmp[0]
for plot_i in tqdm(range(len(factor_unhost_sum_total[0]))):

    mean_unhost.append(np.mean(np.log10(factor_unhost_sum_total[:,plot_i])))#
    percentile_tmp = np.percentile(np.log10(factor_unhost_sum_total[:,plot_i]), [0.15, 99.85])
    print(mean_unhost[-1], percentile_tmp, '\n')
    error_bar_unhost2.append(percentile_tmp[1] - mean_unhost[-1])
    error_bar_unhost1.append(mean_unhost[-1] - percentile_tmp[0])





#plt.fill_betweenx([0,len(factor_unhost_sum_total[0])], mean_host - error_bar_host1, mean_host + error_bar_host2, color = 'lightgrey', alpha = 1)
#
#for plot_i in range(len(factor_unhost_sum_total[0])):
#    error_range_unhost_1 = [[error_bar_unhost1[plot_i]], [error_bar_unhost2[plot_i]]]
#    
#    plt.errorbar(mean_unhost[plot_i], plot_i, xerr=error_range_unhost_1, fmt='o:', capsize=3, color=color_set[3], alpha=0.5, label='unhost')
#    
#plt.grid()
#plt.xlabel('$\log$ (Time delay Area)')
#
#plt.savefig('../../results/inconsistence.png', dpi=450)

inconsistence_unhost = []

mean_max = mean_host #np.percentile(np.log10(factor_host_sum_total), 50) #
error_bar_max = error_bar_host1
for plot_i in range(len(factor_unhost_sum_total[0])):
    inconsistence_unhost.append(np.abs(mean_unhost[plot_i] - mean_max) / np.sqrt(error_bar_host1**2 + error_bar_unhost2[plot_i]**2))
inconsistence_host = np.abs(mean_host - mean_max) / np.sqrt(error_bar_max**2 + error_bar_host1**2)

color_set = ['indianred','royalblue','goldenrod', 'black', 'seagreen']

fig, ax = plt.subplots(1,2,figsize = (10,3))



ax[0].fill_betweenx([0,len(factor_unhost_sum_total[0])], mean_host - error_bar_host1, mean_host + error_bar_host2, color = 'lightgrey', alpha = 1)

ax[0].grid()

ax[0].set_xlabel('$\log$ (Time delay Area [kpc$^2$])')

for plot_i in range(len(factor_unhost_sum_total[0])):
    error_range_unhost_1 = [[error_bar_unhost1[plot_i]], [error_bar_unhost2[plot_i]]]
    
    ax[0].errorbar(mean_unhost[plot_i], plot_i, xerr=error_range_unhost_1, fmt='o:', capsize=3, color=color_set[3], alpha=0.5, label='unhost')
    
    # plt.semilogy()


ax[1].hist(inconsistence_unhost, bins=20, color='black', alpha=0.5, label='False host')
ax[1].plot([inconsistence_host, inconsistence_host], [0, 50], color=color_set[0], linewidth=2, label='True host')
ax[1].set_xlabel('Confidence [$\sigma$]')
ax[1].set_ylabel('#')
ax[1].grid()
ax[1].set_ylim(0, 50)
ax[1].legend()
ax[0].axes.yaxis.set_ticklabels([])
plt.savefig('../../results/inconsistence.png', dpi=450)
plt.close()











'''NA revision2 图6，四像区'''
y1_caustic = np.loadtxt('../../data/NA_revision2_data/test/y1.csv', delimiter=',') 
y2_caustic = np.loadtxt('../../data/NA_revision2_data/test/y2.csv', delimiter=',') 
x_half_light = np.loadtxt('../../data/NA_revision2_data/test/x.csv', delimiter=',') / 0.9593002542900622
y_half_light = np.loadtxt('../../data/NA_revision2_data/test/y.csv', delimiter=',') / 0.9593002542900622
source_center_x = 0.005077306327901487
source_center_y = 0.31586878498210424
dis_x_set = []
dis_y_set = []

for j in tqdm(range(100000)):
    # choose a source position
    source_phi=3.2894984342950724
    source_q=0.358966
    Re_circ = 0.246138397138196 / 0.9593002542900622#circularized half-light radius
    Re_maj = Re_circ / np.sqrt(source_q) #major half-light radius
    Re_min = source_q * Re_maj
    dis_x_max = np.max([np.abs(Re_maj * np.cos(source_phi)), np.abs(Re_min * np.sin(source_phi))])+ 1
    dis_y_max = np.max([np.abs(Re_maj * np.sin(source_phi)), np.abs(Re_min * np.cos(source_phi))])+ 1 
    dis_x = random.uniform(-dis_x_max, dis_x_max)
    dis_y = random.uniform(-dis_y_max, dis_y_max)
    
    while (dis_x * np.cos(source_phi) + dis_y * np.sin(source_phi))**2 / Re_maj ** 2 + (-dis_x * np.sin(source_phi) + dis_y * np.cos(source_phi))**2 / Re_min ** 2 > 1:
        
        dis_x_max = np.max([np.abs(Re_maj * np.cos(source_phi)), np.abs(Re_min * np.sin(source_phi))])+ 1
        dis_y_max = np.max([np.abs(Re_maj * np.sin(source_phi)), np.abs(Re_min * np.cos(source_phi))])+ 1 
        dis_x = random.uniform(-dis_x_max, dis_x_max)
        dis_y = random.uniform(-dis_y_max, dis_y_max)
        
    dis_x_set.append(dis_x)
    dis_y_set.append(dis_y)




dis_x_set = np.array(dis_x_set)
dis_y_set = np.array(dis_y_set)
mod = Sersic2D(amplitude=1, r_eff=Re_circ, n=4, x_0=0, y_0=0, ellip=source_q, theta=source_phi)
light_within_half_light = mod(dis_x_set, dis_y_set) 

    
#找出阴影部分
y1_caustic_yinying = []
y2_caustic_yinying = []
x_half_light_yinying = []
y_half_light_yinying = []
for i in range(len(x_half_light)):
    index_in_y1_caustic = np.where((y1_caustic >= x_half_light[i])&((y1_caustic <= x_half_light[i] + 0.01))&(y2_caustic >= 0))[0]
    if len(index_in_y1_caustic) == 0:
        pass
    else:
        if y2_caustic[index_in_y1_caustic[0]] >= y_half_light[i]:
            y1_caustic_yinying.append(y1_caustic[index_in_y1_caustic[-1]])
            y2_caustic_yinying.append(y2_caustic[index_in_y1_caustic[-1]])
            x_half_light_yinying.append(y1_caustic[index_in_y1_caustic[-1]])
            y_half_light_yinying.append(y_half_light[i])


plt.plot(y1_caustic, y2_caustic, color = color_set[0])
plt.plot(x_half_light, y_half_light, color = color_set[1])

plt.scatter(dis_x_set + source_center_x, dis_y_set + source_center_y, c=np.log10(light_within_half_light), s=0.1)

plt.fill_between(x = y1_caustic_yinying, y1 = y2_caustic_yinying, y2 = y_half_light_yinying, color='none', alpha = 1, edgecolor="k", hatch='/////')

plt.grid()
plt.xlabel('$y_1$')
plt.ylabel('$y_2$')
plt.savefig('../../results/quadruple_region.png', dpi=450)
plt.close()


'''NA revision2 图7，能被认证出的四像红移和天区定位'''

z_s_set = np.loadtxt('../../data/NA_revision2_data/PairResult_quadruple/z_s_set.csv', delimiter=',')

max_sky_localization = np.loadtxt('../../data/NA_revision2_data/PairResult_quadruple/max_sky_localization.csv', delimiter=',')

plt.scatter(z_s_set, max_sky_localization)
plt.semilogy([2.1, 2.1], [0, 15], 'k--')
plt.plot([0, 3], [5, 5], 'k--')
plt.xlim(1, 3)
plt.ylim(0, 15)
plt.grid()
plt.xlabel('$z_s$')
plt.ylabel('localization area')
plt.savefig('../../results/z_s_vs_sky_localization_area.png', dpi=450)
plt.close()


'''NA revision2 图8， 能被认证出的四像微透镜'''
'''图1'''
Mz_cut_set = [20]
fig, ax = plt.subplots(1, 1, figsize=(6.5, 2.5))
#gs1 = plt.GridSpec(2, 1)
#gs.update(wspace=0.4)
#plt.figure(figsize=(6.5,6.5))
for i in range(len(Mz_cut_set)):
    Mz_cut = Mz_cut_set[i]
    low_SNR_envelope = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    low_envelope = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_SNR_envelope = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_envelope = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
        
    SNR_tmp_lensed_sum = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_mean_lensed_sum = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/match_mean_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    theory_match_sum_tmp = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/theory_match' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_std_lensed_sum = np.loadtxt('../../data/NA_revision2_data/Result_file4plot_quadruple/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')

    ax.fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='lightgrey', label='Chirp Mass $\geq$' + str(int(Mz_cut)) + '$\mathrm{M}_\odot$')
   
    ax.scatter(SNR_tmp_lensed_sum, match_mean_lensed_sum, c='#0C5DA5', s = 10, marker='p')
    ax.errorbar(SNR_tmp_lensed_sum, match_mean_lensed_sum, yerr = 1.5*np.array(match_std_lensed_sum), elinewidth=1, ecolor='black', capsize=3, capthick=1, linestyle='none')
    
    ax.set_xlabel('SNR')
    ax.set_ylabel('Match')
    ax.grid()
    ax.set_ylim(0.935,1)
    ax.set_xlim(35,500)

    ax.legend()

plt.savefig('../../results/match_quadruple.png', dpi=450)
plt.close()


#'''NA revison2 图9，power law模型的结果'''
#Result_power_law = np.loadtxt('../../data/NA_revision2_data/120_0corner_Res_data_power_law.csv', delimiter=',')
#Result_SIE = np.loadtxt('../../data/NA_revision2_data/120_0corner_Res_data_SIE.csv', delimiter=',')
#labels_new = [r"$t_0$", r"$t_1$", r"$t_2$", r"$t_3$"]
#corner.corner(np.array(Result1), labels=labels_new, show_titles=True, truths=[kwargs_lens[0]['theta_E'], kwargs_lens[0]['gamma'], kwargs_lens[0]['center_x'], kwargs_lens[0]['center_y'], kwargs_lens[0]['e1'], kwargs_lens[0]['e2'], t_days_truths[0], t_days_truths[1], t_days_truths[2], t_days_truths[3], 0])


'''NA revision2 图9， 微透镜几何光学对比图'''
M_L = 0.37804676999999448084
z_L = 0.4016984227892334
M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
kappa = 0.3810082110922153;
gamma = 0.38100891939996284
mu = 1/((1-kappa)**2 - gamma**2)**0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
constant = 2*np.pi*mu/coeffi
cw = coeffi/2/np.pi

#对比一下Paper4中给referee的回复中所得到的结果
with open('../../data/ResultMinimum/root_2.txt', 'r') as f0:
    data = f0.readlines()
    
fermat_pot = []
mag_geo = []
for i in range(len(data)):
    fermat_pot.append(float(data[i].split(' ')[12]))
    mag_geo.append(float(data[i].split(' ')[8]))
    
yi = 0

f1=open('../../data/ResultMinimum/TimeLength_min_' + str(yi) + '.bin',"rb")
TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
f1.close()
TimeLengthFile_min = np.array(TimeLengthFile_min)
fileArea = "../../data/ResultMinimum/adptive_Area_min_" + str(yi) + ".bin"
fileTime = "../../data/ResultMinimum/adptive_Time_min_" + str(yi) + ".bin"

f1=open(fileArea,"rb")
Area=struct.unpack("d"*TimeLengthFile_min[0], f1.read(8*TimeLengthFile_min[0]))
f1.close()
Area = np.array(Area)

f2=open(fileTime,"rb")
time=struct.unpack("d"*TimeLengthFile_min[0], f2.read(8*TimeLengthFile_min[0]))
f2.close()
time = np.array(time)


index_none_zero = np.where(Area>0)[0][0]
time = time[index_none_zero::]
Area = Area[index_none_zero::]



Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
time_raw = time[0:-1]
time_zero = time_raw[0]

time_raw = time_raw - time_zero
time_geo = coeffi * np.array(fermat_pot) - time_zero
#下面两个变量是为了后面与Paper4作比较的
time_raw_new = copy.deepcopy(time_raw)
Ft_raw_new = copy.deepcopy(Ft_raw)


delta_time_raw = time_raw[1] - time_raw[0] 



"""去掉尾巴"""
length = len(time_raw)
time_raw = time_raw[0:length*4//5]
Ft_raw = Ft_raw[0:length*4//5]
inter_time_raw = interp1d(time_raw, Ft_raw)
time_raw_plot = time_raw
Ft_raw_plot = Ft_raw


"""without apodization"""
Ft_raw[1::] = Ft_raw[1::] - constant
Ft_raw[0] = Ft_raw[0] - constant/2

time_new = time_raw
Ft_new = Ft_raw

"""自己写傅里叶变换"""
omegafreq = np.arange(0.1 * 2 * np.pi,2000*2*np.pi,2*np.pi)
freq = omegafreq/2/np.pi
Freal = []
Fimag = []
for i in tqdm(range(len(omegafreq))):
    Freal.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
    Fimag.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
    
Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
Ff = Ff* omegafreq / complex(0,1)*cw
Ffsgn = constant * cw 
Ff = Ff + Ffsgn

   
Ffabs = np.abs(Ff)
ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
# Ffabs = signal.savgol_filter(Ffabs,3,1)

"""理论值"""
TypeINum = 0
TypeIINum = 0
TypeIIINum = 0
FfTheory = 0
for i in range(len(fermat_pot)):
    
    if mag_geo[i] >0:
        
        FfTheory += mag_geo[i]**0.5*np.exp((2*np.pi*freq*time_geo[i])*complex(0,1))
        TypeINum += 1
         
    else:
        FfTheory += np.abs(mag_geo[i])**0.5*np.exp((2*np.pi*freq*time_geo[i] - np.sign(freq)* np.pi/2)*complex(0,1))
        TypeIINum += 1
        
print("TypeI number = ", TypeINum)
print("TypeII number = ", TypeIINum)
print("TypeIII number = ", TypeIIINum)
FfTheoryAbs = np.abs(FfTheory)
FfTheoryTheta = np.real(-complex(0,1)*np.log(FfTheory/FfTheoryAbs))
    
"""以上"""

plt.figure(figsize=(10,3))

grid = plt.GridSpec(4, 8, hspace = 0.5, wspace=1.5)
main_ax = plt.subplot(grid[0:4,0:4])
main_ax.semilogx(time_raw_plot,Ft_raw_plot + constant,color=color_set[1])
main_ax.scatter(time_geo[1::], inter_time_raw(time_geo[1::]), label='Geometric image', color=color_set[0], marker='*')
main_ax.semilogx([time_raw[0],time_raw[-1]+0.1],[constant,constant], '--',color='grey')
main_ax.set_xlabel('t[s]')
main_ax.set_ylabel('F(t)')
main_ax.set_xlim(10**(-6),5*10**(-1))
main_ax.set_ylim(0.5*10**(6),1.5*10**(6))
main_ax.legend()
main_ax.grid()

sub_ax1 = plt.subplot(grid[0:2,4:8])

sub_ax1.semilogx(freq, Ffabs, label='Full diffraction integral',color=color_set[1]) #//FIXME
sub_ax1.semilogx(freq, FfTheoryAbs, label='Geometrical optics limit',color=color_set[0])
# plt.plot(freq, FfabsOld)
# plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
# plt.plot(freq, Ffabs1, label='Cut')
   
sub_ax1.grid()
sub_ax1.legend()
# plt.ylim(1,3)
#sub_ax1.set_xlabel('f[Hz]')
sub_ax1.set_ylabel('$|F(f)|$')
sub_ax1.set_xlim(0.1, 3000)
# plt.savefig("Ffembed.png",dpi=450)
"""相位"""
sub_ax2 = plt.subplot(grid[2:4,4:8])

sub_ax2.semilogx(freq, ThetaF,color=color_set[1]) #//FIXME
sub_ax2.semilogx(freq, FfTheoryTheta,color=color_set[0])
# plt.plot(freq, ThetaF, label='Cut')
# plt.loglog(test1,test2)
# plt.ylim(-0.2,0.2)
sub_ax2.grid()
sub_ax2.legend()
sub_ax2.set_xlabel('f[Hz]')
sub_ax2.set_ylabel(r'$\theta_F$')
sub_ax2.set_xlim(0.1, 3000)
plt.savefig('../../results/micro_image.png', dpi=450)
plt.close()






'''NA revision3 + precessing， 能被认证出的微透镜'''
'''图1'''
Mz_cut_set = [20]
fig, ax = plt.subplots(1, 1, figsize=(6.5, 2.5))
#gs1 = plt.GridSpec(2, 1)
#gs.update(wspace=0.4)
#plt.figure(figsize=(6.5,6.5))
for i in range(len(Mz_cut_set)):
    Mz_cut = Mz_cut_set[i]
    low_SNR_envelope_precessing = np.loadtxt('../../data/Result_file4plot_precessing/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    low_envelope_precessing = np.loadtxt('../../data/Result_file4plot_precessing/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    
    
    
    
    
    low_SNR_envelope = np.loadtxt('../../data/Result_file4plot/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    low_envelope = np.loadtxt('../../data/Result_file4plot/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_SNR_envelope = np.loadtxt('../../data/Result_file4plot/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_envelope = np.loadtxt('../../data/Result_file4plot/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
        
    
    SNR_tmp_lensed_sum = np.loadtxt('../../data/Result_file4plot/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_mean_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_mean_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    theory_match_sum_tmp = np.loadtxt('../../data/Result_file4plot/theory_match' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_std_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    index_match_let_0p95 = np.where(match_mean_lensed_sum < 0.95)
#    match_mean_lensed_sum[index_match_let_0p95] = 0.95 + index_match_let_0p95[0] / 1000
    match_lensed_up = match_mean_lensed_sum + 1.5 * match_std_lensed_sum
#    match_std_lensed_sum[np.where(match_lensed_up>1)] /= 5
    ax.fill_between(low_SNR_envelope_precessing[0:6], low_envelope_precessing[0:6], 1, interpolate=True, color='lightgrey', label='W/ precession')
    ax.fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='lightgrey')
    ax.fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='grey', label='W/o precession')
    
    ax.scatter(SNR_tmp_lensed_sum, match_mean_lensed_sum, c='#0C5DA5', s = 10, marker='p')
    ax.errorbar(SNR_tmp_lensed_sum, match_mean_lensed_sum, yerr = 1.5*np.array(match_std_lensed_sum), elinewidth=1, ecolor='black', capsize=3, capthick=1, linestyle='none')
    
    ax.set_xlabel('SNR')
    ax.set_ylabel('Match')
    ax.grid()
    ax.set_ylim(0.94,1)
    ax.set_xlim(35,500)

    ax.legend()

plt.savefig('../../results/match_precessing.png', dpi=450)
plt.close()










'''图3 sigma'''
Mz_cut = 20
low_SNR_envelope = np.loadtxt('../../data/Result_file4plot/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
low_envelope = np.loadtxt('../../data/Result_file4plot/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
up_SNR_envelope = np.loadtxt('../../data/Result_file4plot/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
up_envelope = np.loadtxt('../../data/Result_file4plot/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    
low_SNR_envelope_precessing = np.loadtxt('../../data/Result_file4plot_precessing/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
low_envelope_precessing = np.loadtxt('../../data/Result_file4plot_precessing/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')


SNR_tmp_lensed_sum = np.loadtxt('../../data/Result_file4plot/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
match_mean_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_mean_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
theory_match_sum_tmp = np.loadtxt('../../data/Result_file4plot/theory_match' + str(int(Mz_cut)) + '.csv', delimiter=',')
match_std_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')

low_envelope_func = interp1d(low_SNR_envelope, low_envelope, bounds_error=False, fill_value=(low_envelope[0], low_envelope[-1]))
low_envelope_func_precessing = interp1d(low_SNR_envelope_precessing, low_envelope_precessing, bounds_error=False, fill_value=(low_envelope_precessing[0], low_envelope_precessing[-1]))


sigma_lensed = []
for i in range(len(SNR_tmp_lensed_sum)):
    low_env_tmp = low_envelope_func(SNR_tmp_lensed_sum[i])
    low_env_tmp_precessing = low_envelope_func_precessing(SNR_tmp_lensed_sum[i])
    
    sigma_lensed.append((low_env_tmp - match_mean_lensed_sum[i])/match_std_lensed_sum[i])

#plt.fill_between([40, 700], [0, 0], [5, 5], color='lightgrey')
plt.scatter(SNR_tmp_lensed_sum, sigma_lensed, marker='o', edgecolors='grey', c='none')
plt.plot([40, 700], [5, 5], 'k--')
plt.ylim(0.5, 130)
plt.xlim(50,700)
plt.semilogy()
plt.grid()
plt.xlabel('SNR')
plt.ylabel('Metric distance [$\sigma$]')
plt.yticks([0.5, 1, 5, 10, 100], [0.5, 1, 5, 10, 100])
plt.text(450, 10, '$85\%$', fontsize=20, color=color_set[0])
plt.savefig('../../results/evidence.png', dpi=450)
plt.close()
np.sum(np.array(sigma_lensed)>5)/np.sum(np.array(sigma_lensed)>0)





'''noisy data match bilby'''
Mz_cut_set = [20]
fig, ax = plt.subplots(1, 1, figsize=(6.5, 2.5))
#gs1 = plt.GridSpec(2, 1)
#gs.update(wspace=0.4)
#plt.figure(figsize=(6.5,6.5))
for i in range(len(Mz_cut_set)):
    Mz_cut = Mz_cut_set[i]
    
    
    
    
    low_SNR_envelope = np.loadtxt('../../data/NA_revision3_data/match_noisy_inject_VS_Bilby/low_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    low_envelope = np.loadtxt('../../data/NA_revision3_data/match_noisy_inject_VS_Bilby/low_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_SNR_envelope = np.loadtxt('../../data/NA_revision3_data/match_noisy_inject_VS_Bilby/up_SNR_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
    up_envelope = np.loadtxt('../../data/NA_revision3_data/match_noisy_inject_VS_Bilby/up_envelope' + str(int(Mz_cut)) + '.csv', delimiter=',')
        
    
    SNR_tmp_lensed_sum = np.loadtxt('../../data/Result_file4plot/SNR_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
    match_mean_lensed_sum = np.loadtxt('../../data/NA_revision3_data/match_noisy_inject_VS_Bilby/micro_match_noise.csv', delimiter=',')
    match_std_lensed_sum = np.loadtxt('../../data/Result_file4plot/match_std_lensed' + str(int(Mz_cut)) + '.csv', delimiter=',')
#    index_match_let_0p95 = np.where(match_mean_lensed_sum < 0.95)
#    match_mean_lensed_sum[index_match_let_0p95] = 0.95 + index_match_let_0p95[0] / 1000
    match_lensed_up = match_mean_lensed_sum + 1.5 * match_std_lensed_sum
#    match_std_lensed_sum[np.where(match_lensed_up>1)] /= 5
    
    ax.fill_between(low_SNR_envelope, low_envelope, 1, interpolate=True, color='grey', label='Chirp Mass $\geq$' + str(int(Mz_cut)) + '$\mathrm{M}_\odot$')
    
    ax.scatter(SNR_tmp_lensed_sum, match_mean_lensed_sum, c='#0C5DA5', s = 10, marker='p')
    ax.errorbar(SNR_tmp_lensed_sum, match_mean_lensed_sum, yerr = 1.5*np.array(match_std_lensed_sum), elinewidth=1, ecolor='black', capsize=3, capthick=1, linestyle='none')
    
    ax.set_xlabel('SNR')
    ax.set_ylabel('Match')
    ax.grid()
    ax.set_ylim(0.65,1)
    ax.set_xlim(45,500)

    ax.legend()

plt.savefig('../../results/match_noisy.png', dpi=450)
plt.close()





'''chi_p'''
with open('../../data/outdir_Lensed_Saddle_with_microprecessing_IMRPhenomPv2/outdir18073_0/H1_L1_V1_result.json') as f1:
    data_lensed = json.load(f1)

chi_p = data_lensed['posterior']['content']['chi_p']

plt.hist(chi_p, density=True, histtype='step')
plt.xlabel('$\chi_p$')
plt.ylabel('PDF')
plt.grid()
plt.savefig('../../results/chi_p.png', dpi=450)
plt.close()





#waveform to referee
h1f_rec_f = np.loadtxt('../../data/To_referee_figure/h1f_rec_f.csv', delimiter=',')
h1f_rec = np.loadtxt('../../data/To_referee_figure/h1f_rec.csv', delimiter=',')
inject_freq = np.loadtxt('../../data/To_referee_figure/inject_freq.csv', delimiter=',')
inject_signal = np.loadtxt('../../data/To_referee_figure/inject_signal.csv', delimiter=',')
freq_h1 = np.loadtxt('../../data/To_referee_figure/freq_h1.csv', delimiter=',')
tilde_signal_h1 = np.loadtxt('../../data/To_referee_figure/tilde_signal_h1.csv', delimiter=',')

signal_h1_sample_times = np.loadtxt('../../data/To_referee_figure/signal_h1_sample_times.csv', delimiter=',')
signal_h1 = np.loadtxt('../../data/To_referee_figure/signal_h1.csv', delimiter=',')
inject_time = np.loadtxt('../../data/To_referee_figure/inject_time.csv', delimiter=',')
inject_signal_time= np.loadtxt('../../data/To_referee_figure/inject_signal_time.csv', delimiter=',')




fig, axs = plt.subplots(4, 1, figsize=(8, 12)) 

axs[0].semilogy(h1f_rec_f, np.abs(h1f_rec), 'grey', label='cWB', linewidth=1.5)
axs[0].semilogy(inject_freq, np.abs(inject_signal), color='#ff7f0e', label='Injected', linewidth=1.5)
# plt.semilogy(inject_freq[index_between_20_400], np.abs(inject_signal[index_between_20_400]/(Ff[index_between_20_400]*np.exp(complex(0,1) *ThetaF[index_between_20_400])))*Ff[33], label='Bilby')
axs[0].semilogy(freq_h1, np.abs(tilde_signal_h1), color=color_set[1], label='Bilby', linewidth=1.5)
axs[0].set_xlim(20, 220)
axs[0].set_ylim(5*10**(-26), 10**(-23.5))
axs[0].legend()
axs[0].grid()
axs[0].set_xlabel('f[Hz]')
axs[0].set_ylabel('strain')

axs[1].plot(signal_h1_sample_times, signal_h1, color=color_set[1], label='Bilby zoom in')
axs[1].plot(inject_time, inject_signal_time, color='#ff7f0e', label='Injected zoom in')
axs[1].set_ylim(-5 * 10**(-23), 5 * 10**(-23))
axs[1].set_xlim(inject_time[0]+31.55, inject_time[0]+31.85)
axs[1].legend()
axs[1].set_ylabel('strain')
axs[1].grid()
# axs[1].plot(inject_time, inject_signal_time/2)
# # axs[1].plot(time_plot, func_inj_time(time_plot)/func_bilby_time(time_plot))
# axs[1].set_ylim(-5 * 10**(-23), 5 * 10**(-23))
# axs[1].set_xlim(inject_time[0]+31.5, inject_time[0]+31.8)

axs[2].plot(signal_h1_sample_times, signal_h1, color=color_set[1], label='Bilby')
axs[2].plot([1017404076, 1017404076], [-5 * 10**(-23), 5 * 10**(-23)], color=color_set[0])
axs[2].set_ylim(-5 * 10**(-23), 5 * 10**(-23))
axs[2].set_xlim(inject_time[0]+18, inject_time[0]+32)
axs[2].legend()
axs[2].set_ylabel('strain')
axs[2].grid()

axs[3].plot(inject_time, inject_signal_time, color='#ff7f0e', label='Injected')
axs[3].plot([1017404076, 1017404076], [-5 * 10**(-23), 5 * 10**(-23)], color=color_set[0])
# axs[1].plot(time_plot, func_inj_time(time_plot)/func_bilby_time(time_plot))
axs[3].set_ylim(-5 * 10**(-23), 5 * 10**(-23))
axs[3].set_xlim(inject_time[0]+18, inject_time[0]+32)
axs[3].legend(loc='upper left')
axs[3].set_xlabel('t[s]')
axs[3].set_ylabel('strain')
axs[3].grid()
plt.savefig('../../results/waveform_freq_time.png', dpi=450)
plt.close()