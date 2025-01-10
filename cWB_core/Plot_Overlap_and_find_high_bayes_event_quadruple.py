import numpy as np
import matplotlib.pyplot as plt


Bayes_Unlens = np.loadtxt('./Bayes_Unlens/overlap.csv', delimiter=',')

with open('./PairResult_quadruple/identified_index.csv', 'r') as f0:
    Selected_index = f0.readlines()
    
with open('./PairResult_quadruple/identified_file_name_not_duplication.csv', 'r') as f1:
    Selected_file_name_not_duplication = f1.readlines()
    
z_s_set = np.loadtxt('./PairResult_quadruple/z_s_set.csv', delimiter=',')

#计算FAP per pair
Bayes_factor_sum = []
for i in range(len(Selected_index)):
    FAP_per_pair = []
    try:
        Bayes_factor = [float(np.loadtxt('./Bayes_Lens_quadruple/' + Selected_index[i].split('\n')[0] + '.csv', delimiter=','))]
    except TypeError:
        Bayes_factor = np.loadtxt('./Bayes_Lens_quadruple/' + Selected_index[i].split('\n')[0] + '.csv', delimiter=',')
    for j in range(len(Bayes_factor)):
        Bayes_factor_sum.append(Bayes_factor[j])
        FAP_per_pair_tmp = np.sum(Bayes_Unlens >= Bayes_factor[j]) / len(Bayes_Unlens)
        if FAP_per_pair_tmp == 0:
            FAP_per_pair.append(10**(-10 + j))
        else:
            FAP_per_pair.append(FAP_per_pair_tmp)
        np.savetxt('./Bayes_Lens_quadruple/FAP_per_pair_' + Selected_index[i].split('\n')[0] + '.csv', FAP_per_pair, delimiter=',')

#画图，x axis是FAP per pair，纵坐标是index，笔记本上的图三就是来自于这里。
FAP_per_pair_sum = []
fig = plt.figure(figsize=(5,15))
ax2 = fig.add_subplot(111)
ax1 = ax2.twiny()
for i in range(len(Selected_index)):
    try:
        FAP_per_pair = [float(np.loadtxt('./Bayes_Lens_quadruple/FAP_per_pair_' + Selected_index[i].split('\n')[0] + '.csv', delimiter=','))]
    except TypeError:
        FAP_per_pair = np.loadtxt('./Bayes_Lens_quadruple/FAP_per_pair_' + Selected_index[i].split('\n')[0] + '.csv', delimiter=',')
    for j in range(len(FAP_per_pair)):
        FAP_per_pair_sum.append(FAP_per_pair[j])
        ax1.scatter(FAP_per_pair[j], i, color='grey', alpha=0)
        if FAP_per_pair[j] * 528 < 0.1:
            ax2.scatter(FAP_per_pair[j] * 528, i, color='red', marker='*', s=100)
        else:
            ax2.scatter(FAP_per_pair[j] * 528, i, color='grey')
            
ax1.semilogx()
ax2.semilogx()
ax1.set_xlim(10**(-11), 10)
ax2.set_xlim(10**(-11) * 528, 10 * 528)
ax1.set_ylim(-0.5, len(Selected_index) - 0.5)
ax2.set_ylim(-0.5, len(Selected_index) - 0.5)
ax1.set_yticks(range(len(Selected_index)))
ax2.set_yticks(range(len(Selected_index)))
ax1.set_xlabel('FAP/pair')
ax2.set_xlabel('FAP/year')
ax2.set_xticks([10**(-1), 10**(0), 10**(1), 10**(2)], ['$<0.1$', '$1$', '$10$', '$100$'])
ax1.set_xticks([10**(-5), 10**(-3), 10**(-1), 10**(1)])
ax2.grid()




       
plt.hist(Bayes_factor_sum, bins=100, cumulative=True, density=True, histtype='step')

plt.hist(Bayes_Unlens, bins=10000, cumulative=True, density=True, histtype='step')
plt.loglog()
