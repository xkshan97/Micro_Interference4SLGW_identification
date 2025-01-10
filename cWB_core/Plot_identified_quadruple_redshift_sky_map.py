import numpy as np
import matplotlib.pyplot as plt
import json


z_s_set = np.loadtxt('./PairResult_quadruple/z_s_set.csv', delimiter=',')

max_sky_localization = np.loadtxt('./PairResult_quadruple/max_sky_localization.csv', delimiter=',')

plt.scatter(z_s_set, max_sky_localization)
plt.semilogy([2.1, 2.1], [0, 15], 'k--')
plt.plot([0, 3], [5, 5], 'k--')
plt.xlim(1, 3)
plt.ylim(0, 15)
plt.grid()
plt.xlabel('$z_s$')
plt.ylabel('localization area')
plt.savefig('./PairResult_quadruple/z_s_vs_sky_localization_area.png', dpi=450)
