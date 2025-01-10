import matplotlib.pyplot as plt
import numpy as np
from scipy import special as sp

# class Sersic_stellar_density():
def K_4():
    return 7.669 #1.0857 * np.exp(0.6950 + np.log(n) - 0.1789 / n)

def M_tot(theta_E, theta_eff):
    return np.pi * theta_E * theta_eff / 2

def Q_4(theta_eff):
    #lower incomplete Gamma function.
    #https://handwiki.org/wiki/Astronomy:Sersic_profile
    k_4 = K_4()
    return 2 * np.pi * 4 * theta_eff ** 2 * k_4 ** (-2 * 4) * sp.gammainc(8, 2 ** (-1/4) * k_4) * sp.gamma(8)

def A_4(theta_E, theta_eff, sigma_v):
    #sigma_v in unit of km/s
    f_DM = 0.8 * np.log10(sigma_v / 100) + 0.05 #Verdoras里面是错的
    return (1 - f_DM) * M_tot(theta_E, theta_eff) / Q_4(theta_eff)

def kappa_s4(x1, x2, q, theta_E, theta_eff, sigma_v):
    omega = np.sqrt(x1 ** 2 + q ** 2 * x2 ** 2)
    #在后面lens采样的时候，用的lenstronomy，然后假设主轴是在y轴上，这个与Kormann 1994 SIE里面的主轴是相融洽的。
    
    return A_4(theta_E, theta_eff, sigma_v) * np.exp(-K_4() * (omega / np.sqrt(q) / theta_eff) ** (1/4))

def kappa_total(x1, x2, q, theta_E):
    omega = np.sqrt(x1 ** 2 + q ** 2 * x2 ** 2)
    #在后面lens采样的时候，用的lenstronomy，然后假设主轴是在y轴上，这个与Kormann 1994 SIE里面的主轴是相融洽的。
    
    return 1/2 * theta_E * np.sqrt(q) / omega


# x1 = np.arange(0.01, 1.8, 0.01)
# x2 = np.arange(0.01, 1.8, 0.01)
# r = np.sqrt(x1 ** 2 + x2 ** 2)
# theta_E = 0.7
# theta_eff = theta_E / 0.576
# q = 0.698
# sigma_v = 204
# kappa_total_test = kappa_total(x1, x2, q, theta_E)
# kappa_s4_test = kappa_s4(x1, x2, q, theta_E, theta_eff, sigma_v)

# plt.plot(r, kappa_total_test, label = '$\kappa_\mathrm{total}$')
# plt.plot(r, kappa_s4_test, '--', label = '$\kappa_*$')

# plt.plot([0.42, 0.42], [0, 4], 'r', label = '$\\theta_\mathrm{eff}$')
# plt.plot([0.7, 0.7], [0, 4], 'grey', label = '$\\theta_\mathrm{Ein}$')
# plt.plot([0.8, 0.8], [0, 4], 'grey', label = '$\\theta_\mathrm{image}$')
# plt.ylim(0, 2)
# plt.grid()
# plt.legend()
# plt.xlabel('$\\theta[\mathrm{arcsec}]$')
# plt.ylabel('$\kappa$')
# plt.savefig('IntermediaPlot/Vernardos.png', dpi=450)
# kappa_s_image_C = kappa_s4_test[np.where(r >= 0.8)[0][0]]
# kappa_total_image_C = kappa_total_test[np.where(r >= 0.8)[0][0]]
# s = 1 - kappa_s_image_C / kappa_total_image_C
