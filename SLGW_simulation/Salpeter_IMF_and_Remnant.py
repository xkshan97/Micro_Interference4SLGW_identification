import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from MySample import sample
from scipy import stats
from scipy.optimize import curve_fit
#Salpeter Initial Mass Function 
#Sepra 2015 Remnant


def Salpeter_IMF(m):
    return m**(-2.35)




def Remnant_Mass_bar(M_zams):
    #下面这些分段函数的系数我都自己微调了，保证首尾相接
    #M_zams: zero age main sequence mass.
    #
    if M_zams < 1.5:
        print('ERROR! ZAMS STILL ALIVE! ')
        
    M_up = 7
    if M_zams < M_up:
        #不是sepra 2015提供的结果，是用Meena 2022中使用的不公开数据拟合的（J. J. Eldridge 2017图17）
        return (1.4 - 0.1) / (7 - 0.1) * (M_zams - 0.1) + 0.1
    
    elif M_zams >= M_up and M_zams <= 13:
        return 1.4
    elif M_zams >13 and M_zams <= 16:
        return 0.2026 * M_zams - 1.2338
    elif M_zams > 16 and M_zams <= 27.166664:
        return (0.041 * M_zams ** 3 - 0.673 * M_zams ** 2 + 2.18 * M_zams + 0.361) / (0.952 * M_zams + 0.15)
    elif M_zams > 27.166664 and M_zams <= 36:
        return (0.0563 * M_zams ** 3 - 1.1 * M_zams ** 2 + 2.49 * M_zams + 0.318) / (0.952 * M_zams + 0.15)
    elif M_zams > 36:
        return 0.9135 * M_zams + 4.6217
    
def Remnant_Mass(M_zams):
    return 0.9 * Remnant_Mass_bar(M_zams)

M_zmas_test = np.arange(0.1, 100, 0.001)
Remnant_test = []
for i in range(len(M_zmas_test)):
    Remnant_test.append(Remnant_Mass(M_zmas_test[i]))


plt.plot(M_zmas_test, Remnant_test)
plt.grid()
plt.xlabel('$M_\mathrm{ZEMS}$')
plt.ylabel('$M_\mathrm{remnant}$')
plt.savefig('./IntermediaPlot/Initial_final_relation.png', dpi=450)


#采样Salpeter IMF，得到100万个样本
m_min = 1.5
Salpeter_sample = sample(Salpeter_IMF, 1000000, lower_bd=m_min, upper_bd=100, guess=m_min)

#Salpeter IMF解析公式结果

m_test = np.arange(m_min, 100, 0.001)

Salpeter_test = Salpeter_IMF(m_test)
Salpeter_test = Salpeter_test / np.sum(Salpeter_test) / 0.001

plt.hist(Salpeter_sample, bins=10000, density=True, label='Sampled')
plt.loglog(m_test, Salpeter_test, label = 'Salpeter IMF Analytic')
plt.xlabel('m')
plt.ylabel('dn/dm')
plt.grid()
plt.legend()
plt.savefig('./IntermediaPlot/m_min = ' + str(m_min) + 'Salpeter_sample4remnant.png', dpi=450)

#将上面100万个样本转换到remnant，看一下remnant的质量函数
Remnant_MF = []
for i in range(len(Salpeter_sample)):
    Remnant_MF.append(Remnant_Mass(Salpeter_sample[i]))

Remnant_MF = np.array(Remnant_MF)

BH_frac = np.sum(np.array(Remnant_MF)>3) / len(Remnant_MF)
Gtr1_frac = np.sum(np.array(Remnant_MF)>1) / len(Remnant_MF)
SMBH_frac = np.sum(np.array(Remnant_MF)>25) / len(Remnant_MF)

print('BH_frac: ', BH_frac)
print('SMBH_frac: ', SMBH_frac)
print('> 1 $M_\odot$: ', Gtr1_frac)

ln_Remnant_MF = np.log(Remnant_MF)
ln_Remnant_MF_pdf = stats.relfreq(ln_Remnant_MF, numbins=100) # numbins 是统计一次的间隔(步长)是多大
ln_Remnant_MF_pdf_value = ln_Remnant_MF_pdf.frequency / ln_Remnant_MF_pdf.binsize
ln_Remnant_M = ln_Remnant_MF_pdf.lowerlimit + np.linspace(0, ln_Remnant_MF_pdf.binsize * ln_Remnant_MF_pdf.frequency.size, ln_Remnant_MF_pdf.frequency.size)

Remnant_M = np.exp(ln_Remnant_M)
Remnant_MF_pdf_value = ln_Remnant_MF_pdf_value / Remnant_M


plt.hist(Remnant_MF, bins=100, density=True)
plt.semilogy(Remnant_M, Remnant_MF_pdf_value)

plt.loglog()
plt.xlim(1, 105)
plt.xlabel('remnant m')
plt.ylabel('dn/dm')
plt.savefig('./IntermediaPlot/Salpeter2Remnant_MF')

np.savetxt('IntermediaOut/Remnant_MF.csv',[Remnant_M, Remnant_MF_pdf_value], delimiter=',')
# def Curve_fit_ln_Remnant_PDF(M, a, b, c, d):
#     return a * np.log(M) ** 3 + b * np.log(M) ** 2 + c * np.log(M) + d

# popt, pcov = curve_fit(Curve_fit_ln_Remnant_PDF, Remnant_M, np.log(Remnant_MF_pdf_value))

# def Fitted_Curve_fit_ln_Remnant_PDF(M, a = popt[0], b = popt[1], c = popt[2], d = popt[3]):
#     return  a * np.log(M) ** 3 + b * np.log(M) ** 2 + c * np.log(M) + d

# Fitted_Remnant_PDF = np.exp(Fitted_Curve_fit_ln_Remnant_PDF(Remnant_M))
# plt.semilogy(Remnant_M, Remnant_MF_pdf_value)
# plt.plot(Remnant_M, Fitted_Remnant_PDF, '--')
# plt.xlim(1, 105)
