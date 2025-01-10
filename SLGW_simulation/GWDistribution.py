import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
import astropy.units as u
from astropy.constants import c
from MySample import sample
import random
from scipy import interpolate as sci_inter
from scipy import special as sp
from tqdm import tqdm



#下面来自MD14，metallicity = 0.3（NA修改版本还是没加这个，因为加了以后跟徐菲的不一样，因为用了以后跟徐菲正文里结果不一样）, t_d = 50 Myr===================
def MetallicityCut(frac_meta, z):
    #frac_meta is the matallicity fraction of solar metallicity.
    #z is the redshift.
    #credit Langer & Norman 2006 file name : Metallicity Cut at paper4
    
    return sp.gammainc(0.84, frac_meta**2 * 10**(0.3 * z)) 

def MDSFR(z):
    #z is the redshift
    # in unit of M_sun/Mpc^(3)/yr^{1}
    #credit Fei Xu 2022
    return 0.015 * (1 + z) ** 2.7 / (1 + ((1 + z)/2.9)**5.6 )


def FindLookBackRedshift(z, MinDelayTime):
    #find the minimum redshift of formation.
    #z is the merger redshift
    #MinDelayTime is the minimum time between formation and merger, in unit of Myr^{-1}
    z_tmp = z
    LookBackTime = Planck18.lookback_time(z_tmp).value * (u.Gyr).to(u.Myr) - Planck18.lookback_time(z).value * (u.Gyr).to(u.Myr)
    while LookBackTime < MinDelayTime:
        z_tmp += 0.0001
        LookBackTime = Planck18.lookback_time(z_tmp).value * (u.Gyr).to(u.Myr) - Planck18.lookback_time(z).value * (u.Gyr).to(u.Myr)
        
    return z_tmp 




def MergeRateZero(MinDelayTime, frac_meta):
    #merge rate at redshift 0 in order to nomalize MergeRate function.
    
    MinFormationz = FindLookBackRedshift(0, MinDelayTime)
    Formationz_set = np.exp(np.linspace(np.log(MinFormationz), np.log(1100),1000))
    phi = MetallicityCut(frac_meta, Formationz_set)
    MDrho = MDSFR(Formationz_set) 
    phi = 1#(phi[1::] + phi[0:-1])/2 #1
    MDrho = (MDrho[1::] + MDrho[0:-1])/2
    time_at_0 = Planck18.lookback_time(0).value * (u.Gyr).to(u.Myr) 
    DelayTime_set = Planck18.lookback_time(Formationz_set).value * (u.Gyr).to(u.Myr)  - time_at_0
    Formationz_set = (Formationz_set[1::] + Formationz_set[0:-1])/2
    return np.sum(phi * MDrho / (1 + Formationz_set) / np.log(DelayTime_set[-1]/MinDelayTime) * 2 / (DelayTime_set[1::] + DelayTime_set[0:-1]) * (DelayTime_set[1::] - DelayTime_set[0:-1])) 
    

def MergeRate(z, MinDelayTime, frac_meta):
    #merge rate at z, R0 is from ligo's constrain in unit of M_sun/Mpc^{3}/yr
    #in source frame (no 1/(1+z))
    #e(z) Xu 2022 ApJ (B1, B2)
    MinFormationz = FindLookBackRedshift(z, MinDelayTime)
    Formationz_set = np.exp(np.linspace(np.log(MinFormationz), np.log(1100),1000))
    phi = MetallicityCut(frac_meta, Formationz_set)
    MDrho = MDSFR(Formationz_set) 
    phi = 1#(phi[1::] + phi[0:-1])/2 #1
    MDrho = (MDrho[1::] + MDrho[0:-1])/2
    time_at_z = Planck18.lookback_time(z).value * (u.Gyr).to(u.Myr) 
    DelayTime_set = Planck18.lookback_time(Formationz_set).value * (u.Gyr).to(u.Myr)  - time_at_z 
    Formationz_set = (Formationz_set[1::] + Formationz_set[0:-1])/2
    return np.sum(phi * MDrho / (1 + Formationz_set) / np.log(DelayTime_set[-1]/MinDelayTime) * 2 / (DelayTime_set[1::] + DelayTime_set[0:-1]) * (DelayTime_set[1::] - DelayTime_set[0:-1]))  / MergeRateZero(MinDelayTime, frac_meta)
    
    

def DifferentialMergeRateMD14(z, MinDelayTime, frac_meta):
    R0 = 65 / 10 ** 9
    ComovingDis = Planck18.comoving_distance(z).value
    Hz = Planck18.H(z).value
    ckm = c.value * (u.m/u.s).to(u.km/u.s) #千米每秒单位下的光速
    return R0 * MergeRate(z, MinDelayTime, frac_meta) * 4 * np.pi * ckm * ComovingDis ** 2 / Hz

z_test = np.arange(0, 15, 0.1)
dNdz_test = []
for i in tqdm(range(len(z_test))):
    dNdz_test.append(DifferentialMergeRateMD14(z_test[i], 50, 0.3))


plt.semilogy(z_test, dNdz_test, label = 'MD14, $\Delta t_\mathrm{min} = 50 \mathrm{Myr}$')
plt.grid()
plt.ylim(1, 10 ** 5)
plt.xlabel('z')
plt.ylabel('dN/dz')
plt.savefig('./IntermediaPlot/dN_dt.png', dpi=450)

TotalNum = np.sum(dNdz_test) * 0.1

def RedshiftMD14(z):
    #假设最大的红移是14，这样的话可以和haris2018一样。
    #z_set = np.linspace(0, 14, number)
    # test = DifferentialMergeRateMD14(z_set)
    func_inter = sci_inter.interp1d(z_test, dNdz_test / TotalNum)
    return func_inter(z)

PDF_redshift = RedshiftMD14(z_test)
plt.semilogy(z_test, PDF_redshift, label = 'PDF')
plt.semilogy(z_test, dNdz_test / TotalNum, '--', label = 'dNdz/Total')
plt.grid()
plt.legend()
plt.xlabel('z')
plt.ylabel('PDF')
plt.savefig('./IntermediaPlot/PDF_z_test.png', dpi=450)


# #=========================
# #More & More 2021 MNRAS ref: Oguri 2018
# def MoreR(z):
#     R_0 = 22 / 10 ** 9 # /Mpc^3/yr
#     a_1 = 6.6 * 10 ** 3
#     a_2 = 1.6
#     a_3 = 2.1
#     a_4 = 30
#     return R_0 * (a_1 * np.exp(a_2 * z)) / (np.exp(a_3 * z) + a_4)

# def MoredNdz(z):
    
#     ComovingDis = Planck18.comoving_distance(z).value
#     Hz = Planck18.H(z).value
#     ckm = c.value * (u.m/u.s).to(u.km/u.s) #千米每秒单位下的光速
#     return MoreR(z) * 4 * np.pi * ckm * ComovingDis ** 2 / Hz
    


# Moretest = MoredNdz(z_test)
# plt.semilogy(z_test, dNdz_test, 'b')
# plt.semilogy(z_test, Moretest, 'r')
# plt.xlabel(r'$z_s$')
# plt.ylabel(r'$\frac{\mathrm{d}N}{\mathrm{d}z\mathrm{d}t}$')
# plt.xlim(0,20)
# plt.grid()

# #=======================



l_peak = 0.1 #
alpha = 2.63 #
m_min = 4.59 #
delta_m = 4.82 #
m_max = 86.22 #
mu_m = 33.07 #
sigma_m = 5.69 #
beta_q = 1.26 #
def S(m):
    if m < m_min:
        return 0
    elif m >= m_min and m < m_min + delta_m:
        return (np.exp(delta_m/(m - m_min) + delta_m/(m - m_min - delta_m))+1)**(-1)
    else:
        return 1
    
def Model2Pm1(m1):

    B = m1**(-alpha)*(1-alpha)/(m_max**(1-alpha) - m_min**(1-alpha))
    G = 1/(np.sqrt(2*np.pi)*sigma_m)*np.exp(-(m1-mu_m)**2/2/sigma_m**2)
    
    P_m1 = ((1 - l_peak)*B + l_peak*G)*S(m1)
    return P_m1

def Model2Pm2(m1, m2):
    P_m2 = (m2/m1)**beta_q*S(m2)
    return P_m2



def Spin1(chi1):
    #0~0.99之间均匀分布
    return 1/0.99

def Spin2(chin2):
    #0～0.99之间均匀分布
    return 1/0.99

#inclination和polarization取这个形式的分布，
# 只是为了让角动量在空间中均匀分布，也就是按照面积大小来决定概率。
def Inclination(inc):
    #定义域0～pi
    return np.sin(inc)/2

def Polarization(psi):
    #定义域0～pi
    return 1/np.pi

def RightAscension(ra):
    #定义域0～2pi
    return 1/2/np.pi

def Declination(dec):
    #定义域-pi/2~pi/2
    return np.cos(dec)/2

def GPSTime(time):
    #从10**9开始，间隔一年31536000，均匀分布。
    return 1/31536000
    
    

OutPut = [[] for i in range(10)]
number = int(TotalNum)
OutPut[0] = sample(RedshiftMD14, number, lower_bd=0.00001, upper_bd=14.9, guess=2, chebyshev = True)
OutPut[1] = sample(Model2Pm1, number, lower_bd=m_min + 0.1, upper_bd=m_max - 1, guess=6)
for i in range(len(OutPut[1])):
    # if 0.1 * OutPut[1][i] > m_min + 0.05:
    #     lowlim = 0.1 * OutPut[1][i]
    # else:
    #     lowlim = m_min + 0.05
    lowlim = m_min + 0.01
    tmp = sample(Model2Pm2, 1, lower_bd=lowlim, upper_bd=OutPut[1][i], guess=(lowlim + OutPut[1][i])/2, args=(OutPut[1][i],))[0]
    while(tmp<lowlim or tmp > OutPut[1][i]):
        tmp = sample(Model2Pm2, 1, lower_bd=lowlim, upper_bd=OutPut[1][i], guess=(lowlim + OutPut[1][i])/2, args=(OutPut[1][i],))[0]
    OutPut[2].append(tmp)
OutPut[3] = sample(Spin1, number, lower_bd=0, upper_bd=0.99, guess=0.5) 
OutPut[4] = sample(Spin2, number, lower_bd=0, upper_bd=0.99, guess=0.5)
OutPut[5] = sample(Inclination, number, lower_bd=0, upper_bd=np.pi, guess=np.pi/2)
OutPut[6] = sample(Polarization, number, lower_bd=0, upper_bd=np.pi, guess=np.pi/2)
OutPut[7] = sample(RightAscension, number, lower_bd=0, upper_bd=2*np.pi, guess=np.pi/2)
OutPut[8] = sample(Declination, number, lower_bd=-np.pi/2, upper_bd=np.pi/2, guess=0)
OutPut[9] = [random.uniform(10**9, 10**9 + 31536000) for i in range(number)]


np.savetxt('./SampleResult/SampleParameter.csv', OutPut, delimiter=',')

GW_distribution = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')


#PDF
fig, ax = plt.subplots(2, 5, figsize=(12, 5.5))

ax[0,0].hist(GW_distribution[0], density=True, bins=100, label='z')
ax[0,0].semilogy(z_test, dNdz_test / TotalNum)
ax[0,0].legend()


ax[0,1].hist(GW_distribution[1], density=True, bins=100, label='$m_1$')
m1_test = np.arange(m_min, m_max, 0.1)
pdf_m1_test = []
for i in range(len(m1_test)):
    pdf_m1_test.append(Model2Pm1(m1_test[i]))
pdf_m1_test = np.array(pdf_m1_test) / np.sum(pdf_m1_test) / 0.1
ax[0,1].plot(m1_test, pdf_m1_test)
ax[0,1].legend()
ax[0,1].set_xlim(0,100)

ax[0,2].hist(GW_distribution[2], density=True, bins=100, label='$m_2$')
ax[0,2].legend()

ax[0,3].hist(GW_distribution[3], density=True, bins=100, label='$a_1$')
ax[0,3].legend()

ax[0,4].hist(GW_distribution[4], density=True, bins=100, label='$a_2$')
ax[0,4].legend()

ax[1,0].hist(GW_distribution[5], density=True, bins=100, label='$\\theta_{jn}$')
inc_test = np.arange(0, np.pi, 0.1)
ax[1,0].plot(inc_test, Inclination(inc_test))
ax[1,0].legend()


ax[1,1].hist(GW_distribution[6], density=True, bins=100, label='$\\Psi$')
ax[1,1].legend()

ax[1,2].hist(GW_distribution[7], density=True, bins=100, label='RA')
ax[1,2].legend()

ax[1,3].hist(GW_distribution[8], density=True, bins=100, label='DEC')
ax[1,3].legend()

ax[1,4].hist(GW_distribution[9], density=True, bins=100, label='$t_s$')
ax[1,4].legend()


plt.savefig('./IntermediaPlot/GW_sample_result.png', dpi=450)