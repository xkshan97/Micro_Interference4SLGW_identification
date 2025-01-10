'''
Author: your name
Date: 2021-10-16 16:04:33
LastEditTime: 2021-10-22 17:34:46
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /code/mock_/fromC2py.py
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import scipy.fft
from numpy import polyfit, poly1d
from tqdm import tqdm
import struct
from  multiprocessing import Process,Pool


with open('./root_2.txt', 'r') as f0:
    data = f0.readlines()
    
fermat_pot = []
mag_geo = []
for i in range(len(data)):
    fermat_pot.append(float(data[i].split(' ')[12]))
    mag_geo.append(float(data[i].split(' ')[8]))
    

AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')

magnification = np.loadtxt('../Paper4_CE_Modify/SampleResult/magnification.csv', delimiter=',')
timedelay = np.loadtxt('../Paper4_CE_Modify/SampleResult/timedelayDays.csv', delimiter=',')
kappa = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('../Paper4_CE_Modify/SampleResult/gamma.csv', delimiter=',')
imagenum = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
lens_z_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')

i_ = 9
i = int(AcceptLensIndex[i_]) #在SamplePama中的位置。
imag_num = int(imagenum[i_])
save_index = 0
IndexSumimage = int(np.sum(imagenum[0:i_]))
z_L = lens_z_set[i_]

M_L = 0.38
M_sun = 2 * 10**30
G = 6.67 * 10**(-11)
c = 3 * 10**8
mu = 1/((1-kappa[IndexSumimage])**2 - gamma[IndexSumimage]**2)**0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
constant = 2*np.pi*mu/coeffi
cw = coeffi/2/np.pi

yi = 0

f1=open("./ResultMinimum_manual/Total" + str(i) + "_" + str(save_index) + "ReadAreaMinimum.bin","rb")
Area=struct.unpack("l"*1167458, f1.read(8*1167458))
f1.close()
Area = np.array(Area)

f2=open("./ResultMinimum_manual/Total" + str(i) + "_" + str(save_index) + "ReadTimeMinimum.bin","rb")
time=struct.unpack("d"*1167458, f2.read(8*1167458))
f2.close()
time = np.array(time)


x_step = 0.005
x1 = np.arange(-81.6227766017,81.6227766017, x_step)
# x1 = np.arange(-30.05,30.05,x_step)
x2 = x1
Ft_raw = Area[0:-1] * (x1[-1] - x1[0])**2 / len(x1)**2 / (time[1::] - time[0:-1])
time_raw = time[0:-1]

index_tmp = np.where(Ft_raw > 0)[0][0]
Ft_raw = Ft_raw[index_tmp:]
time_raw = time_raw[index_tmp:]

time_geo = coeffi * np.array(fermat_pot) - time_raw[0]
# time_geo = time_geo[1::]
time_raw = time_raw - time_raw[0]

inter_time_raw = interp1d(time_raw, Ft_raw)

delta_time_raw = time_raw[1] - time_raw[0]
index_max = np.where(Ft_raw == np.max(Ft_raw))[0][0]
time_max = time_raw[index_max]
time_geo = time_geo - time_geo[2] + time_max

plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
plt.semilogx(time_raw[0:len(time_raw)*4//5],Ft_raw[0:len(time_raw)*4//5], label='micro')
plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant], '--',label='macro')
plt.scatter(time_geo[1::], inter_time_raw(time_geo[1::]), label='geometric image', color='grey')
plt.xlabel('$t$[s]')
plt.ylabel('$F(t)$')
plt.grid()
plt.legend()
plt.xlim(2.5*10**(-5),1)
plt.ylim(0.8 * 10**6, 1.5 * 10 **6)
plt.savefig('to_referee/Ft_png', dpi=450)



"""去掉尾巴"""
length = len(time_raw)
time_raw = time_raw[0:length*4//5]
Ft_raw = Ft_raw[0:length*4//5]


"""以上"""

# """设置保护区"""
# time_raw_2 = np.arange(time_raw[-1] + delta_time_raw, time_raw[-1] + 1,delta_time_raw)
# Ft_raw_2 = np.array([constant]*len(time_raw_2))
# time_raw = np.append(time_raw, time_raw_2)
# Ft_raw = np.append(Ft_raw, Ft_raw_2)

plt.semilogx(time_raw, Ft_raw)
plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant])
plt.grid()


"""without apodization"""
Ft_raw[1::] = Ft_raw[1::] - constant
Ft_raw[0] = Ft_raw[0] - constant/2



plt.semilogx(time_raw, Ft_raw,label='$F(t)-F_{smooth}(t)$')
plt.xlabel('time')
plt.ylabel('F(t)')
plt.legend()
plt.grid()
"""以上"""

time_new = time_raw
Ft_new = Ft_raw
"""加窗归零"""
lengthnew = len(time_new)
window = np.hanning(lengthnew*2//5)
try:
    Ft_new[lengthnew*4//5+1::] = Ft_new[lengthnew*4//5+1::] * window[lengthnew*1//5::] 
except ValueError:
    Ft_new[lengthnew*4//5::] = Ft_new[lengthnew*4//5::] * window[lengthnew*1//5::] 
    

np.savetxt('./timenewMinimum_manual/time_new_Total' + str(i) + "_" + str(save_index) + '.csv',time_new,delimiter=',')
np.savetxt('./timenewMinimum_manual/Ft_new_Total' + str(i) + "_" + str(save_index) + '.csv',Ft_new,delimiter=',')
    

#加窗归零
# lengthnew = len(time_new)
# window = np.hanning(lengthnew*2//10)
# try:
    
#     Ft_new[lengthnew*9//10::] = Ft_new[lengthnew*9//10::] * window[lengthnew*1//10::] 
# except ValueError:
#     Ft_new[lengthnew*9//10+1::] = Ft_new[lengthnew*9//10+1::] * window[lengthnew*1//10+1::]  
    

plt.semilogx(time_new, Ft_new)
# # plt.ylim(-8000,2000)
# plt.grid()
    
"""自己写傅里叶变换"""
omegafreq = np.arange(0.1*2*np.pi,500*2*np.pi, 0.1 * 2*np.pi)
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



"""以上"""


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
    


plt.figure(1)

plt.semilogx(freq, Ffabs,label = 'Diffraction Result') #//FIXME
plt.semilogx(freq, FfTheoryAbs, label='Geo Approximation')
plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro Magnification')

plt.grid()
plt.legend()
plt.ylim(1,3)
plt.xlabel('$f$[Hz]')
plt.ylabel('$|F(f)|$')
plt.xlim(0.1, 850)
plt.savefig("./to_referee/Ff.png",dpi=450)
"""相位"""
plt.figure(2)

plt.semilogx(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME
plt.semilogx(freq, FfTheoryTheta, label = 'geo approximation')
# plt.loglog(test1,test2)
# plt.ylim(-0.2,0.2)
plt.grid()
plt.legend()
plt.xlabel('f[Hz]')
plt.ylabel(r'$\theta_F$')


np.savetxt('./ResultMinimum/Ffabs'+file+'.csv',Ffabs,delimiter=',')
np.savetxt('./ResultMinimum/ThetaF'+file + '.csv',ThetaF,delimiter=',')
