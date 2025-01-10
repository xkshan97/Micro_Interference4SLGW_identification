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
from scipy.optimize import curve_fit
import struct
from  multiprocessing import Process,Pool


M_L = 0.38
z_L = 1.271432154745443688e+00
M_sun = 2 * 10**30
G = 6.67 * 10**(-11)
c = 3 * 10**8
kappa = 5.984599014269287309e+00
gamma = 5.984599014269287309e+00
mu = abs(1/((1-kappa)**2 - gamma**2))**0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
constant = mu/coeffi
cw = coeffi/2/np.pi
x_step = 0.005 #像平面分辨率
delta_time_raw = 10**(-6) #时间分辨率




#60 边界




f1=open("./ResultSaddle/Total4421_1ReadAreaSaddle.bin","rb")
Area=struct.unpack("l"*13980387, f1.read(8*13980387))
f1.close()
Area = np.array(Area)

f2=open("./ResultSaddle/Total4421_1ReadTimeSaddle.bin","rb")
time=struct.unpack("d"*13980387, f2.read(8*13980387))
f2.close()
time = np.array(time)

# index_ = np.where(Area>0)
# time = time[index_]
# Area = Area[index_]

Ft_raw = Area[0:-1] * x_step**2 / (time[1::] - time[0:-1]) #根据像素点数得到dS/dt
time_raw = time[0:-1]
# time_raw = time_raw - time_raw[np.where(Ft_raw==np.max(Ft_raw))]
delta_time_raw = np.diff(time_raw)[0]

plt.plot(time_raw,Ft_raw)
# plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant])
plt.xlabel('time')
plt.ylabel('dA/dt')
plt.grid()


index0 = np.where(time_raw ==0)[0]
if len(index0) != 0:
    print(index0)
    time_raw = (time_raw[1::] + time_raw[0:-1])/2
    Ft_raw = Ft_raw[0:-1]
    plt.plot(time_raw, Ft_raw)
else:
    pass


length = len(time_raw)
time_raw = time_raw[length//10:length//10*9] #去掉尾巴处的下降部分
Ft_raw = Ft_raw[length//10:length//10*9]


x10 = 429.891
x20 = 129.799
mur = 1 - kappa + gamma
mut = kappa + gamma - 1
Axisa = np.sqrt(1/coeffi/mur)
Axisb = np.sqrt(1/coeffi/mut)
index1 = np.where(time_raw < 0)
index2 = np.where(time_raw >= 0)
F1 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisa) + np.log(np.abs(time_raw[index1])) - 2 * np.log(x10 + np.sqrt(2 * Axisa**2 * np.abs(time_raw[index1]) + x10**2)))
F2 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisb) + np.log(np.abs(time_raw[index2])) - 2 * np.log(x20 + np.sqrt(2 * Axisb**2 * np.abs(time_raw[index2]) + x20**2))) 
Ft_theory = np.append(F1, F2)

Ft_subtract = Ft_raw - Ft_theory

time_new = time_raw
"""加窗归零"""
lengthnew = len(time_new)
window = np.hanning(lengthnew*4//10)
try:
    Ft_subtract[0:lengthnew*2//10] = Ft_subtract[0:lengthnew*2//10] * window[0:lengthnew*2//10]
    Ft_subtract[lengthnew*8//10+1::] = Ft_subtract[lengthnew*8//10+1::] * window[lengthnew*2//10::] 
except ValueError:
    Ft_subtract[0:lengthnew*2//10] = Ft_subtract[0:lengthnew*2//10] * window[0:lengthnew*2//10]
    Ft_subtract[lengthnew*8//10::] = Ft_subtract[lengthnew*8//10::] * window[lengthnew*2//10::] 
    
    
np.savetxt('./ResultSaddle/60time_new'+file+'.csv',time_new,delimiter=',')
np.savetxt('./ResultSaddle/60Ft_new'+file+'.csv',Ft_subtract,delimiter=',')




"""自己写傅里叶变换"""
omegafreq = np.arange(0.1*2*np.pi,2000*2*np.pi,0.5*2*np.pi)
freq = omegafreq/2/np.pi
Freal = []
Fimag = []
for i in tqdm(range(len(omegafreq))):
    Freal.append(np.sum(Ft_subtract*np.cos(omegafreq[i]*time_new))*delta_time_raw)
    Fimag.append(np.sum(Ft_subtract*np.sin(omegafreq[i]*time_new))*delta_time_raw)
    
Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
Ff = Ff* omegafreq / complex(0,1)*cw
Ffsgn = - np.sign(omegafreq)* 2*np.pi*constant * cw * complex(0,1)
Ff = Ff + Ffsgn

'''
"""下面做FFT"""
N = len(time_new)
freq = np.linspace(0, 1/(2*(delta_time_raw)), N//2) 
freq = freq[0:N//2:1]
omegafreq = freq * 2 * np.pi #得到圆频率
Ff = scipy.fft.fft(Ft_subtract)

"""符号函数修正"""
Ffsgn = - 2*np.pi*constant * cw * complex(0,1)

# Ffreal = np.real(-Ff[-1:-N//2 - 1 :-1]* delta_time_raw * omegafreq * complex(0,1)*cw)
# Ffimag = np.imag(-Ff[-1:-N//2 - 1 :-1]* delta_time_raw * omegafreq * complex(0,1)*cw) + Ffsgn
Ff = Ff[-1:-N//2-1:-1] * delta_time_raw * omegafreq / complex(0,1)*cw #Ff[-1:-N//2:-1]是为了取傅里叶变换中负频率部分
#上面之所以取Ff的后半部分，是因为从理论公式上来看，要想从dS/dt得到F(f)，如果做的是傅里叶变换，则需要取负频率部分，
# 如果是做逆傅里叶变换，则取正频率部分。
# 同理Ffsgn也是如此
"""用来翻转符号的"""
# Ffabs = np.abs(Ff)
Ff[0:len(Ff):2] = - Ff[0:len(Ff):2] #把奇数位置的频率变负，等价于把相位减去pi
# test1 = np.real(-complex(0,1)*np.log(Ff[2]/Ffabs[2]))
# test2 = np.real(-complex(0,1)*np.log(Ff[3]/Ffabs[3]))
# if(test1>=0 and test2 <=0):
#     Ff[1:len(Ff):2] = - Ff[1:len(Ff):2]
# elif(test2>=0 and test1 <=0):
#     Ff[0:len(Ff):2] = - Ff[0:len(Ff):2]
# else:
#     print(test1, test2)
# Ff[0::2] = -Ff[0::2] """嗯？"""
Ff = Ff + Ffsgn
'''

Ffabs = np.abs(Ff) #取模   
ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs)) #计算相位
# Ffabs = signal.savgol_filter(Ffabs, 9, 1) #这个是把振幅做平滑平均
# ThetaF = signal.savgol_filter(ThetaF, 119, 1)
Ff = Ffabs * np.exp(complex(0,1)*ThetaF) #用平滑后的振幅重组复数Ff
# freq = freq[1::] #去除0频率
# Ff = Ff[1::]
Ffabs = np.abs(Ff)



plt.figure(4)
plt.semilogx(freq, Ffabs, label = 'Sign') #//FIXME
# plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
# plt.plot(freq, Ffabs1, label='Cut')
# plt.plot(freq, FfTheoryAbs, '--', color = 'red',label='Geometric')

# plt.xlim(1,2000)
# plt.loglog(test1,test2)
# plt.ylim(2.5,3.5)
plt.grid()
plt.legend()
plt.xlabel('f[Hz]')
plt.ylabel('|F(f)|')
# plt.savefig("Ffembed.png",dpi=450)

plt.figure(5)
plt.semilogx(freq, ThetaF, label = 'Sign') #//FIXME

# plt.xlim(1,2000)
# plt.loglog(test1,test2)
# plt.ylim(-2.3,-1)
plt.grid()
plt.legend()
plt.xlabel('f[Hz]')
plt.ylabel(r'$\theta_F$')



np.savetxt('./ResultSaddle/60Ffabs'+file+'.csv',Ffabs,delimiter=',')
np.savetxt('./ResultSaddle/60ThetaF'+file+'.csv',ThetaF,delimiter=',')
np.savetxt('./ResultSaddle/60Geoabs'+file+'.csv',FfTheoryAbs,delimiter=',')
np.savetxt('./ResultSaddle/60GeoThetaF'+file+'.csv',FfTheoryTheta,delimiter=',')


