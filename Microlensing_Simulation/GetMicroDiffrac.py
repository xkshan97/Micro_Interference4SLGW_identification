import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import scipy.fft
from numpy import polyfit, poly1d
from tqdm import tqdm
import struct
from  multiprocessing import Process,Pool
import gc

class DiffracMicro(object):
    def __init__(self):
        self.M_sun = 2 * 10**30

        self.G = 6.67*10**(-11)

        self.c = 3*10**8
        

    def Minimum(self, M_L, z_L, kappa, gamma, fileArea, fileTime, binlength, x_step, freq):

        mu = 1/((1-kappa)**2 - gamma**2)**0.5 #sqrt(mu)
        coeffi = 4 * self.G * M_L * self.M_sun * (1 + z_L) / self.c**3
        constant = 2*np.pi*mu/coeffi
        cw = coeffi/2/np.pi
        
        f1=open(fileArea,"rb")
        Area=struct.unpack("l"*binlength, f1.read(8*binlength))
        f1.close()
        Area = np.array(Area)

        f2=open(fileTime,"rb")
        time=struct.unpack("d"*binlength, f2.read(8*binlength))
        f2.close()
        time = np.array(time)

        index_none_zero = np.where(Area>0)[0][0]
        time = time[index_none_zero::]
        Area = Area[index_none_zero::]
        
    
        Ft_raw = Area[0:-1] * x_step**2 / (time[1::] - time[0:-1])
        time_raw = time[0:-1]
    
        time_raw = time_raw - time_raw[0]
        

        delta_time_raw = time_raw[1] - time_raw[0]
        
        # plt.figure(100)
        plt.semilogx(time_raw,Ft_raw)
        plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
        # plt.xlabel('time')
        # plt.ylabel('dA/dt')
        
        # # plt.xlim(10**(-2),10**(-1))
        # plt.legend()
        
        
        """去掉尾巴"""
        length = len(time_raw)
        time_raw = time_raw[0:length*7//10]
        Ft_raw = Ft_raw[0:length*7//10]
        


        """without apodization"""
        Ft_raw[1::] = Ft_raw[1::] - constant
        Ft_raw[0] = Ft_raw[0] - constant/2
        
        
        
        
        """以上"""

        time_new = time_raw
        Ft_new = Ft_raw
        del time_raw, Ft_raw
        gc.collect()
        #保存时域
        savetmp1 = fileArea.split('Read')[0]
        savetmp2 = savetmp1.split('/')[-1]
        np.savetxt('/disk1/home/shanxk/work/Paper4_MicroLens_Modify/FtnewMinimum_precessing/Ftnew_' + savetmp2 + '.csv', Ft_new, delimiter=',')
        np.savetxt('/disk1/home/shanxk/work/Paper4_MicroLens_Modify/timenewMinimum_precessing/timenew_' + savetmp2 + '.csv', time_new, delimiter=',')
        """加窗归零"""
        lengthnew = len(time_new)
        window = np.hanning(lengthnew*2//10)
        try:
            Ft_new[lengthnew*9//10+1::] = Ft_new[lengthnew*9//10+1::] * window[lengthnew*1//10::] 
        except ValueError:
            Ft_new[lengthnew*9//10::] = Ft_new[lengthnew*9//10::] * window[lengthnew*1//10::] 
            
        plt.semilogx(time_new, Ft_new,label='$F(t)-F_{smooth}(t)$')
        plt.xlabel('time')
        plt.ylabel('F(t)')
        plt.legend()
        plt.grid()
        plt.savefig('/disk1/home/shanxk/work/Paper4_CE_Modify/FigMicro/FigFt_Minimum/Ft_' + fileArea.split('/')[-1].split('.')[0] + '_precessing.png', dpi=450)
        plt.close()
            
        """自己写傅里叶变换"""
        omegafreq = freq * 2 * np.pi
        Freal = []
        Fimag = []
        for i in range(len(omegafreq)):
            Freal.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
            Fimag.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
            
        Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
        Ff = Ff* omegafreq / complex(0,1)*cw
        Ffsgn = constant * cw 
        Ff = Ff + Ffsgn
        del Ft_new, time_new
        gc.collect()
    
        # """以上"""
        Ffabs = np.abs(Ff)
        ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
        # plt.figure(1)
        
        plt.semilogx(freq, Ffabs,label = 'Numerical Result') #//FIXME
    
        plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
    
        plt.grid()
        plt.legend()
    
        plt.xlabel('f[Hz]')
        plt.ylabel('|F(f)|')
        plt.savefig('/disk1/home/shanxk/work/Paper4_CE_Modify/FigMicro/FigFf_Minimum/Ff_' + fileArea.split('/')[-1].split('.')[0] + '_precessing.png', dpi=450)
        plt.close()
        # """相位"""
        # plt.figure(2)
        
        plt.semilogx(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME
        plt.plot([freq[0], freq[-1]], [0, 0], label='Macro magnification')
        plt.grid()
        plt.legend()
        plt.xlabel('f[Hz]')
        plt.ylabel('$\\theta_F$')
        plt.savefig('/disk1/home/shanxk/work/Paper4_CE_Modify/FigMicro/FigFf_Minimum/ThetaF_' + fileArea.split('/')[-1].split('.')[0] + '_precessing.png', dpi=450)
        plt.close()
        
        return Ffabs, ThetaF
    
    def Saddle(self, M_L, z_L, kappa, gamma, fileArea, fileTime, binlength, x_step, x10, x20, freq):
        
        mu = abs(1/((1-kappa)**2 - gamma**2))**0.5
        coeffi = 4 * self.G * M_L * self.M_sun * (1 + z_L) / self.c**3
        constant = mu/coeffi
        cw = coeffi/2/np.pi
        
        f1=open(fileArea, "rb")
        Area=struct.unpack("l"*binlength, f1.read(8*binlength))
        f1.close()
        Area = np.array(Area)

        f2=open(fileTime, "rb")
        time=struct.unpack("d"*binlength, f2.read(8*binlength))
        f2.close()
        time = np.array(time)
        

        Ft_raw = Area[0:-1] * x_step**2 / (time[1::] - time[0:-1]) #根据像素点数得到dS/dt
        time_raw = time[0:-1]
    
        delta_time_raw = np.diff(time_raw)[0]

        # plt.plot(time_raw,Ft_raw)
        # # plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant])
        # plt.xlabel('time')
        # plt.ylabel('dA/dt')
        # plt.grid()
        # plt.savefig('./test.png')
       
        #防止出现t=0，因为ln|t|。
        index0 = np.where(time_raw ==0)[0]
        if len(index0) != 0:
            print(index0)
            time_raw = (time_raw[1::] + time_raw[0:-1])/2
            Ft_raw = Ft_raw[0:-1]
            plt.plot(time_raw, Ft_raw)
        else:
            pass
        
        
        length = len(time_raw)
        time_raw = time_raw[length*5//20:length//20*15] #去掉尾巴处的下降部分
        Ft_raw = Ft_raw[length*5//20:length//20*15]
    
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
        
        del time_raw, Ft_raw, Ft_theory
        gc.collect()
        
        plt.subplot(211)
        plt.plot(time_new, Ft_subtract)
        # plt.xlabel('time')
        plt.ylabel('F_subtract(t)')
        plt.grid()
        plt.subplot(212)
        plt.semilogx(time_new, Ft_subtract)
        plt.xlabel('time')
        plt.ylabel('F_subtract(t)')
        # plt.legend()
        plt.grid()
        plt.savefig('/disk1/home/shanxk/work/Paper4_CE_Modify/FigMicro/FigFt_Saddle/Ft_' + fileArea.split('/')[-1].split('.')[0] + '_precessing.png', dpi=450)
        plt.close()
        
        #保存时域
        savetmp1 = fileArea.split('Read')[0]
        savetmp2 = savetmp1.split('/')[-1]
        np.savetxt('/disk1/home/shanxk/work/Paper4_MicroLens_Modify/FtnewSaddle_precessing/Ftnew_' + savetmp2 + '.csv', Ft_subtract, delimiter=',')
        np.savetxt('/disk1/home/shanxk/work/Paper4_MicroLens_Modify/timenewSaddle_precessing/timenew_' + savetmp2 + '.csv', time_new, delimiter=',')
        """加窗归零"""
        lengthnew = len(time_new)
        window = np.hanning(lengthnew*2//10)
        try:
            Ft_subtract[0:lengthnew*1//10] = Ft_subtract[0:lengthnew*1//10] * window[0:lengthnew*1//10]
            Ft_subtract[lengthnew*9//10+1::] = Ft_subtract[lengthnew*9//10+1::] * window[lengthnew*1//10::] 
        except ValueError:
            Ft_subtract[0:lengthnew*1//10] = Ft_subtract[0:lengthnew*1//10] * window[0:lengthnew*1//10]
            Ft_subtract[lengthnew*9//10::] = Ft_subtract[lengthnew*9//10::] * window[lengthnew*1//10::] 
            


        """自己写傅里叶变换"""
        omegafreq = 2 * np.pi * freq
        Freal = []
        Fimag = []
        for i in range(len(omegafreq)):
            Freal.append(np.sum(Ft_subtract*np.cos(omegafreq[i]*time_new))*delta_time_raw)
            Fimag.append(np.sum(Ft_subtract*np.sin(omegafreq[i]*time_new))*delta_time_raw)
            
        Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
        Ff = Ff* omegafreq / complex(0,1)*cw
        Ffsgn = - np.sign(omegafreq)* 2*np.pi*constant * cw * complex(0,1)
        Ff = Ff + Ffsgn
        
        del Ft_subtract, time_new
        gc.collect()


        Ffabs = np.abs(Ff) #取模   
        ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs)) #计算相位
        
        plt.semilogx(freq, Ffabs,label = 'Numerical Result') #//FIXME
    
        plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
    
        plt.grid()
        plt.legend()
    
        plt.xlabel('f[Hz]')
        plt.ylabel('|F(f)|')
        plt.savefig('/disk1/home/shanxk/work/Paper4_CE_Modify/FigMicro/FigFf_Saddle/Ff_' + fileArea.split('/')[-1].split('.')[0] + '_precessing.png', dpi=450)
        plt.close()
        # """相位"""
        # plt.figure(2)
        
        plt.semilogx(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME
        plt.plot([freq[0], freq[-1]], [-np.pi/2,-np.pi/2], label='macro magnification')
        plt.grid()
        plt.legend()
        plt.xlabel('f[Hz]')
        plt.ylabel('$\\theta_F$')
        plt.savefig('/disk1/home/shanxk/work/Paper4_CE_Modify/FigMicro/FigFf_Saddle/ThetaF_' + fileArea.split('/')[-1].split('.')[0] + '_precessing.png', dpi=450)
        plt.close()
        
        return Ffabs, ThetaF
    
    def Maximum(self, M_L, z_L, kappa, gamma, fileArea, fileTime, binlength, x_step, freq):
       
        mu = 1/((1-kappa)**2 - gamma**2)**0.5
        coeffi = 4 * self.G * M_L * self.M_sun * (1 + z_L) / self.c**3
        constant = 2*np.pi*mu/coeffi
        cw = coeffi/2/np.pi
        
        f1=open(fileArea,"rb")
        Area=struct.unpack("l"*binlength, f1.read(8*binlength))
        f1.close()
        Area = np.array(Area)

        f2=open(fileTime,"rb")
        time=struct.unpack("d"*binlength, f2.read(8*binlength))
        f2.close()
        time = np.array(time)
        
        Ft_raw = Area[0:-1] * x_step**2 / (time[1::] - time[0:-1])
        time_raw = time[0:-1]
      
        time_raw = time_raw - time_raw[-1]

        delta_time_raw = time_raw[1] - time_raw[0]
        
        plt.semilogx(-time_raw,Ft_raw)
        # plt.plot([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
        plt.xlabel('time')
        plt.ylabel('dA/dt')
        # plt.xlim(10**(-2),10**(-1))
        plt.legend()
        
        """去掉尾巴"""
        length = len(time_raw)
        time_raw = time_raw[length*1//5::]
        Ft_raw = Ft_raw[length*1//5::]
      

        """without apodization"""
        Ft_raw[0:-1] = Ft_raw[0:-1] - constant
        Ft_raw[-1] = Ft_raw[-1] - constant/2
        
        
        time_new = time_raw
        Ft_new = Ft_raw
        """加窗"""
        lengthnew = len(time_new)
        window = np.hanning(lengthnew*2//10)
        try:
            Ft_new[0:lengthnew*1//10] = Ft_new[0:lengthnew*1//10] * window[0:lengthnew*1//10]
            
        except ValueError:
            Ft_new[0:lengthnew*1//10] = Ft_new[0:lengthnew*1//10] * window[0:lengthnew*1//10]
            

        
        """自己写傅里叶变换"""
        omegafreq = 2*np.pi*freq
        Freal = []
        Fimag = []
        for i in range(len(omegafreq)):
            Freal.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
            Fimag.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
            
        Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
        Ff = Ff* omegafreq / complex(0,1)*cw
        Ffsgn = -constant * cw 
        Ff = Ff + Ffsgn
        
        Ffabs = np.abs(Ff)
        ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
        
        return Ffabs, ThetaF