import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import scipy.fft
from numpy import polyfit, poly1d
import struct
from tqdm import tqdm
import pandas as pd
import seaborn as sns
"""
新的数据放在/data/pubdata/Data3_cxc_to_sxk里了:
还是原来的参数和分辨率
lenspos_848.bin :  星的位置,848*2
map_266505625.bin : X->Y映射, 266505625*6, 像平面23231*23231个点的位置(0 1),对应的源的位置(2 3),放大率(4),引力势(5)
mag_for_each_src.bin : 1000000 int, 源平面上1000*1000个源对应的像的数目
Roots目录下data_i=i1_j=j1.bin,  记录编号i=i1,  j=j1的源成像的信息,  (N_image+1)*5,  最后一行记录该源的坐标和总放大率,  即data[N_image][0] data[N_image][1]  data[N_image][2] 
圆圈是5.5倍的像平面分辨率，会丢掉放大率0.01以下的像
"""


'''
#read map
M_L = 1
z_L = 0.5
M_sun = 2 * 10**30
G = 6.67 * 10**(-11)
c = 3 * 10**8
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
npix=20046**2
column = 1
f=open("./Phi.dat","rb")
data=struct.unpack("d"*npix*column, f.read(8*npix*column))
f.close()
# data = np.array(data).reshape(-1,6)
PhiMine = np.array(data)
data = 0

PhiMicroMine = PhiMine 
PhiMicroMine = PhiMicroMine.reshape(-1, 20046)
# Phi = np.flip(Phi,axis=0)
PhiMine = 0




plt.imshow(PhiMicroMine)
plt.colorbar()








#taylor
npix=20046**2
column = 1
f=open("./TaylorPhi.dat","rb")
data=struct.unpack("d"*npix*column, f.read(8*npix*column))
f.close()
# data = np.array(data).reshape(-1,6)
PhiMineTaylor = np.array(data)
data = 0

# print(len(PhiMine))
PhiMineTaylor = PhiMineTaylor.reshape(-1, 20046) 
PhiMicroMineTaylor = PhiMineTaylor - (1/2*0.4*(X1**2 + Y1**2) - 1/2 * 0.4 *(Y1**2 - X1**2))
PhiMineTaylor = 0
plt.imshow(PhiMicroMineTaylor)
plt.colorbar()


DeltaPhiMine = (PhiMicroMineTaylor - PhiMicroMine)/PhiMicroMine
plt.imshow(DeltaPhiMine)
plt.colorbar()
'''
kappa_star_set = [0.01, 0.03, 0.05, 0.07, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 2.00, 4.00, 6.00, 8.00]
for i in range(len(kappa_star_set)):
    kappa_star = '%.2f'%kappa_star_set[len(kappa_star_set) - 1 - i]
    #read Mnimum map
    npix=20000**2
    column = 1
    f=open('stellar_and_remnant_data_Phi/' + kappa_star + ".bin","rb")
    data=struct.unpack("d"*npix*column, f.read(8*npix*column))
    f.close()
    # data = np.array(data).reshape(-1,6)
    # ImagPosx = np.array(data[0::5])
    # ImagPosy = np.array(data[1::5])
    # SourPosx = np.array(data[2::5])
    # SourPosy = np.array(data[3::5])
    PhiMicro = np.array(data)
    del data
    PhiMicro = PhiMicro.reshape(-1, 20000).T
    # ImagPosx = ImagPosx.reshape(-1, 20046)
    # plt.figure(i)
    plt.imshow(PhiMicro)
    plt.colorbar()
    

    plt.savefig('./IntermediaPlot/kappa_star_' + kappa_star + '.png', dpi=450)
    plt.close()
    del PhiMicro
   


# #read Saddle map
# npix=20000 * 20000
# column = 1
# f=open("/data/pubdata/Data_cxc_to_sxk/test.bin","rb")
# data=struct.unpack("d"*npix*column, f.read(8*npix*column))
# f.close()
# # data = np.array(data).reshape(-1,6)
# # ImagPosx = np.array(data[0::5])
# # ImagPosy = np.array(data[1::5])
# # SourPosx = np.array(data[2::5])
# # SourPosy = np.array(data[3::5])
# CutPhiMicro = np.array(data)
# data = 0
# CutPhiMicro = CutPhiMicro.reshape(-1, 20000) 
# # ImagPosx = ImagPosx.reshape(-1, 20046)

# plt.imshow(CutPhiMicro)
# plt.colorbar()
# plt.savefig('Cut1.png', dpi=450)





''''
PhiMicroDis = [[]]*9622
for i in tqdm(range(9622)):
    for j in range(9622):
        PhiMicroDis[i].append(PhiMicro[i][j] - PhiMicroTest[i][j])
    # PhiMicroDis[i] = np.array(PhiMicroDis[i])
PhiMicroDis = np.array(PhiMicroDis)

plt.imshow(PhiMicro, vmin = -42, vmax = 35)
plt.colorbar()
# with open('potentialmap.bin','wb') as fp:
#     for x in Phi:
#         a = struct.pack('d',x)
#         fp.write(a)
    
# fp.close()


DeltaPhi = np.abs((PhiMicro - PhiMicroMine)/PhiMicro)
plt.imshow(np.log(DeltaPhi))
plt.colorbar()

plt.imshow(np.log(np.abs(DeltaPhi)))
plt.colorbar()
'''