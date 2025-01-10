import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import matplotlib.pyplot as pp
import pycbc.noise
import pycbc.psd
import pylab
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import Detector
from astropy.cosmology import Planck18
from tqdm import tqdm
from pycbc import frame
import pycbc.filter
from pycbc.filter import highpass_fir, matched_filter
from pycbc.psd import welch, interpolate
from  multiprocessing import Process,Pool
from scipy.interpolate import interp1d
import scipy
import struct
import sys
sys.path.append(r"../Paper4_MicroLens_CE/")
from GetMicroDiffrac import DiffracMicro

# Generate a PSD using an analytic expression for 
# the full design Advanced LIGO noise curve
f_lower = 30
duration = 4096
sample_rate = 4096
lens_freq = np.arange(-sample_rate//2, sample_rate//2+1, 1)
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)


name_sample = '1YearFeiXuMD14'
number = 37360
multidetector = True
if multidetector:
    name = '1YearFeiXuMD14_multi'
else:
    name = '1YearFeiXuMD14_ligo' 
SampleParam = np.loadtxt('./SampleResult/SampleParameter' + name_sample + str(number) + '.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

apx = 'IMRPhenomPv2'


# Now, let's generate noise that has the same spectrum
htilde_l1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=1000)
htilde_h1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=2000)

AcceptLensIndex = np.loadtxt('./SampleResult/AcceptLensIndex' + name + str(number) + '.csv', delimiter=',')

magnification = np.loadtxt('./SampleResult/magnification.csv', delimiter=',')
timedelay = np.loadtxt('./SampleResult/timedelayDays.csv', delimiter=',')
kappa = np.loadtxt('./SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('./SampleResult/gamma.csv', delimiter=',')
imagenum = np.loadtxt('./SampleResult/imagenumber.csv', delimiter=',')
lens_z_set = np.loadtxt('./SampleResult/lens_z.csv', delimiter=',')


SNR_network = np.loadtxt('./Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')
Index_SNR_gtr_12 = np.where(SNR_network>12)[0]



num_min = 10
num_sad = 10
f1=open('../Paper4_MicroLens_Modify/ResultMinimum/TimeLength.bin',"rb")
TimeLengthFileMinimum=struct.unpack("l"*num_min, f1.read(8*num_min))
f1.close()
TimeLengthFileMinimum = np.array(TimeLengthFileMinimum)

f2=open('../Paper4_MicroLens_Modify/ResultSaddle/TimeLength.bin',"rb")
TimeLengthFileSaddle=struct.unpack("l"*num_sad, f2.read(8*num_sad))
f2.close()
TimeLengthFileSaddle = np.array(TimeLengthFileSaddle)

f3=open('../Paper4_MicroLens_Modify/ResultSaddle/X10.bin',"rb")
X10File=struct.unpack("d"*num_sad, f3.read(8*num_sad))
f3.close()
X10File = np.array(X10File)

f4=open('../Paper4_MicroLens_Modify/ResultSaddle/X20.bin',"rb")
X20File=struct.unpack("d"*num_sad, f4.read(8*num_sad))
f4.close()
X20File = np.array(X20File)