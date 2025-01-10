import numpy as np
import matplotlib.pyplot as pp
import pycbc.noise
import pycbc.psd
from pycbc import frame
import pycbc.filter
from pycbc.psd import welch, interpolate
import pylab

# Generate a PSD using an analytic expression for 
# the full design Advanced LIGO noise curve
f_lower = 10
duration = 4096
sample_rate = 4096
tsamples = sample_rate * duration
fsamples = tsamples // 2 + 1
df = 1.0 / duration
psd_l1_h1 = pycbc.psd.CosmicExplorerP1600143(fsamples, df, f_lower)


htilde_l1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=10000)
htilde_h1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=20000)
htilde_v1 = pycbc.noise.frequency_noise_from_psd(psd_l1_h1, seed=30000)


# Equivelantly in the time domain
hoft_h1 = htilde_h1.to_timeseries()
hoft_l1 = htilde_l1.to_timeseries()
hoft_v1 = htilde_v1.to_timeseries()

"""=================================="""
"""对时域的噪音hoft用welch方法重新估计psd，并与原始的psd对比"""

# Use Welch's method with 4s segments
psd_h1_estimated = interpolate(welch(hoft_h1), 1.0 / hoft_h1.duration)
psd_l1_estimated = interpolate(welch(hoft_l1), 1.0 / hoft_l1.duration)
psd_v1_estimated = interpolate(welch(hoft_v1), 1.0 / hoft_v1.duration)

'''
pylab.loglog(psd_h1_estimated.sample_frequencies, psd_h1_estimated, label='Estimated')
pylab.loglog(psd_l1_h1.sample_frequencies, psd_l1_h1, label='Original')
pylab.xlim(5, 1024)
pylab.ylim(1e-48, 1e-45)

'''
"""================================="""


start_time = 1126259642.413#SampleParam[9][i]
declination = -1.2108#SampleParam[8][i]
right_ascension = 1.375#SampleParam[7][i]
polarization = 2.659#SampleParam[6][i]

hoft_h1.start_time += start_time
hoft_l1.start_time += start_time
hoft_v1.start_time += start_time

hoft_h1.start_time -= 2048
hoft_l1.start_time -= 2048
hoft_v1.start_time -= 2048

#++++++++++++++++++++++
#调节开始时间为整数
hoft_h1.start_time = hoft_h1.start_time.gpsSeconds
hoft_l1.start_time = hoft_l1.start_time.gpsSeconds
hoft_v1.start_time = hoft_v1.start_time.gpsSeconds

#===============================


frame.write_frame('./Sim_PSD_Data/H1/test.gwf','H1',hoft_h1)
frame.write_frame('./Sim_PSD_Data/L1/test.gwf','L1',hoft_l1)
frame.write_frame('./Sim_PSD_Data/V1/test.gwf','V1',hoft_v1)

# pylab.plot(hoft_l1_h1.sample_times, hoft_l1_h1,'b')
# pp.plot(signal_h1.sample_times, signal_h1)
# pp.ylabel('Strain')
# pp.xlabel('Time (s)')
# pp.legend()
# pp.show()