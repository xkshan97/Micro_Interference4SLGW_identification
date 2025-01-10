import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import bilby
from gwpy.timeseries import TimeSeries
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import gwpy

sample_RATE = '4KHZ'
sample_rate = 4096

logger = bilby.core.utils.logger
H1_snr = np.loadtxt('./Sim_GW_Data/H1/H1_snr.csv', delimiter=',')
L1_snr = np.loadtxt('./Sim_GW_Data/L1/L1_snr.csv', delimiter=',')
V1_snr = np.loadtxt('./Sim_GW_Data/V1/V1_snr.csv', delimiter=',')

trigger_time_H1 = np.loadtxt('./Sim_GW_Data/trigger_time_H1.csv', delimiter=',')
trigger_time_L1 = np.loadtxt('./Sim_GW_Data/trigger_time_L1.csv', delimiter=',')
trigger_time_V1 = np.loadtxt('./Sim_GW_Data/trigger_time_V1.csv', delimiter=',')



gw_duration = np.loadtxt('./Sim_GW_Data/gw_duration.csv', delimiter=',')

gw_duration_20 = np.loadtxt('./Sim_GW_Data/gw_duration_20.csv', delimiter=',')

dt_h1_from_earth_center = np.loadtxt('./Sim_GW_Data/dt_h1_from_earth_center.csv', delimiter=',')
dt_l1_from_earth_center = np.loadtxt('./Sim_GW_Data/dt_l1_from_earth_center.csv', delimiter=',')      
dt_v1_from_earth_center = np.loadtxt('./Sim_GW_Data/dt_v1_from_earth_center.csv', delimiter=',')      


SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec, trigger_time


#===============================
#得到联合snr大于12的事例（这里查一下，12这个阈值是对于几个探测起来说的）
SNR_network = np.sqrt(H1_snr**2 + L1_snr**2 + V1_snr ** 2)

SNR_select_index = np.loadtxt('../Paper5_cWB_Modify/Data_start_time_and_snr_and_match_and_outdir_name/SNR_select_index_unlens.csv', delimiter=',')
'''
#找到所有采样点中最小的chirp mass, luminosity distance
chirp_mass_set = []
luminosity_distance_set = []
for i in range(len(SampleParam[0][0:300])):
    luminosity_distance = Planck15.luminosity_distance(SampleParam[0][i]).value
    luminosity_distance_set.append(luminosity_distance) 
    #计算一下输入信号的chirp mass，为的是在后面把先验的范围改一下。
    Mz = (SampleParam[0][i] + 1) * (SampleParam[1][i] * SampleParam[2][i]) ** (3/5) / (SampleParam[1][i] + SampleParam[2][i]) ** (1/5) 
    chirp_mass_set.append(Mz)
    
chirp_mass_minimum = np.min(chirp_mass_set)
luminosity_distance_maximum = np.max(luminosity_distance_set)
#在遍历1811个能被探测到的事例以后，最小的chirp mass=6.38, 最大192.8，最大的光度距离为84872，最小1715

np.savetxt('./chirp_mass_and_Dl.csv',[chirp_mass_minimum, luminosity_distance_maximum], delimiter=',')
#=================================#
'''
for i in SNR_select_index[-39:-38]:
    i = int(i)
    #找到文件所对应的inject time
    file_name_list_tmp = os.listdir('./Sim_GW_Data/H1')
    for pp in range(len(file_name_list_tmp)):
        if '-4096_' + str(i) + '.gwf' in file_name_list_tmp[pp]:
            file_name_tmp = file_name_list_tmp[pp]
    Data_start_time = file_name_tmp.split('-')[2]
    #=======================
    #injection parameters
    luminosity_distance = Planck15.luminosity_distance(SampleParam[0][i]).value
    Mz = (SampleParam[0][i] + 1) * (SampleParam[1][i] * SampleParam[2][i]) ** (3/5) / (SampleParam[1][i] + SampleParam[2][i]) ** (1/5)
    injection_parameters = dict(
    mass_ratio = SampleParam[2][i]/SampleParam[1][i], chirp_mass = Mz, a_1=SampleParam[3][i], 
    a_2=SampleParam[4][i], tilt_1=0, tilt_2=0,
    phi_12=0, phi_jl=0, luminosity_distance=luminosity_distance, theta_jn=SampleParam[5][i], psi=SampleParam[6][i],
    phase=0, ra=SampleParam[7][i], dec=SampleParam[8][i], geocent_time = trigger_time_H1[i] - dt_h1_from_earth_center[i])
    #注意geocent_time改成了H1_time。
    #=======================
    
    
    outdir = 'outdir_unlensed/outdir' + str(i)
    label = 'H1_L1_V1'


    # Data set up
    geoc_trigger_time = SampleParam[9][i]
    roll_off = 0.4  # Roll off duration of tukey window in seconds, default is 0.4s
    duration = 4  # Analysis segment duration
    #+++++++++++++++++++++
    """++++++++++++++++++"""
    maximum_frequency = 512 #这个频率怎么选比较好？？？？？？？？？？？？？？？？？？
    minimum_frequency = 20


    # We now use gwpy to obtain analysis and psd data and create the ifo_list
    ifo_list = bilby.gw.detector.InterferometerList([])

    fname_list = ['./Sim_GW_Data/H1/H-H1_GWOSC_' + sample_RATE + '_R1-' + Data_start_time + '-4096_' + str(i) + '.gwf' , './Sim_GW_Data/L1/L-L1_GWOSC_' + sample_RATE + '_R1-' + Data_start_time + '-4096_' + str(i) + '.gwf', './Sim_GW_Data/V1/V-V1_GWOSC_' + sample_RATE + '_R1-' + Data_start_time + '-4096_' + str(i) + '.gwf']
    fchannel_list = ['H1:GWOSC_' + sample_RATE + '_R1_STRAIN', 'L1:GWOSC_' + sample_RATE + '_R1_STRAIN', 'V1:GWOSC_' + sample_RATE + '_R1_STRAIN']
    psd_list = ['./Sim_PSD_Data/H1/test.gwf' , './Sim_PSD_Data/L1/test.gwf', './Sim_PSD_Data/V1/test.gwf']
    psd_channel = ['H1', 'L1', 'V1']

    gw_duration_20[i] = gw_duration_20[i] #不要延长1秒，为的是把信号完全包含在里面
    
    print(' duration 20 Hz = ' + str(gw_duration_20[i]))
    print(' duration 10 Hz = ' + str(gw_duration[i]))
    
    
    dt_from_earth_center_min = np.min([dt_h1_from_earth_center[i], dt_l1_from_earth_center[i], dt_v1_from_earth_center[i]])
    dt_from_earth_center_max = np.max([dt_h1_from_earth_center[i], dt_l1_from_earth_center[i], dt_v1_from_earth_center[i]])
        
    for det_index , det in enumerate(["H1", "L1", "V1"]):#(["H1", "L1", "V1", "K1"]):
        ifo = bilby.gw.detector.get_empty_interferometer(det)
        data = TimeSeries.read(fname_list[det_index], channel = fchannel_list[det_index])
        #这下面的data.pass应不应该加？在我过去的版本中是没有加入的
        # data.highpass(minimum_frequency)
        # data.lowpass(maximum_frequency)
        if det == "H1":
            # data_start_index = np.where(np.array(data.times) >= int(SampleParam[9][i]) + dt_from_earth_center_min)[0][0]
            # data_end_index = np.where(np.array(data.times) <= int(SampleParam[9][i]) + dt_from_earth_center_max + gw_duration_20[i])[0][-1]
            data_start_index = np.where(np.array(data.times) >= int(SampleParam[9][i]) - gw_duration_20[i])[0][0]
            data_end_index = np.where(np.array(data.times) <= int(SampleParam[9][i]) + dt_from_earth_center_max + gw_duration_20[i])[0][-1]
            print("H1 cutted data length: ", (data_end_index - data_start_index)/sample_rate)
            data_start = data.times[data_start_index]
            data_end = data.times[data_end_index]
            print("H1 data stat time: ", data_start)
            print("H1 data end time: ", data_end)
        elif det == "L1":
            # data_start_index = np.where(np.array(data.times) >= int(SampleParam[9][i]) + dt_from_earth_center_min)[0][0]
            # data_end_index = np.where(np.array(data.times) <= int(SampleParam[9][i]) + dt_from_earth_center_max + gw_duration_20[i])[0][-1]
            data_start_index = np.where(np.array(data.times) >= int(SampleParam[9][i]) - gw_duration_20[i])[0][0]
            data_end_index = np.where(np.array(data.times) <= int(SampleParam[9][i]) + dt_from_earth_center_max + gw_duration_20[i])[0][-1]
            print("L1 cutted data length: ", (data_end_index - data_start_index)/sample_rate)
            data_start = data.times[data_start_index]
            data_end = data.times[data_end_index]
            print("L1 data stat time: ", data_start)
            print("L1 data end time: ", data_end)
            
        elif det == "V1":
            # data_start_index = np.where(np.array(data.times) >= int(SampleParam[9][i]) + dt_from_earth_center_min)[0][0]
            # data_end_index = np.where(np.array(data.times) <= int(SampleParam[9][i]) + dt_from_earth_center_max + gw_duration_20[i])[0][-1]
            data_start_index = np.where(np.array(data.times) >= int(SampleParam[9][i]) - gw_duration_20[i])[0][0]
            data_end_index = np.where(np.array(data.times) <= int(SampleParam[9][i]) + dt_from_earth_center_max + gw_duration_20[i])[0][-1]
            print("V1 cutted data length: ", (data_end_index - data_start_index)/sample_rate)
            data_start = data.times[data_start_index]
            data_end = data.times[data_end_index]
            print("V1 data stat time: ", data_start)
            print("V1 data end time: ", data_end)
            
             
        data = data[data_start_index: data_end_index]
        ifo.strain_data.set_from_gwpy_timeseries(data)

        psd_data = TimeSeries.read(psd_list[det_index] , channel = psd_channel[det_index])
        
        psd_alpha = 2 * roll_off / duration
        psd = psd_data.psd(
            fftlength=duration,
            overlap=0,
            window=("tukey", psd_alpha),
            method="median"
        )
        ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
            frequency_array=psd.frequencies.value, psd_array=psd.value)
        ifo.maximum_frequency = maximum_frequency
        ifo.minimum_frequency = minimum_frequency

        ifo_list.append(ifo)

    logger.info("Saving data plots to {}".format(outdir))
    bilby.core.utils.check_directory_exists_and_if_not_mkdir(outdir)
    ifo_list.plot_data(outdir=outdir, label=label)

    # We now define the prior.
    # We have defined our prior distribution in a local file, GW150914.prior
    # The prior is printed to the terminal at run-time.
    # You can overwrite this using the syntax below in the file,
    # or choose a fixed value by just providing a float value as the prior.


    priors = bilby.gw.prior.BBHPriorDict()
    priors.pop('mass_1')
    priors.pop('mass_2')

    priors['chirp_mass'] = bilby.prior.Uniform(
    name='chirp_mass', latex_label='$M$', minimum=6, maximum=300.0,
    unit='$M_{\\odot}$')
    
    

    priors['luminosity_distance'] = bilby.prior.Uniform(
        minimum=400, maximum = 150000, name='luminosity_distance', latex_label='$d_L$', unit="Mpc"
        )
        
    
    # # Add the geocent time prior
    # priors["geocent_time"] = bilby.core.prior.Uniform(
    #     injection_parameters['geocent_time'] - 0.1, injection_parameters['geocent_time'] + 0.1, name="geocent_time"
    # )
    print("Merger time of H1 = ", trigger_time_H1[i])
    print("Merger time of L1 = ", trigger_time_L1[i])
    print("Merger time of V1 = ", trigger_time_V1[i])
    print("Merger time of Delta t H1L1= ", trigger_time_H1[i] - trigger_time_L1[i])
    print("Merger time of Delta t H1V1= ", trigger_time_H1[i] - trigger_time_V1[i])
    
    #新加的===============
    if H1_snr[i] >= L1_snr[i] and H1_snr[i] >= V1_snr[i]:
        priors["H1_time"] = bilby.core.prior.Uniform(
            minimum=trigger_time_H1[i] - 0.1,
            maximum=trigger_time_H1[i] + 0.1,
            name="H1_time",
            latex_label="$t_H$",
            unit="$s$",
        )
        time_reference_name = "H1"
        
    elif L1_snr[i] >= H1_snr[i] and L1_snr[i] >= V1_snr[i]:
        priors["L1_time"] = bilby.core.prior.Uniform(
            minimum=trigger_time_L1[i] - 0.1,
            maximum=trigger_time_L1[i] + 0.1,
            name="L1_time",
            latex_label="$t_L$",
            unit="$s$",
        )
        time_reference_name = "L1"
        
    elif V1_snr[i] >= H1_snr[i] and V1_snr[i] >= L1_snr[i]:
        priors["V1_time"] = bilby.core.prior.Uniform(
            minimum=trigger_time_V1[i] - 0.1,
            maximum=trigger_time_V1[i] + 0.1,
            name="V1_time",
            latex_label="$t_V$",
            unit="$s$",
        )
        time_reference_name = "V1"
        
    print("time reference name: = ", time_reference_name)
    
    del priors["ra"], priors["dec"]
    priors["zenith"] = bilby.core.prior.Sine(latex_label="$\\kappa$")
    priors["azimuth"] = bilby.core.prior.Uniform(
        minimum=0, maximum=2 * np.pi, latex_label="$\\epsilon$", boundary="periodic"
    )
    # priors["theta_jn"] = injection_parameters["theta_jn"]
    # priors["luminosity_distance"] = injection_parameters["luminosity_distance"]
    # priors["psi"] = injection_parameters["psi"]
    

    #===================

    # In this step we define a `waveform_generator`. This is the object which
    # creates the frequency-domain strain. In this instance, we are using the
    # `lal_binary_black_hole model` source model. We also pass other parameters:
    # the waveform approximant and reference frequency and a parameter conversion
    # which allows us to sample in chirp mass and ratio rather than component mass
    waveform_generator = bilby.gw.WaveformGenerator(
        # start_time = geoc_trigger_time,
        sampling_frequency=sample_rate, duration= gw_duration_20[i], #注意注意这个地方还有个duration
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments={'waveform_approximant': 'IMRPhenomPv2',#'catch_waveform_errors': True,
                            'reference_frequency': 50, 'minimum_frequency':minimum_frequency})





    # In this step, we define the likelihood. Here we use the standard likelihood
    # function, passing it the data and the waveform generator.
    # likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    #     ifo_list, waveform_generator, priors=priors, time_marginalization=True,
    #     phase_marginalization=True, distance_marginalization=True)
    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        ifo_list, waveform_generator, priors=priors, distance_marginalization=True, time_marginalization=True,
        phase_marginalization=True, jitter_time = True, reference_frame="H1L1V1", time_reference=time_reference_name,)






    # Finally, we run the sampler. This function takes the likelihood and prior
    # along with some options for how to do the sampling and how to save the data

    result = bilby.run_sampler(
        likelihood, priors, sample="acceptance-walk", sampler='dynesty', outdir=outdir, label=label,
        injection_parameters=injection_parameters, check_point_delta_t=3600 * 5, n_check_point=100000, check_point_plot=True,
        nlive=1500, naccept=60, maxmcmc = 7000,
        conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,result_class=bilby.gw.result.CBCResult,
        npool = 50)
    
        # Plot the inferred waveform superposed on the actual data.
    result.plot_waveform_posterior(n_samples=1000, start_time=-1.5, end_time=0.5)

    result.plot_corner()


#Bilby 15参数的运行结果
# result = bilby.run_sampler(
#     likelihood=likelihood,
#     priors=priors,
#     sampler="dynesty",
#     nlive=1000,
#     walks=20,
#     nact=50,
#     maxmcmc=2000,
#     npool=4,
#     injection_parameters=injection_parameters,
#     outdir=outdir,
#     label=label,
#     conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,
#     result_class=bilby.gw.result.CBCResult,
# )
