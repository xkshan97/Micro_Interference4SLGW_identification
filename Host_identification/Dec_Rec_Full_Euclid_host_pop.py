'''
This notebooks provides examples in how to use the lenstronomy.SimulationAPI modules in simulating (realistic) mock lenses taylored to a specific observation and instrument.

The module enables to use the astronomical magnitude conventions and can translate those into the lenstronomy core module configurations.
'''
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import copy
import numpy as np
import corner
from PIL import Image
import matplotlib.pyplot as plt
import lenstronomy.Util.data_util as data_util
import lenstronomy.Util.util as util
import lenstronomy.Plots.plot_util as plot_util
from lenstronomy.Plots.plot_util import coordinate_arrows, scale_bar
from lenstronomy.SimulationAPI.sim_api import SimAPI
from lenstronomy.LightModel.Profiles.gaussian import GaussianEllipse
gauss = GaussianEllipse()
from astropy.constants import c
import astropy.units as u
import lenstronomy.Util.param_util as param_util
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from  multiprocessing import Process,Pool

from astropy.cosmology import Planck18
cosmo = Planck18
from lenstronomy.Workflow.fitting_sequence import FittingSequence


Euclid_file = 'Euclid_host_pop_long_time'

SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')

magnification = np.loadtxt('../Paper4_CE_Modify/SampleResult/magnification.csv', delimiter=',')
timedelay = np.loadtxt('../Paper4_CE_Modify/SampleResult/timedelayDays.csv', delimiter=',')
kappa = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa.csv', delimiter=',')
gamma = np.loadtxt('../Paper4_CE_Modify/SampleResult/gamma.csv', delimiter=',')
imagenum = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
lens_z = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')
Sigma_v_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/sigma_v.csv', delimiter=',')
Q_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/axis_ratio.csv', delimiter=',')
y1_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/y1.csv', delimiter=',')
y2_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/y2.csv', delimiter=',')
x_image_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/x1_image.csv', delimiter=',')
y_image_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/x2_image.csv', delimiter=',')

index = 120
z_lens = lens_z[index]
z_source = SampleParam[0][int(AcceptLensIndex[index])]
sigma_v = Sigma_v_set[index]
q = Q_set[index]
y1 = y1_set[index]
y2 = y2_set[index]

#source pop
mag_Res = np.loadtxt('./SampleResult_GWGalaxy/mag_Res_' + str(index) + '.csv', delimiter=',')
GW_Lens_Galaxy_Res = np.loadtxt('./SampleResult_GWGalaxy/GW_Lens_Galaxy_' + str(index) + '.csv', delimiter=',')
# lensing quantities
ckm = c.value * (u.m/u.s).to(u.km/u.s) #千米每秒单位下的光速
D_ds = Planck18.angular_diameter_distance_z1z2(z_lens, z_source).value
D_s = Planck18.angular_diameter_distance(z_source).value
theta_E = 4 * np.pi * sigma_v ** 2 / ckm ** 2 * D_ds / D_s * 3600 * 180 / np.pi #以角秒为单位
e1, e2 = param_util.phi_q2_ellipticity(phi=np.pi/2, q=q)

source_x = y1 * theta_E
source_y = y2 * theta_E

'''
Define camera and observations

As an example, we define the camera and observational settings of a LSST-like observation. We define one camera setting and three different observations corresponding th g,r,i imaging.

For the complete list of possible settings, we refer to the SimulationAPI.observation_api classes. There are pre-configured settings which approximately mimic observations from current and future instruments. Be careful using those and check whether they are sufficiently accurate for your specific science case!
'''
# Instrument setting from pre-defined configurations

for index_tmp in range(len(mag_Res)): #宿主星系population
    if mag_Res[index_tmp][0] < 25.5: 
        #mag_min_Res = [mag_min, axis_ratio_min, position_angle_min, R_e_min, sersic_n_min, np.sum(np.array(mag_weight_by_SFR_50)<25.5) / len(mag_weight_by_SFR_50), np.sum(np.array(mag_weight_by_SFR_100)<25.5) / len(mag_weight_by_SFR_100)]
        from lenstronomy.SimulationAPI.ObservationConfig.Euclid import Euclid
        Euclid_VIS = Euclid(band='VIS', psf_type='GAUSSIAN', coadd_years=6)
        kwargs_VIS_band = Euclid_VIS.kwargs_single_band()
        # kwargs_VIS_band['exposure_time'] = 3600 
        kwargs_VIS_band['num_exposures'] = 40

        '''
        Define model settings

        The model settings are handled by the SimulationAPI.model_api ModelAPI class. The role is to return instances of the lenstronomy LightModel, LensModel, PointSource modules according to the options chosen by the user. Currently, all other model choices are equivalent to the ones provided by LightModel, LensModel, PointSource. The current options of the class instance only describe a subset of possibilities and we refer to the specific class instances for details about all the possibilities.

        For this example, we chose a single lens plane and a single source plane, elliptical Sersic profiles and an additional lensed point source.
        '''

        kwargs_model = {'lens_model_list': ['SIE'],  # list of lens models to be used
                        'lens_light_model_list': ['SERSIC_ELLIPSE'],  # list of unlensed light models to be used
                        'source_light_model_list': ['SERSIC_ELLIPSE'],  # list of extended source models to be used
            }

        '''
        Generate SimAPI instance for the different observations

        Here we make an instance of the SimAPI class and execute the major tasks accessible as the interface to the ImSim core module.
        '''

        numpix = 60  # number of pixels per axis of the image to be modelled

        # here we define the numerical options used in the ImSim module. 
        # Have a look at the ImageNumerics class for detailed descriptions.
        # If not further specified, the default settings are used.
        kwargs_numerics = {'supersampling_factor': 1, 'supersampling_convolution': False}

        sim_VIS = SimAPI(numpix=numpix, kwargs_single_band=kwargs_VIS_band, kwargs_model=kwargs_model)

        # return the ImSim instance. With this class instance, you can compute all the
        # modelling accessible of the core modules. See class documentation and other notebooks.
        imSim_VIS = sim_VIS.image_model_class(kwargs_numerics)




        '''
        Brightness definitions in magnitude space

        One core feature is the support of light profile amplitudes in astronomical magnitude space (at least for few selected well defined brightness profiles).

        We first define all parameters in magnitude space and then use the SimAPI routine to translate the arguments into lenstronomy conventions used by the ImSim module. The second model of each light component we defined as 'INTERPOL', which sets an interpolation grid given an image. This can be used to past real galaxies as lenses or sources into lenstronomy.
        '''

        # VIS-band

        # lens light
        kwargs_lens_light_mag_VIS = [{'magnitude': GW_Lens_Galaxy_Res[0], 'R_sersic': GW_Lens_Galaxy_Res[1], 'n_sersic': 4, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}]
        # source light
        e1_s, e2_s = param_util.phi_q2_ellipticity(phi=mag_Res[index_tmp][2] * np.pi / 180, q=mag_Res[index_tmp][1])
        kwargs_source_mag_VIS = [{'magnitude': mag_Res[index_tmp][0], 'R_sersic': mag_Res[index_tmp][3], 'n_sersic': mag_Res[index_tmp][4], 'e1': e1_s, 'e2': e2_s, 'center_x': source_x, 'center_y': source_y}]

        # turn magnitude kwargs into lenstronomy kwargs
        kwargs_lens_light_VIS, kwargs_source_VIS, kwargs_ps_VIS = sim_VIS.magnitude2amplitude(kwargs_lens_light_mag_VIS, kwargs_source_mag_VIS)


        '''
        setting lens model parameters

        And finally we need a lens model. The default lensing units are in observed deflection angle (arc seconds) mapping the image to the source plane (reduced angles). In a single lens and single source plane model, this is all we need to specify and no futher cosmology is required.
        '''

        kwargs_lens = [
            {'theta_E': theta_E, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}  # SIE model
        ]

        '''
        simluate image

        Finally we can simulate the images with the ImageModel class instance and the lenstronomy parameters for the different bands. Note that in the specific example we included a point source (representing a quasar) in the center of the lensed galaxy. The SimulationAPI supports various options. Further down below we simulate multi-lens multi-source plane configurations too.
        '''

        image_VIS = imSim_VIS.image(kwargs_lens, kwargs_source_VIS, kwargs_lens_light_VIS)

        # add noise
        image_VIS += sim_VIS.noise_for_model(model=image_VIS)


        # and plot it

        img = np.zeros((image_VIS.shape[0], image_VIS.shape[1], 1), dtype=float)
        img[:,:,0] = plot_util.sqrt(image_VIS)

        plot_size = numpix * kwargs_VIS_band['pixel_scale']
        f, ax = plt.subplots(1, 1, figsize=(4, 4))

        ax.imshow(img, aspect='equal', origin='lower', extent=[0, plot_size, 0, plot_size])
        ax.axis('off')
        scale_bar(ax, d=plot_size, dist=1., text='1"', color='w', font_size=15, flipped=False)
        plt.title('source_mag = ' + str(round(mag_Res[index_tmp][0], 2)) + 'SFR = ' + str(round(mag_Res[index_tmp][5], 2)))
        plt.savefig('./' + Euclid_file +'/' + str(index) + '_' + str(index_tmp) + '_Sim_fig.png', dpi=450)
        plt.close()

        kwargs_data_VIS = {'background_rms': sim_VIS.kwargs_data['background_rms'], 'exposure_time': sim_VIS.kwargs_data['exposure_time'], 'ra_at_xy_0': sim_VIS.kwargs_data['ra_at_xy_0'], 'dec_at_xy_0': sim_VIS.kwargs_data['dec_at_xy_0'], 'transform_pix2angle': sim_VIS.kwargs_data['transform_pix2angle'], 'image_data': image_VIS}

        image_band = [kwargs_data_VIS, sim_VIS.kwargs_psf, kwargs_numerics]
        multi_band_list = [image_band]


        # lens models
        fixed_lens = []
        kwargs_lens_init = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []

        fixed_lens.append({})  # for this example, we fix the power-law index of the lens model to be isothermal
        kwargs_lens_init.append({'theta_E': 0.7, 'e1': 0., 'e2': 0.,
                                'center_x': 0., 'center_y': 0.})
        kwargs_lens_sigma.append({'theta_E': .2, 'e1': 0.05, 'e2': 0.05,
                                'center_x': 0.05, 'center_y': 0.05})
        kwargs_lower_lens.append({'theta_E': 0.01, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_lens.append({'theta_E': 10., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10})

        fixed_lens.append({'ra_0': 0, 'dec_0': 0})

        lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]


        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []

        fixed_source.append({})
        kwargs_source_init.append({'R_sersic': 0.2, 'n_sersic': 1, 'e1': 0.3, 'e2': 0, 'center_x': 0., 'center_y': 0})
        kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.1, 'e1': 0.05, 'e2': 0.05, 'center_x': 0.2, 'center_y': 0.2})
        kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})

        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]


        fixed_lens_light = []
        kwargs_lens_light_init = []
        kwargs_lens_light_sigma = []
        kwargs_lower_lens_light = []
        kwargs_upper_lens_light = []


        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'R_sersic': 0.5, 'n_sersic': 2, 'e1': 0, 'e2': 0, 'center_x': 0., 'center_y': 0})
        kwargs_lens_light_sigma.append({'n_sersic': 1, 'R_sersic': 0.3, 'e1': 0.05, 'e2': 0.05, 'center_x': 0.1, 'center_y': 0.1})
        kwargs_lower_lens_light.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001, 'n_sersic': .5, 'center_x': -10, 'center_y': -10})
        kwargs_upper_lens_light.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10})

        lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]

        kwargs_params = {'lens_model': lens_params,
                        'source_model': source_params,
                        'lens_light_model': lens_light_params}



        kwargs_likelihood = {'source_marg': False,
                            }

        kwargs_data_joint = {'multi_band_list': multi_band_list, 
                            'multi_band_type': 'multi-linear'  # 'multi-linear': every imaging band has independent solutions of the surface brightness, 'joint-linear': there is one joint solution of the linear coefficients demanded across the bands.
                            }


        #    'joint_source_with_source':list [[i_source, k_source, ['param_name1', 'param_name2', ...]], [...], ...],
        #    joint parameter between two source surface brightness models, the second adopts the value of the first

        # here we add (optional) constraints between the light profile of the two different imaging bands. We demand the center ellipticity to be the same.
        # effectively this leaves the half-light radius and Sersic index free between the colors.
        kwargs_constraints = {'linear_solver': True}  # optional, if 'linear_solver': False, lenstronomy does not apply a linear inversion of the 'amp' parameters during fitting but instead samples them.


        index_source_light_model_list = [[0], [1]]
        kwargs_model_fit = {'lens_model_list': ['SIE'],  # list of lens models to be used
                        'lens_light_model_list': ['SERSIC_ELLIPSE'],  # list of unlensed light models to be used
                        'source_light_model_list': ['SERSIC_ELLIPSE'],  # list of extended source models to be used
            }


        fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model_fit, kwargs_constraints, kwargs_likelihood, kwargs_params, mpi=False)

        fitting_kwargs_list = [['PSO', {'sigma_scale': 1., 'n_particles': 200, 'n_iterations': 200}],
                            ['MCMC', {'n_burn': 200, 'n_run': 600, 'n_walkers': 200, 'sigma_scale': .1}]
                ]

        chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
        kwargs_result = fitting_seq.best_fit()



        #Plot
        from lenstronomy.Plots.model_plot import ModelPlot

        modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat")
            
        f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

        modelPlot.data_plot(ax=axes[0,0])
        modelPlot.model_plot(ax=axes[0,1])
        modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=100, v_min=-5, v_max=-3)
        modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
        modelPlot.magnification_plot(ax=axes[1, 2])
        f.tight_layout()
        f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        plt.savefig('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'rec_obs.png', dpi=450)
        plt.close()

        f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

        modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True)
        modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True)
        modelPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True)
        modelPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
        modelPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True, unconvolved=True)
        modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        f.tight_layout()
        f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        plt.savefig('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'decomposition.png', dpi=450)
        plt.close()
        print(kwargs_result)


        from lenstronomy.Plots import chain_plot

        for i in range(len(chain_list)):
            chain_plot.plot_chain_list(chain_list, i)
            plt.savefig('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'chain.png', dpi=450)
            
        sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]

        param_class = fitting_seq.param_class
        param_truths = param_class.kwargs2args(kwargs_lens=kwargs_lens, kwargs_source=kwargs_source_VIS, kwargs_lens_light=kwargs_lens_light_VIS)

        print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
        print("parameters in order: ", param_mcmc)
        print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])

        if not samples_mcmc == []:
                n, num_param = np.shape(samples_mcmc)
                plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True, truths=param_truths)
                plt.savefig('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'corner_lens_full.png', dpi=450)
                plt.close() 
                plot = corner.corner(samples_mcmc[:,8:], labels=param_mcmc[8:], show_titles=True, truths=param_truths[8:])
                plt.savefig('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'corner_lens_sub.png', dpi=450)
                plt.close() 
                

        # here we plot a subset of the 'processed' MCMC posteriors

        # import the parameter handling class #
        from lenstronomy.Sampling.parameters import Param
        from lenstronomy.Analysis.td_cosmography import TDCosmography
        from tqdm import tqdm
        point_source_list = ['LENSED_POSITION']
        kwargs_model_td = {'lens_model_list': ['SIE'],
                                    'source_light_model_list': ['SERSIC_ELLIPSE'],
                                    'lens_light_model_list': ['SERSIC_ELLIPSE'],
                                    'point_source_model_list': point_source_list,
                                    'fixed_magnification_list': [False],# list of bools (same length as point_source_type_list). If True, magnification ratio of point sources is fixed to the one given by the lens model 
                                    }
        td_cosmo = TDCosmography(z_lens, z_source, kwargs_model_td, cosmo_fiducial=cosmo)
        lens_model_class = LensModel(lens_model_list=['SIE'], z_lens=z_lens, z_source=z_source, cosmo=cosmo)

        lensEquationSolver = LensEquationSolver(lens_model_class)
        # make instance of parameter class with given model options, constraints and fixed parameters #

        param = Param(kwargs_model, fixed_lens, fixed_source, fixed_lens_light, 
                    kwargs_lens_init=kwargs_result['kwargs_lens'], **kwargs_constraints)


        def Generate_mcmc_td(IndexStart, IndexEnd):
            theta_E_result = []
            center_x_result = []
            center_y_result = []
            e1_result = []
            e2_result = []
            t_days_0_result = []
            t_days_1_result = []
            t_days_2_result = []
            t_days_3_result = []
            for i in tqdm(range(IndexStart, IndexEnd)):
                # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
                
                kwargs_result_tmp = param.args2kwargs(samples_mcmc[i])
                theta_E_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['theta_E']
                center_x_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['center_x']
                center_y_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['center_y']
                e1_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['e1']
                e2_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['e2']
                
                x_source_result_tmp = kwargs_result_tmp['kwargs_source'][0]['center_x'] 
                y_source_result_tmp = kwargs_result_tmp['kwargs_source'][0]['center_y'] 
                
                kwargs_sie_result_tmp = {'theta_E': theta_E_result_tmp, 'center_x': center_x_result_tmp, 'center_y': center_y_result_tmp, 'e1': e1_result_tmp, 'e2': e2_result_tmp}  # parameters of the deflector lens model
                kwargs_lens_result_tmp = [kwargs_sie_result_tmp]
                
                x_image_result_tmp, y_image_result_tmp = lensEquationSolver.findBrightImage(x_source_result_tmp, y_source_result_tmp, kwargs_lens_result_tmp, numImages=4,
                                                                min_distance=0.5*kwargs_VIS_band['pixel_scale'], search_window=numpix * kwargs_VIS_band['pixel_scale'])
                ps_amp_result_tmp = [1] * len(x_image_result_tmp) 
                if len(x_image_result_tmp) == 4:
                    kwargs_ps_result_tmp = [{'ra_image': x_image_result_tmp, 'dec_image': y_image_result_tmp,
                                'point_amp': ps_amp_result_tmp}]  # quasar point source position in the source plane and intrinsic brightness
                    
                    t_days = td_cosmo.time_delays(kwargs_lens_result_tmp, kwargs_ps_result_tmp, kappa_ext=0)
                    
                    theta_E_result.append(theta_E_result_tmp)
                    center_x_result.append(center_x_result_tmp)
                    center_y_result.append(center_y_result_tmp)
                    e1_result.append(e1_result_tmp)
                    e2_result.append(e2_result_tmp)
                    t_days_0_result.append(t_days[0])
                    t_days_1_result.append(t_days[1])
                    t_days_2_result.append(t_days[2])
                    t_days_3_result.append(t_days[3])
                
            
            return theta_E_result, center_x_result, center_y_result, e1_result, e2_result, t_days_0_result, t_days_1_result, t_days_2_result, t_days_3_result 


        # Generate_file(0,1)
        if __name__=='__main__':
            # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
            threadcount = 200
            res_z = [0 for i in range(threadcount)]
            length = len(samples_mcmc)#int(threadcount * 500) #
            
            
            percount = length//threadcount

            print("percount = ", percount)
        
            # t0 = time.time()
            pool = Pool(threadcount) 
            # for i in range(len(m1_set)):
            for i in range(threadcount):
                if i < threadcount - 1:
                    res = pool.apply_async(func=Generate_mcmc_td, args=(percount * i, percount * (i + 1),))
                    res_z[i] = res
                    # print("0",i)
                else:
                    # print("1",i)
                    res = pool.apply_async(func=Generate_mcmc_td, args=(percount * i, length,)) 
                    res_z[i] = res

            pool.close()
            pool.join() 
            
        Result1, Result2, Result3, Result4, Result5, Result6, Result7, Result8, Result9  = list(res_z[0].get()[0]), list(res_z[0].get()[1]), list(res_z[0].get()[2]), list(res_z[0].get()[3]), list(res_z[0].get()[4]), list(res_z[0].get()[5]), list(res_z[0].get()[6]), list(res_z[0].get()[7]), list(res_z[0].get()[8])
        for i in range(threadcount-1):
            Result1.extend(list(res_z[i+1].get()[0]))
            Result2.extend(list(res_z[i+1].get()[1]))
            Result3.extend(list(res_z[i+1].get()[2]))
            Result4.extend(list(res_z[i+1].get()[3]))
            Result5.extend(list(res_z[i+1].get()[4]))
            Result6.extend(list(res_z[i+1].get()[5]))
            Result7.extend(list(res_z[i+1].get()[6]))
            Result8.extend(list(res_z[i+1].get()[7]))
            Result9.extend(list(res_z[i+1].get()[8]))
                

        mcmc_new_list = []
        for i in range(len(Result1)):
            mcmc_new_list.append([Result1[i], Result2[i], Result3[i], Result4[i], Result5[i], Result6[i], Result7[i], Result8[i], Result9[i]])
                
        x_image_truths, y_image_truths = lensEquationSolver.findBrightImage(source_x, source_y, kwargs_lens, numImages=4,
                                                            min_distance=0.1*kwargs_VIS_band['pixel_scale'], search_window=10 * numpix * kwargs_VIS_band['pixel_scale'])
        kwargs_ps_truths = [{'ra_image': x_image_truths, 'dec_image': y_image_truths,
                        'point_amp': [1] * len(x_image_truths)}]  # quasar point source position in the source plane and intrinsic brightness

        t_days_truths = td_cosmo.time_delays(kwargs_lens, kwargs_ps_truths, kappa_ext=0)

        labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3"]
        try:
            plot = corner.corner(np.array(mcmc_new_list), labels=labels_new, show_titles=True, truths=[kwargs_lens[0]['theta_E'], kwargs_lens[0]['center_x'], kwargs_lens[0]['center_y'], kwargs_lens[0]['e1'], kwargs_lens[0]['e2'], t_days_truths[0], t_days_truths[1], t_days_truths[2], t_days_truths[3]])
            plt.savefig('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'corner_Res.png', dpi=450)
            plt.close()
            np.savetxt('./' + Euclid_file + '/' + str(index) + '_' + str(index_tmp) + 'corner_Res_data.csv', mcmc_new_list, delimiter=',')
        except AssertionError:
            pass