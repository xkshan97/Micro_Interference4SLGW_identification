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
import random
from astropy.cosmology import Planck18
cosmo = Planck18
from lenstronomy.Workflow.fitting_sequence import FittingSequence


JWST_file = 'JWST_host_pop_power_law'

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
angular_diameter_distance_to_lens = Planck18.angular_diameter_distance(z_lens).to(u.kpc).value
z_source = SampleParam[0][int(AcceptLensIndex[index])]
sigma_v = Sigma_v_set[index]
q = Q_set[index]
y1 = y1_set[index]
y2 = y2_set[index]

#source pop
JWST_mag_Res = np.loadtxt('./SampleResult_GWGalaxy/JWST_mag_Res_' + str(index) + '.csv', delimiter=',')
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
# for i in range(1):
def Main_func(index_tmp):
    # index_tmp = 25
    index_tmp = int(index_tmp)
    #mag_min_Res = [mag_min, axis_ratio_min, position_angle_min, R_e_min, sersic_n_min, np.sum(np.array(mag_weight_by_SFR_50)<25.5) / len(mag_weight_by_SFR_50), np.sum(np.array(mag_weight_by_SFR_100)<25.5) / len(mag_weight_by_SFR_100)]
    from lenstronomy.SimulationAPI.ObservationConfig.JWST import JWST
    JWST_F200W = JWST(band='F200W', psf_type='PIXEL', coadd_years=None)
    kwargs_F200W_band = JWST_F200W.kwargs_single_band()
    kwargs_F200W_band['seeing'] = 0.064 #https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions#NIRCamPointSpreadFunctions-PSFFWHM
    kwargs_F200W_band['psf_type'] = 'GAUSSIAN'
    # kwargs_F200W_band['exposure_time'] = 3600 
    # kwargs_F200W_band['num_exposures'] = 24

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

    numpix = 150  # number of pixels per axis of the image to be modelled

    # here we define the numerical options used in the ImSim module. 
    # Have a look at the ImageNumerics class for detailed descriptions.
    # If not further specified, the default settings are used.
    kwargs_numerics = {'supersampling_factor': 1, 'supersampling_convolution': False}

    sim_F200W = SimAPI(numpix=numpix, kwargs_single_band=kwargs_F200W_band, kwargs_model=kwargs_model)

    # return the ImSim instance. With this class instance, you can compute all the
    # modelling accessible of the core modules. See class documentation and other notebooks.
    imSim_F200W = sim_F200W.image_model_class(kwargs_numerics)




    '''
    Brightness definitions in magnitude space

    One core feature is the support of light profile amplitudes in astronomical magnitude space (at least for few selected well defined brightness profiles).

    We first define all parameters in magnitude space and then use the SimAPI routine to translate the arguments into lenstronomy conventions used by the ImSim module. The second model of each light component we defined as 'INTERPOL', which sets an interpolation grid given an image. This can be used to past real galaxies as lenses or sources into lenstronomy.
    '''

    # F200W-band

    # lens light
    kwargs_lens_light_mag_F200W = [{'magnitude': GW_Lens_Galaxy_Res[0], 'R_sersic': GW_Lens_Galaxy_Res[1] / angular_diameter_distance_to_lens * u.rad.to(u.arcsec), 'n_sersic': 4, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}]
    # source light
    e1_s, e2_s = param_util.phi_q2_ellipticity(phi=JWST_mag_Res[index_tmp][2] * np.pi / 180, q=JWST_mag_Res[index_tmp][1])
    kwargs_source_mag_F200W = [{'magnitude': JWST_mag_Res[index_tmp][0], 'R_sersic': JWST_mag_Res[index_tmp][3], 'n_sersic': JWST_mag_Res[index_tmp][4], 'e1': e1_s, 'e2': e2_s, 'center_x': source_x, 'center_y': source_y}]

    # turn magnitude kwargs into lenstronomy kwargs
    kwargs_lens_light_F200W, kwargs_source_F200W, kwargs_ps_F200W = sim_F200W.magnitude2amplitude(kwargs_lens_light_mag_F200W, kwargs_source_mag_F200W)


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

    image_F200W = imSim_F200W.image(kwargs_lens, kwargs_source_F200W, kwargs_lens_light_F200W)

    # add noise
    image_F200W += sim_F200W.noise_for_model(model=image_F200W)


    # and plot it

    img = np.zeros((image_F200W.shape[0], image_F200W.shape[1], 1), dtype=float)
    img[:,:,0] = plot_util.sqrt(image_F200W)

    plot_size = numpix * kwargs_F200W_band['pixel_scale']
    f, ax = plt.subplots(1, 1, figsize=(4, 4))

    ax.imshow(img, aspect='equal', origin='lower', extent=[0, plot_size, 0, plot_size])
    ax.axis('off')
    scale_bar(ax, d=plot_size, dist=1., text='1"', color='w', font_size=15, flipped=False)
    plt.title('source_mag = ' + str(round(JWST_mag_Res[index_tmp][0], 2)) + 'SFR = ' + str(round(JWST_mag_Res[index_tmp][5], 2)))
    plt.savefig('./' + JWST_file +'/' + str(index) + '_' + str(index_tmp) + '_Sim_fig.png', dpi=450)
    plt.close()

    kwargs_data_F200W = {'background_rms': sim_F200W.kwargs_data['background_rms'], 'exposure_time': sim_F200W.kwargs_data['exposure_time'], 'ra_at_xy_0': sim_F200W.kwargs_data['ra_at_xy_0'], 'dec_at_xy_0': sim_F200W.kwargs_data['dec_at_xy_0'], 'transform_pix2angle': sim_F200W.kwargs_data['transform_pix2angle'], 'image_data': image_F200W}

    image_band = [kwargs_data_F200W, sim_F200W.kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]


    # lens models
    fixed_lens = []
    kwargs_lens_init = []
    kwargs_lens_sigma = []
    kwargs_lower_lens = []
    kwargs_upper_lens = []

    fixed_lens = [{}]  # for this example, we fix the power-law index of the lens model to be isothermal
    kwargs_lens_init = [{'theta_E': theta_E, 'e1': e1, 'e2': e2, 'gamma': 2,
                            'center_x': 0., 'center_y': 0.}]
    kwargs_lens_sigma = [{'theta_E': 0.2, 'e1': 0.05, 'e2': 0.05, 'gamma': 0.05,
                            'center_x': 0.05, 'center_y': 0.05}]
    kwargs_lower_lens = [{'theta_E': 0.01, 'e1': -0.5, 'e2': -0.5, 'gamma': 1.5, 'center_x': -1, 'center_y': -1}]
    kwargs_upper_lens = [{'theta_E': 2., 'e1': 0.5, 'e2': 0.5, 'gamma': 2.5, 'center_x': 1, 'center_y': 1}]

    fixed_lens.append({'ra_0': 0, 'dec_0': 0})

    lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]


    fixed_source = []
    kwargs_source_init = []
    kwargs_source_sigma = []
    kwargs_lower_source = []
    kwargs_upper_source = []

    fixed_source.append({})
    kwargs_source_init.append({'R_sersic': JWST_mag_Res[index_tmp][3], 'n_sersic': JWST_mag_Res[index_tmp][4], 'e1': e1_s, 'e2': e2_s, 'center_x': source_x, 'center_y': source_y})
    kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.1, 'e1': 0.05, 'e2': 0.05, 'center_x': 0.2, 'center_y': 0.2})
    kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.1, 'center_x': -1, 'center_y': -1})
    kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 1, 'n_sersic': 7., 'center_x': 1, 'center_y': 1})
    
    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]


    fixed_lens_light = []
    kwargs_lens_light_init = []
    kwargs_lens_light_sigma = []
    kwargs_lower_lens_light = []
    kwargs_upper_lens_light = []


    fixed_lens_light.append({})
    kwargs_lens_light_init.append({'R_sersic': GW_Lens_Galaxy_Res[1] / angular_diameter_distance_to_lens * u.rad.to(u.arcsec), 'n_sersic': 3, 'e1': e1, 'e2': e2, 'center_x': 0., 'center_y': 0})
    kwargs_lens_light_sigma.append({'n_sersic': 1, 'R_sersic': 0.3, 'e1': 0.05, 'e2': 0.05, 'center_x': 0.1, 'center_y': 0.1})
    kwargs_lower_lens_light.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': .5, 'center_x': -1, 'center_y': -1})
    kwargs_upper_lens_light.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 5, 'n_sersic': 5., 'center_x': 1, 'center_y': 1})
    
    lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]

    kwargs_params = {'lens_model': lens_params,
                    'source_model': source_params,
                    'lens_light_model': lens_light_params}


    # we have the option to set priors in the model parameters.
    # prior_lens = [[0, 'e1', 0, 0.2], [0, 'e2', 0, 0.2]]#  [index_model, param_name, mean, 1-sigma priors]
    # have a look in the LikelihoodModule for a complete description of implemented priors


    kwargs_likelihood = {'check_bounds': True,
                     'force_no_add_image': False,
                     'source_marg': False,
                     'image_position_uncertainty': 0.004,
                     'check_matched_source_position': True,
                     'source_position_tolerance': 0.001,
                     'source_position_sigma': 0.001
                    #  'prior_lens': prior_lens
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
    kwargs_model_fit = {'lens_model_list': ['EPL'], #['SIE'],  # list of lens models to be used
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

    modelPlot = ModelPlot(multi_band_list, kwargs_model_fit, kwargs_result, arrow_size=0.02, cmap_string="gist_heat")
        
    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.data_plot(ax=axes[0,0])
    modelPlot.model_plot(ax=axes[0,1])
    modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
    modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=400, v_min=-4, v_max=-1)
    modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
    modelPlot.magnification_plot(ax=axes[1, 2])
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'rec_obs.png', dpi=450)
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
    plt.savefig('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'decomposition.png', dpi=450)
    plt.close()
    print(kwargs_result)


    from lenstronomy.Plots import chain_plot

    for i in range(len(chain_list)):
        chain_plot.plot_chain_list(chain_list, i)
        plt.savefig('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'chain.png', dpi=450)
        
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[1]

    param_class = fitting_seq.param_class
    kwargs_lens[0]['gamma'] = 2
    param_truths = param_class.kwargs2args(kwargs_lens=kwargs_lens, kwargs_source=kwargs_source_F200W, kwargs_lens_light=kwargs_lens_light_F200W)

    print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
    print("parameters in order: ", param_mcmc)
    print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])

    if not samples_mcmc == []:
        n, num_param = np.shape(samples_mcmc)
        plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True, truths=param_truths)
        plt.savefig('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'corner_lens_full.png', dpi=450)
        plt.close() 
        plot = corner.corner(samples_mcmc[:,8:], labels=param_mcmc[8:], show_titles=True, truths=param_truths[8:])
        plt.savefig('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'corner_lens_sub.png', dpi=450)
        plt.close() 
    
    np.savetxt('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'samples_mcmc.csv', samples_mcmc, delimiter=',')
    np.save('./' + JWST_file + '/' + str(index) + '_' + str(index_tmp) + 'kwargs_result.npy', kwargs_result)            


# Generate_file(0,1)
if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = len(JWST_mag_Res)
    length = threadcount#int(threadcount * 150) #
    
    
    percount = length//threadcount

    print("percount = ", percount)

    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
    
        res = pool.apply_async(func=Main_func, args=(i,))
        
    pool.close()
    pool.join()