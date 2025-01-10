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
from astropy.modeling.models import Sersic2D

JWST_file = 'JWST_unhost_26_little_sigma'

Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images = np.loadtxt('SampleResult_Galaxy/Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images.csv', delimiter=',')
galaxy_redshift = np.loadtxt('./SampleResult_Galaxy/GalaxySourceLensedRedshift.csv', delimiter=',')
JWST_apparent_mag_set = np.loadtxt('./SampleResult_Galaxy/JWST_ApparentMag.csv', delimiter=',')
source_axis_ratio = np.loadtxt('./SampleResult_Galaxy/GalaxySourceAxis_ratio.csv', delimiter=',')
source_position_angle = np.loadtxt('./SampleResult_Galaxy/GalaxySourcePosition_angle.csv', delimiter=',')
source_Re_circ = np.loadtxt('./SampleResult_Galaxy/GalaxySourceRe_circ.csv', delimiter=',')
source_Re_circ_arc_sec = np.loadtxt('./SampleResult_Galaxy/R_e_arc_sec.csv', delimiter=',')
source_Sersic_n = np.loadtxt('./SampleResult_Galaxy/GalaxySourceSersic_n.csv', delimiter=',')



timedelay = np.loadtxt('./SampleResult_Galaxy/timedelayDays.csv', delimiter=',')
imagenum = np.loadtxt('./SampleResult_Galaxy/imagenumber.csv', delimiter=',')
lens_z = np.loadtxt('./SampleResult_Galaxy/lens_z.csv', delimiter=',')
Sigma_v_set = np.loadtxt('./SampleResult_Galaxy/sigma_v.csv', delimiter=',')
Q_set = np.loadtxt('./SampleResult_Galaxy/axis_ratio.csv', delimiter=',')
Lens_light_mag = np.loadtxt('./SampleResult_Galaxy/Lens_light_mag.csv', delimiter=',')
Lens_R_e = np.loadtxt('./SampleResult_Galaxy/Lens_R_e.csv', delimiter=',')
Lens_R_e_arc_sec = np.loadtxt('./SampleResult_Galaxy/R_e_Lens_arc_sec.csv', delimiter=',')
theta_E_set = np.loadtxt('./SampleResult_Galaxy/theta_E.csv', delimiter=',')

y1_set = np.loadtxt('./SampleResult_Galaxy/y1.csv', delimiter=',')
y2_set = np.loadtxt('./SampleResult_Galaxy/y2.csv', delimiter=',')
x_image_set = np.loadtxt('./SampleResult_Galaxy/x1_image.csv', delimiter=',')
y_image_set = np.loadtxt('./SampleResult_Galaxy/x2_image.csv', delimiter=',')
from lenstronomy.SimulationAPI.ObservationConfig.JWST import JWST
JWST_F200W = JWST(band='F200W', psf_type='PIXEL', coadd_years=None)
kwargs_F200W_band = JWST_F200W.kwargs_single_band()
kwargs_F200W_band['seeing'] = 0.064 #https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions#NIRCamPointSpreadFunctions-PSFFWHM
kwargs_F200W_band['psf_type'] = 'GAUSSIAN'
# kwargs_F200W_band['exposure_time'] = 3600 
# kwargs_F200W_band['num_exposures'] = 24
numpix = 150


for index in Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images:
    index = int(index)
    samples_mcmc = np.loadtxt('./' + JWST_file + '/' + str(index) + 'samples_mcmc.csv', delimiter=',')
    kwargs_result = np.load('./' + JWST_file + '/' + str(index) + 'kwargs_result.npy', allow_pickle=True).item()

    index = int(index)
    z_lens = lens_z[index]
    z_source = galaxy_redshift[index]
    sigma_v = Sigma_v_set[index]
    q = Q_set[index]
    y1 = y1_set[index]
    y2 = y2_set[index]
    # lensing quantities
    ckm = c.value * (u.m/u.s).to(u.km/u.s) #千米每秒单位下的光速
    D_ds = Planck18.angular_diameter_distance_z1z2(z_lens, z_source).value
    D_s = Planck18.angular_diameter_distance(z_source).value
    theta_E = 4 * np.pi * sigma_v ** 2 / ckm ** 2 * D_ds / D_s * 3600 * 180 / np.pi #以角秒为单位
    e1, e2 = param_util.phi_q2_ellipticity(phi=np.pi/2, q=q)
    
    source_x = y1 * theta_E
    source_y = y2 * theta_E
    kwargs_lens = [
        {'theta_E': theta_E, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}  # SIE model
    ]
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
    kwargs_model = {'lens_model_list': ['SIE'],  # list of lens models to be used
                    'lens_light_model_list': ['SERSIC_ELLIPSE'],  # list of unlensed light models to be used
                    'source_light_model_list': ['SERSIC_ELLIPSE'],  # list of extended source models to be used
        }
    fixed_lens = []
    fixed_lens.append({})
    fixed_lens.append({'ra_0': 0, 'dec_0': 0})
    fixed_source = []
    fixed_source.append({})
    fixed_lens_light = []
    fixed_lens_light.append({})
    kwargs_constraints = {'linear_solver': True}
    param = Param(kwargs_model, fixed_lens, fixed_source, fixed_lens_light, 
                kwargs_lens_init=kwargs_result['kwargs_lens'], **kwargs_constraints)


    def TD_consistent(index_start, index_end):
        SFR_set = []
        x_source_result = []
        y_source_result = []
        theta_E_result = []
        center_x_result = []
        center_y_result = []
        e1_result = []
        e2_result = []
        t_days_0_result = []
        t_days_1_result = []
        t_days_2_result = []
        t_days_3_result = []
        for index_in_mcmc in tqdm(range(index_start, index_end)):
            # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
            dis_x_set = []
            dis_y_set = []
            
            kwargs_result_tmp = param.args2kwargs(samples_mcmc[index_in_mcmc])
            '''透镜质量参数'''
            theta_E_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['theta_E']
            center_x_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['center_x']
            center_y_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['center_y']
            e1_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['e1']
            e2_result_tmp = kwargs_result_tmp['kwargs_lens'][0]['e2']
            kwargs_sie_result_tmp = {'theta_E': theta_E_result_tmp, 'center_x': center_x_result_tmp, 'center_y': center_y_result_tmp, 'e1': e1_result_tmp, 'e2': e2_result_tmp}  # parameters of the deflector lens model
            kwargs_lens_result_tmp = [kwargs_sie_result_tmp]
            '''源参数'''
            e1_source_light_result_tmp = kwargs_result_tmp['kwargs_source'][0]['e1']
            e2_source_light_result_tmp = kwargs_result_tmp['kwargs_source'][0]['e2']
            n_sersic_result_tmp = kwargs_result_tmp['kwargs_source'][0]['n_sersic'] 
            R_cir_result_tmp = kwargs_result_tmp['kwargs_source'][0]['R_sersic']
            x_source_center_result_tmp = kwargs_result_tmp['kwargs_source'][0]['center_x'] 
            y_source_center_result_tmp = kwargs_result_tmp['kwargs_source'][0]['center_y'] 
            
            source_phi, source_q = param_util.ellipticity2phi_q(e1_source_light_result_tmp, e2_source_light_result_tmp)
            mod = Sersic2D(amplitude=1, r_eff=R_cir_result_tmp, n=n_sersic_result_tmp, x_0=0, y_0=0,
               ellip=source_q, theta=source_phi)
            
            Re_maj_result_tmp = R_cir_result_tmp / np.sqrt(source_q) #major half-light radius
            Re_min_result_tmp = source_q * Re_maj_result_tmp
            dis_x_max = np.max([np.abs(Re_maj_result_tmp * np.cos(source_phi)), np.abs(Re_min_result_tmp * np.sin(source_phi))])+ 1
            dis_y_max = np.max([np.abs(Re_maj_result_tmp * np.sin(source_phi)), np.abs(Re_min_result_tmp * np.cos(source_phi))])+ 1 
            for tmp_i in range(1000):
                dis_x = random.uniform(-dis_x_max, dis_x_max)
                dis_y = random.uniform(-dis_y_max, dis_y_max)
                while (dis_x * np.cos(source_phi) + dis_y * np.sin(source_phi))**2 / Re_maj_result_tmp ** 2 + (-dis_x * np.sin(source_phi) + dis_y * np.cos(source_phi))**2 / Re_min_result_tmp ** 2 > 1:
                    dis_x = random.uniform(-dis_x_max, dis_x_max)
                    dis_y = random.uniform(-dis_y_max, dis_y_max)
            
                x_source_result_tmp, y_source_result_tmp = x_source_center_result_tmp + dis_x, y_source_center_result_tmp + dis_y
            
            
                        
                x_image_result_tmp, y_image_result_tmp = lensEquationSolver.findBrightImage(x_source_result_tmp, y_source_result_tmp, kwargs_lens_result_tmp, numImages=4,
                                                                min_distance=0.5*kwargs_F200W_band['pixel_scale'], search_window=numpix * kwargs_F200W_band['pixel_scale'])
                ps_amp_result_tmp = [1] * len(x_image_result_tmp) 
                if len(x_image_result_tmp) == 4:
                    dis_x_set.append(dis_x)
                    dis_y_set.append(dis_y)
                    kwargs_ps_result_tmp = [{'ra_image': x_image_result_tmp, 'dec_image': y_image_result_tmp,
                                'point_amp': ps_amp_result_tmp}]  # quasar point source position in the source plane and intrinsic brightness
                    
                    t_days = td_cosmo.time_delays(kwargs_lens_result_tmp, kwargs_ps_result_tmp, kappa_ext=0)
                    
                    x_source_result.append(x_source_result_tmp)
                    y_source_result.append(y_source_result_tmp)
                    theta_E_result.append(theta_E_result_tmp)
                    center_x_result.append(center_x_result_tmp)
                    center_y_result.append(center_y_result_tmp)
                    e1_result.append(e1_result_tmp)
                    e2_result.append(e2_result_tmp)
                    t_days_0_result.append(t_days[0])
                    t_days_1_result.append(t_days[1])
                    t_days_2_result.append(t_days[2])
                    t_days_3_result.append(t_days[3])
                    SFR_set.append(mod(dis_x, dis_y))

        mcmc_new_list = []
        for i in range(len(theta_E_result)):
            mcmc_new_list.append([theta_E_result[i], center_x_result[i], center_y_result[i], e1_result[i], e2_result[i], t_days_0_result[i], t_days_1_result[i], t_days_2_result[i], t_days_3_result[i], SFR_set[i]])
        
        return mcmc_new_list
                    


    # Generate_file(0,1)
    if __name__=='__main__':
        
        threadcount = 400
        res_z = [0 for i in range(threadcount)]
        length = len(samples_mcmc) // 20
        
        
        percount = length//threadcount

        print("percount = ", percount)
    
        # t0 = time.time()
        pool = Pool(threadcount) 
        # for i in range(len(m1_set)):
        for i in range(threadcount):
            if i < threadcount - 1:
                res = pool.apply_async(func=TD_consistent, args=(percount * i, percount * (i + 1),))
                res_z[i] = res
                # print("0",i)
            else:
                # print("1",i)
                res = pool.apply_async(func=TD_consistent, args=(percount * i, length,)) 
                res_z[i] = res

        pool.close()
        pool.join() 
        # Result1 = list(res_z[0].get())
        # Result1 = list(res_z[0].get())
        Result1  = list(res_z[0].get())
        for i in range(threadcount-1):
            Result1.extend(list(res_z[i+1].get()))
            
        

    x_image_truths, y_image_truths = lensEquationSolver.findBrightImage(source_x, source_y, kwargs_lens, numImages=4,
                                                        min_distance=0.01*kwargs_F200W_band['pixel_scale'], search_window=numpix * kwargs_F200W_band['pixel_scale'])
    kwargs_ps_truths = [{'ra_image': x_image_truths, 'dec_image': y_image_truths,
                    'point_amp': [1] * len(x_image_truths)}]  # quasar point source position in the source plane and intrinsic brightness

    t_days_truths = td_cosmo.time_delays(kwargs_lens, kwargs_ps_truths, kappa_ext=0)

    labels_new = [r"$theta_E$", r"$center_x$", r"$center_y$", r"$e_1$", r"$e_2$", r"$t_0$", r"$t_1$", r"$t_2$", r"t_3", r"SFR"]
    try:
        plot = corner.corner(np.array(Result1), labels=labels_new, show_titles=True, truths=[kwargs_lens[0]['theta_E'], kwargs_lens[0]['center_x'], kwargs_lens[0]['center_y'], kwargs_lens[0]['e1'], kwargs_lens[0]['e2'], t_days_truths[0], t_days_truths[1], t_days_truths[2], t_days_truths[3], 0])
        plt.savefig('./' + JWST_file + '/' + str(index) + 'corner_Res.png', dpi=450)
        plt.close()
        np.savetxt('./' + JWST_file + '/' + str(index) + 'corner_Res_data.csv', Result1, delimiter=',')
    except AssertionError:
        pass
