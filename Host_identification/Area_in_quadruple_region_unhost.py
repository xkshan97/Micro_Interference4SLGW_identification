"""前半段计算宿主星系的半光半径面积以及能成四像区域的面积，单位是kpc^2，注意转换，因为JWST_Res中储存的是角秒"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
import astropy.units as u
from astropy.constants import c
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
import astropy.cosmology.units as cu #//FIXME
from tqdm import tqdm
from MySample import sample
from scipy import special as sp
import random
import lenstronomy.Util.param_util as param_util
from  multiprocessing import Process,Pool

cosmo = Planck18


"""下面计算非宿主星系的面积"""
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

# S_unhost_quadruple = []
def cal_area(index):
# for index in Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images:
    index = int(index)
    count_in_elliptical = 0
    count_in_quadruple = 0
    z_lens = lens_z[index]
    z_source = galaxy_redshift[index]
    angular_diameter_distance_to_unhost = Planck18.angular_diameter_distance(z_source).to(u.kpc).value
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

    source_x_0 = y1 * theta_E
    source_y_0 = y2 * theta_E

    lens_model_list = ['SIE']
    #initialize classes with the specific lens models
    lensModel = LensModel(lens_model_list, cosmo=cosmo, z_lens=z_lens, z_source=z_source)
    lensEquationSolver = LensEquationSolver(lensModel)
    # chose a lens model parameterization
    kwargs_sie = [{'theta_E': theta_E, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}]
    
    
    
    source_phi=source_position_angle[index] * np.pi / 180
    source_q=source_axis_ratio[index]
    Re_circ = source_Re_circ_arc_sec[index] #circularized half-light radius
    Re_maj = Re_circ / np.sqrt(source_q) #major half-light radius
    Re_min = source_q * Re_maj
    S_elliptical = np.pi * Re_maj * Re_min * angular_diameter_distance_to_unhost ** 2 * (u.arcsec.to(u.rad)) ** 2 #以kpc^为单位。
    dis_x_max = np.max([np.abs(Re_maj * np.cos(source_phi)), np.abs(Re_min * np.sin(source_phi))])+ 0.1
    dis_y_max = np.max([np.abs(Re_maj * np.sin(source_phi)), np.abs(Re_min * np.cos(source_phi))])+ 0.1 
    for test_i in tqdm(range(1000)):
        dis_x = random.uniform(-dis_x_max, dis_x_max)
        dis_y = random.uniform(-dis_y_max, dis_y_max)
        test_in_ellip = (dis_x * np.cos(source_phi) + dis_y * np.sin(source_phi))**2 / Re_maj ** 2 + (-dis_x * np.sin(source_phi) + dis_y * np.cos(source_phi))**2 / Re_min ** 2
        if test_in_ellip < 1:
            count_in_elliptical += 1
            x_source, y_source = source_x_0 + dis_x, source_y_0 + dis_y
            x_image, y_image = lensEquationSolver.image_position_from_source(kwargs_lens=kwargs_sie, sourcePos_x=x_source, 
                                                                                sourcePos_y=y_source, min_distance=0.005, search_window=5, 
                                                                                    precision_limit=10**(-10), num_iter_max=100)
            if len(x_image) == 4:
                count_in_quadruple += 1
    
    S_quadruple = count_in_quadruple / count_in_elliptical * S_elliptical
    # S_unhost_quadruple.append(S_quadruple)
    return S_quadruple


# Generate_file(0,1)
if __name__=='__main__':
    
    threadcount = len(Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images)
    res_z = [0 for i in range(threadcount)]
    length = threadcount
    
    percount = length//threadcount

    print("percount = ", percount)

    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    i = 0
    for index in Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images:
    
        res = pool.apply_async(func=cal_area, args=(index,))
        res_z[i] = res
        i += 1

    pool.close()
    pool.join() 
    
    S_unhost_quadruple = []
    for i in range(threadcount):
        S_unhost_quadruple.append(res_z[i].get())
    
np.savetxt('./SampleResult_Galaxy/S_unhost_quadruple.csv', S_unhost_quadruple, delimiter=',')

