#将引力波的宿主星系中心位置随机摆放，保证中心在半光半径之内。
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from tqdm import tqdm
from astropy.cosmology import Planck18
import astropy.units as u
from astropy.constants import c
import astropy.cosmology.units as cu #//FIXME
from tqdm import tqdm
from MySample import sample
from scipy import special as sp
import random
import lenstronomy.Util.param_util as param_util
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from  multiprocessing import Process,Pool
from Lens_light_population import Sample_u_R2mag_i

cosmo = Planck18
random.seed(1000)

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

source_x_0 = y1 * theta_E
source_y_0 = y2 * theta_E

lens_model_list = ['SIE']
#initialize classes with the specific lens models
lensModel = LensModel(lens_model_list, cosmo=cosmo, z_lens=z_lens, z_source=z_source)
lensEquationSolver = LensEquationSolver(lensModel)
# chose a lens model parameterization
kwargs_sie = [{'theta_E': theta_E, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}]

# for source_index in range(len(JWST_mag_Res)):

dis_x_set = []
dis_y_set = []
source_x = []
source_y = []
t_days_set = []
test_num = 0
source_index = 25
for j in tqdm(range(1000)):
    sum_area = [0]*len(JWST_mag_Res) #用来存放引力波的位置是否都在40个宿主星系的半光半径内
    # choose a source position
    source_phi=JWST_mag_Res[source_index][2] * np.pi / 180
    source_q=JWST_mag_Res[source_index][1]
    Re_circ = JWST_mag_Res[source_index][3] #circularized half-light radius
    Re_maj = Re_circ / np.sqrt(source_q) #major half-light radius
    Re_min = source_q * Re_maj
    dis_x_max = np.max([np.abs(Re_maj * np.cos(source_phi)), np.abs(Re_min * np.sin(source_phi))])+ 1
    dis_y_max = np.max([np.abs(Re_maj * np.sin(source_phi)), np.abs(Re_min * np.cos(source_phi))])+ 1 
    dis_x = random.uniform(-dis_x_max, dis_x_max)
    dis_y = random.uniform(-dis_y_max, dis_y_max)
    source_index_i = source_index
    
    while (dis_x * np.cos(source_phi) + dis_y * np.sin(source_phi))**2 / Re_maj ** 2 + (-dis_x * np.sin(source_phi) + dis_y * np.cos(source_phi))**2 / Re_min ** 2 > 1:
        source_phi=JWST_mag_Res[source_index][2] * np.pi / 180
        source_q=JWST_mag_Res[source_index][1]
        Re_circ = JWST_mag_Res[source_index][3] #circularized half-light radius
        Re_maj = Re_circ / np.sqrt(source_q) #major half-light radius
        Re_min = source_q * Re_maj
        dis_x_max = np.max([np.abs(Re_maj * np.cos(source_phi)), np.abs(Re_min * np.sin(source_phi))])+ 1
        dis_y_max = np.max([np.abs(Re_maj * np.sin(source_phi)), np.abs(Re_min * np.cos(source_phi))])+ 1 
        dis_x = random.uniform(-dis_x_max, dis_x_max)
        dis_y = random.uniform(-dis_y_max, dis_y_max)
        
    dis_x_set.append(dis_x)
    dis_y_set.append(dis_y)
    x_source, y_source = source_x_0 + dis_x, source_y_0 + dis_y
    x_image, y_image = lensEquationSolver.image_position_from_source(kwargs_lens=kwargs_sie, sourcePos_x=x_source, 
                                                                        sourcePos_y=y_source, min_distance=0.005, search_window=5, 
                                                                            precision_limit=10**(-10), num_iter_max=100)
    if len(x_image) < 4:
        test_num += 1
    
    else:
        mag_ = lensModel.magnification(x_image, y_image, kwargs_sie)
        t_days = lensModel.arrival_time(x_image, y_image, kwargs_sie)
        kappa = lensModel.kappa(x_image, y_image, kwargs_sie)
        gamma1, gamma2 = lensModel.gamma(x_image, y_image, kwargs_sie)
        gamma = np.sqrt(gamma1**2+ gamma2**2)
        
        source_x.append(x_source)
        source_y.append(y_source)
        t_days_set.append(t_days)
        
np.savetxt('./test/source_x.csv', source_x, delimiter=',')
np.savetxt('./test/source_y.csv', source_y, delimiter=',') 
# np.savetxt('Random_pin_down_GW/x_image_' + str(source_index) + '_' + str(j) + '_.csv', x_image)
# np.savetxt('Random_pin_down_GW/y_image_' + str(source_index) + '_' + str(j) + '_.csv', y_image)
# np.savetxt('Random_pin_down_GW/mag_' + str(source_index) + '_' + str(j) + '_.csv', mag_)
np.savetxt('./test/t_days.csv', t_days_set, delimiter=',')

# np.savetxt('Random_pin_down_GW/dis_x_set_' + str(source_index) + '_.csv', dis_x_set)
# np.savetxt('Random_pin_down_GW/dis_y_set_' + str(source_index) + '_.csv', dis_y_set)



#caustic
phi = np.arange(0, 2*np.pi, 0.01)
delta_phi = np.sqrt(np.cos(phi)**2 + q**2 * np.sin(phi)**2)
f_prime = np.sqrt(1 - q**2)
y1 = np.sqrt(q) / delta_phi * np.cos(phi) - np.sqrt(q) / f_prime * np.arcsinh(f_prime/q*np.cos(phi))
y2 = np.sqrt(q) / delta_phi * np.sin(phi) - np.sqrt(q) / f_prime * np.arcsin(f_prime*np.sin(phi))


#半光半径 
theta = np.arange(0, 2*np.pi, 0.01)
r = np.sqrt(1 / ((np.cos(theta) * np.cos(source_phi) + np.sin(theta) * np.sin(source_phi))**2 / Re_maj ** 2 + (-np.cos(theta) * np.sin(source_phi) + np.sin(theta) * np.cos(source_phi))**2 / Re_min ** 2))
x = source_x_0 + r * np.cos(theta)
y = source_y_0 + r * np.sin(theta)
np.savetxt('./test/y1.csv', y1, delimiter=',')
np.savetxt('./test/y2.csv', y2, delimiter=',') 
np.savetxt('./test/x.csv', x, delimiter=',')
np.savetxt('./test/y.csv', y, delimiter=',') 