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
cosmo = Planck18
from  multiprocessing import Process,Pool

"""下面计算宿主星系的面积"""
JWST_file = 'JWST_host_pop_little_sigma'

SampleParam = np.loadtxt('../Paper4_CE_Modify/SampleResult/SampleParameter.csv', delimiter=',')
#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')
lens_z = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')
Sigma_v_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/sigma_v.csv', delimiter=',')
Q_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/axis_ratio.csv', delimiter=',')
y1_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/y1.csv', delimiter=',')
y2_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/y2.csv', delimiter=',')

index = 120
z_lens = lens_z[index]
angular_diameter_distance_to_lens = Planck18.angular_diameter_distance(z_lens).to(u.kpc).value
z_source = SampleParam[0][int(AcceptLensIndex[index])]
angular_diameter_distance_to_host = Planck18.angular_diameter_distance(z_source).to(u.kpc).value
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


def cal_area(source_index_i):
# for source_index_i in range(len(JWST_mag_Res)):
    count_in_elliptical = 0
    count_in_quadruple = 0
    source_phi=JWST_mag_Res[source_index_i][2] * np.pi / 180
    source_q=JWST_mag_Res[source_index_i][1]
    Re_circ = JWST_mag_Res[source_index_i][3] #circularized half-light radius
    Re_maj = Re_circ / np.sqrt(source_q) #major half-light radius
    Re_min = source_q * Re_maj
    S_elliptical = np.pi * Re_maj * Re_min * angular_diameter_distance_to_host ** 2 * (u.arcsec.to(u.rad)) ** 2 #以kpc^为单位。
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
    # S_host_quadruple.append(S_quadruple)
    return S_quadruple



# Generate_file(0,1)
if __name__=='__main__':
    
    threadcount = len(JWST_mag_Res)
    res_z = [0 for i in range(threadcount)]
    length = threadcount
    
    percount = length//threadcount

    print("percount = ", percount)

    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
    
        res = pool.apply_async(func=cal_area, args=(i,))
        res_z[i] = res

    pool.close()
    pool.join() 
    
    S_host_quadruple = []
    for i in range(threadcount):
        S_host_quadruple.append(res_z[i].get())
        
np.savetxt('./SampleResult_GWGalaxy/S_host_quadruple.csv', S_host_quadruple, delimiter=',')




