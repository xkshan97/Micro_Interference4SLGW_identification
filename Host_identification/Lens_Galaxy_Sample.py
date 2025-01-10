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


hdu=fits.open('JWST_catalogue/JADES_SF_mock_r1_v1.2.fits')
name = hdu[1].data.names
galaxy_redshift = hdu[1].data['redshift']


tau = []
acceptIndex = []
randomset = []
for i in range(len(galaxy_redshift)):
    #按照徐菲的文章，把SIS的多像optical depth转换成了SIE的 
    tau.append((1 - np.exp(-4.17*10**(-6)*(Planck18.comoving_distance(galaxy_redshift[i]).to(u.Gpc).value)**3)) * 2.69)

np.random.seed(1000000)

apparent_mag_set = []
JWST_apparent_mag_set = []
for j in tqdm(range(300)): #上面的星表是11 * 11 arcmin^2的，转换成10平方度内的透镜星系，所以乘以300，注意这里的7平方度是挑选出来的3个四像引力波事例定位最大的天区。
    for i in range(len(tau)):
        tmp = np.random.random()
        if tau[i] > tmp and i not in acceptIndex: #按道理说这里应该是没有i not in acceptIndex吧
            acceptIndex.append(i)
            randomset.append(tmp)
            
            HST_F606W_fnu_tmp = hdu[1].data['HST_F606W_fnu'][i]
            NRC_F200W_fnu_tmp = hdu[1].data['NRC_F200W_fnu'][i] 

            HST_flux_tmp = HST_F606W_fnu_tmp * 1.0000000000000002e-32 #转成erg/s/cm^2/Hz
            apparent_mag_tmp = -2.5 * np.log10(HST_flux_tmp) - 48.6
            apparent_mag_set.append(apparent_mag_tmp)

            JWST_flux_tmp = NRC_F200W_fnu_tmp * 1.0000000000000002e-32 #转成erg/s/cm^2/Hz
            JWST_apparent_mag_tmp = -2.5 * np.log10(JWST_flux_tmp) - 48.6
            JWST_apparent_mag_set.append(JWST_apparent_mag_tmp)
            
        else:
            pass

where_are_inf = np.isinf(apparent_mag_set)
apparent_mag_set = np.array(apparent_mag_set)
apparent_mag_set[where_are_inf] = 100


where_are_inf = np.isinf(JWST_apparent_mag_set)
JWST_apparent_mag_set = np.array(JWST_apparent_mag_set)
JWST_apparent_mag_set[where_are_inf] = 100


np.savetxt('./SampleResult_Galaxy/AcceptLensIndex.csv', acceptIndex, delimiter=',')
np.savetxt('./SampleResult_Galaxy/tau.csv', tau, delimiter=',')
np.savetxt('./SampleResult_Galaxy/AcceptRandomValue.csv', randomset, delimiter=',')
np.savetxt('./SampleResult_Galaxy/ApparentMag.csv', apparent_mag_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/JWST_ApparentMag.csv', JWST_apparent_mag_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/GalaxySourceLensedRedshift.csv', galaxy_redshift[acceptIndex], delimiter=',')
np.savetxt('./SampleResult_Galaxy/GalaxySourceAxis_ratio.csv', hdu[1].data['axis_ratio'][acceptIndex], delimiter=',')
np.savetxt('./SampleResult_Galaxy/GalaxySourcePosition_angle.csv', hdu[1].data['position_angle'][acceptIndex], delimiter=',')
np.savetxt('./SampleResult_Galaxy/GalaxySourceRe_circ.csv', hdu[1].data['Re_circ'][acceptIndex], delimiter=',')
np.savetxt('./SampleResult_Galaxy/GalaxySourceSersic_n.csv', hdu[1].data['Sersic_n'][acceptIndex], delimiter=',')






def Lens_Dc(x):
    return 30 * x**2 * (1 - x)**2

def SampleLens_Dc():
    return sample(Lens_Dc, 1, lower_bd=0, upper_bd=1, guess=0.5)[0]
    


def GetLens_z(Dc_s):
    r =SampleLens_Dc()
    Dc_l = r * Dc_s
    Dc_l = Dc_l * u.Mpc
    z_l = Dc_l.to(cu.redshift, cu.redshift_distance(Planck18, kind="comoving"))
    return z_l.value

def Parameter_a(x):
    #范围只需要从0到2.5就行
    alpha = 2.32
    beta = 2.67
    return x ** (alpha - 1) * np.exp(-x**beta) * beta / sp.gamma(alpha/beta)

def Parameter_b(x, a):
    #经过验证，a越大，分布越宽，但是x最大为2.5也足够了。
    #范围是0到0.8，这是haris文章中要求的，附录中的第四步，until we get a sample b < 0.8。
    s = 0.38 - 0.09177 * a #参考文献：Beyond the Detector Horizon: Forecasting Gravitational-Wave Strong Lensing
    return x/s**2 * np.exp(-x**2/2/s**2)

def Get_axis_ratio_and_sigma_v():
    a = sample(Parameter_a, 1, lower_bd=0, upper_bd=2.5, guess=0.7)[0]
    b = sample(Parameter_b, 1, lower_bd=0, upper_bd=0.8, guess=0.15, args=(a,))[0]
    while b < 0 or b > 0.8:
        a = sample(Parameter_a, 1, lower_bd=0, upper_bd=2.5, guess=0.7)[0]
        b = sample(Parameter_b, 1, lower_bd=0, upper_bd=0.8, guess=0.15, args=(a,))[0]  
    sigma_v = 161 * a
    q = 1 - b
    return sigma_v, q

def Parameter_y1(q):
    if q == 1:
        return 0
    else:
        return random.uniform(-(q/(1 - q**2))**0.5 * np.arccosh(1/q), (q/(1 - q**2))**0.5 * np.arccosh(1/q))

def Parameter_y2(q):
    if q > 0.3942:
        if q == 1:
            return 0
        else:
            return random.uniform(-(q/(1 - q**2))**0.5 * np.arccos(q), (q/(1 - q**2))**0.5 * np.arccos(q))
    else:
        return random.uniform((q/(1 - q**2))**0.5 * np.arccos(q) - (1/q)**0.5, (1/q)**0.5 - (q/(1 - q**2))**0.5 * np.arccos(q))

#=======================================
#lenstronomy求解

def Solve_Lens_equation_lenstronomy(y1, y2, q, z_s, z_l, sigma_v):
    # chose a lens model (list of parameterized lens models)
    z_lens = z_l
    z_source = z_s
    ckm = c.value * (u.m/u.s).to(u.km/u.s) #千米每秒单位下的光速
    D_ds = Planck18.angular_diameter_distance_z1z2(z_l, z_s).value
    D_s = Planck18.angular_diameter_distance(z_s).value
    theta_E = 4 * np.pi * sigma_v ** 2 / ckm ** 2 * D_ds / D_s * 3600 * 180 / np.pi #以角秒为单位
    cosmo = Planck18
    
    lens_model_list = ['SIE']
    e1, e2 = param_util.phi_q2_ellipticity(phi=np.pi/2, q=q)
    
    #initialize classes with the specific lens models
    lensModel = LensModel(lens_model_list, cosmo=cosmo, z_lens=z_lens, z_source=z_source)
    lensEquationSolver = LensEquationSolver(lensModel)
    # chose a lens model parameterization
    kwargs_sie = [{'theta_E': theta_E, 'e1': e1, 'e2': e2, 'center_x': 0, 'center_y': 0}]

    # choose a source position
    x_source, y_source = y1 * theta_E, y2 * theta_E
    x_image, y_image = lensEquationSolver.image_position_from_source(kwargs_lens=kwargs_sie, sourcePos_x=x_source, 
                                                                         sourcePos_y=y_source, min_distance=0.005, search_window=5, 
                                                                            precision_limit=10**(-10), num_iter_max=100)
    
    mag_ = lensModel.magnification(x_image, y_image, kwargs_sie)
    t_days = lensModel.arrival_time(x_image, y_image, kwargs_sie)
    kappa = lensModel.kappa(x_image, y_image, kwargs_sie)
    gamma1, gamma2 = lensModel.gamma(x_image, y_image, kwargs_sie)
    gamma = np.sqrt(gamma1**2+ gamma2**2)
    # print('kappa: ', kappa)
    # print('gamma: ', gamma)
    # print('time delays: ', t_days)
    return x_image, y_image, mag_, t_days, kappa, gamma, theta_E
    



def LensCal(start_index, end_index):
    

    Lens_z_set_tmp = []
    Q_set_tmp = []
    Sigma_v_set_tmp = []
    y1_set_tmp = []
    y2_set_tmp = []
    mag_set_tmp = []
    t_days_set_tmp = []
    kappa_set_tmp = []
    gamma_set_tmp = []
    x_image_set_tmp = []
    y_image_set_tmp = []
    theta_E_tmp = []
    Lens_light_mag_tmp = []
    Lens_R_e_tmp = []
    
    for i_ in tqdm(range(start_index, end_index)):
        z_s = galaxy_redshift[int(acceptIndex[i_])]
        
        comoving_dis_zs = Planck18.comoving_distance(z_s).value
        lens_z = GetLens_z(comoving_dis_zs)
        sigma_v, q = Get_axis_ratio_and_sigma_v()
        y1 = Parameter_y1(q)
        y2 = Parameter_y2(q)
        
        x_image, y_image, mag_, t_days, kappa, gamma, theta_E = Solve_Lens_equation_lenstronomy(y1, y2, q, z_s, lens_z, sigma_v)
        while len(x_image) == 0 or len(x_image) == 1 or len(x_image) == 3:
            lens_z = GetLens_z(comoving_dis_zs)
            sigma_v, q = Get_axis_ratio_and_sigma_v()
            y1 = Parameter_y1(q)
            y2 = Parameter_y2(q)
            x_image, y_image, mag_, t_days, kappa, gamma, theta_E = Solve_Lens_equation_lenstronomy(y1, y2, q, z_s, lens_z, sigma_v)
        
        
        Lens_z_set_tmp.append(lens_z)
        Q_set_tmp.append(q)
        Sigma_v_set_tmp.append(sigma_v)
        y1_set_tmp.append(y1)
        y2_set_tmp.append(y2)
        theta_E_tmp.append(theta_E)
        
        mag_set_tmp.append(list(mag_))
        t_days_set_tmp.append(list(t_days))
        kappa_set_tmp.append(list(kappa))
        gamma_set_tmp.append(list(gamma))
        x_image_set_tmp.append(list(x_image))
        y_image_set_tmp.append(list(y_image))
        
        Sample_u_R2mag_i_tmp = Sample_u_R2mag_i(sigma_v, lens_z)
        Lens_light_mag_tmp.append(Sample_u_R2mag_i_tmp[0].value)
        Lens_R_e_tmp.append(Sample_u_R2mag_i_tmp[1])
        
    
    return Lens_z_set_tmp, Q_set_tmp, Sigma_v_set_tmp, y1_set_tmp, y2_set_tmp, mag_set_tmp, t_days_set_tmp, kappa_set_tmp, gamma_set_tmp, x_image_set_tmp, y_image_set_tmp, theta_E_tmp, Lens_light_mag_tmp, Lens_R_e_tmp
    


# Generate_file(0,1)
if __name__=='__main__':
    # pool = Pool(len(m1_set)) #创建一个5个进程的进程池
    threadcount = 380
    res_z = [0 for i in range(threadcount)]
    length = len(acceptIndex)
    
     
    percount = length//threadcount

    print("percount = ", percount)
   
    # t0 = time.time()
    pool = Pool(threadcount) 
    # for i in range(len(m1_set)):
    for i in range(threadcount):
        if i < threadcount - 1:
            res = pool.apply_async(func=LensCal, args=(percount * i, percount * (i + 1),))
            res_z[i] = res
            # print("0",i)
        else:
            # print("1",i)
            res = pool.apply_async(func=LensCal, args=(percount * i, length,)) 
            res_z[i] = res

    pool.close()
    pool.join()
    

Result1, Result2, Result3, Result4, Result5, Result6, Result7, Result8, Result9, Result10, Result11, Result12, Result13, Result14  = list(res_z[0].get()[0]), list(res_z[0].get()[1]), list(res_z[0].get()[2]), list(res_z[0].get()[3]), list(res_z[0].get()[4]), list(res_z[0].get()[5]), list(res_z[0].get()[6]), list(res_z[0].get()[7]), list(res_z[0].get()[8]), list(res_z[0].get()[9]), list(res_z[0].get()[10]), list(res_z[0].get()[11]), list(res_z[0].get()[12]), list(res_z[0].get()[13])
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
    Result10.extend(list(res_z[i+1].get()[9]))
    Result11.extend(list(res_z[i+1].get()[10]))
    Result12.extend(list(res_z[i+1].get()[11]))
    Result13.extend(list(res_z[i+1].get()[12]))
    Result14.extend(list(res_z[i+1].get()[13]))
    

Lens_z_set = Result1
Q_set = Result2
Sigma_v_set = Result3
y1_set = Result4
y2_set = Result5
mag_set = Result6
t_days_set = Result7
kappa_set = Result8
gamma_set = Result9
x_image_set = Result10
y_image_set = Result11
theta_E_set = Result12
Lens_light_mag = Result13
Lens_R_e = Result14

image_num = []
for i in range(len(mag_set)):
    image_num.append(len(mag_set[i]))
    
mag_set = [y for x in mag_set for y in x]
t_days_set = [y for x in t_days_set for y in x]
kappa_set = [y for x in kappa_set for y in x]
gamma_set = [y for x in gamma_set for y in x]
x_image_set = [y for x in x_image_set for y in x]
y_image_set = [y for x in y_image_set for y in x]



np.savetxt('./SampleResult_Galaxy/magnification.csv', mag_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/timedelayDays.csv', t_days_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/kappa.csv', kappa_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/gamma.csv', gamma_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/imagenumber.csv', image_num, delimiter=',')
np.savetxt('./SampleResult_Galaxy/lens_z.csv', Lens_z_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/axis_ratio.csv', Q_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/sigma_v.csv', Sigma_v_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/y1.csv', y1_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/y2.csv', y2_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/x1_image.csv', x_image_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/x2_image.csv', y_image_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/theta_E.csv', theta_E_set, delimiter=',')
np.savetxt('./SampleResult_Galaxy/Lens_light_mag.csv', Lens_light_mag, delimiter=',')
np.savetxt('./SampleResult_Galaxy/Lens_R_e.csv', Lens_R_e, delimiter=',')



#根据上面的结果画图（Collett 2015的图）
import matplotlib.pyplot as plt
theta_E_set = np.loadtxt('./SampleResult_Galaxy/theta_E.csv', delimiter=',')
Lens_z_set = np.loadtxt('./SampleResult_Galaxy/lens_z.csv', delimiter=',')
apparent_mag_set = np.loadtxt('./SampleResult_Galaxy/ApparentMag.csv', delimiter=',')
JWST_apparent_mag_set = np.loadtxt('./SampleResult_Galaxy/JWST_ApparentMag.csv', delimiter=',')
Sigma_v_set = np.loadtxt('./SampleResult_Galaxy/sigma_v.csv', delimiter=',')
acceptIndex_tmp = np.loadtxt('./SampleResult_Galaxy/AcceptLensIndex.csv', delimiter=',')
acceptIndex = []
for i in range(len(acceptIndex_tmp)):
    acceptIndex.append(int(acceptIndex_tmp[i]))
acceptIndex = np.array(acceptIndex)
source_galaxy_redshift_lensed = np.loadtxt('./SampleResult_Galaxy/GalaxySourceLensedRedshift.csv', delimiter=',')
image_num = np.loadtxt('./SampleResult_Galaxy/imagenumber.csv', delimiter=',')

index_mag_less_27_and_sigma_v_gtr_100 = np.where((apparent_mag_set < 27)&(Sigma_v_set>100))[0]


# #PDF
fig, ax = plt.subplots(1, 3, figsize=(12, 4))

ax[0].hist(theta_E_set[index_mag_less_27_and_sigma_v_gtr_100], density=True, bins=50)
ax[0].set_xlim(0, 3)
ax[0].set_xlabel('$\\theta_E$')
ax[0].grid()
ax[0].semilogy()


ax[1].hist(Lens_z_set[index_mag_less_27_and_sigma_v_gtr_100], density=True, bins=100)
ax[1].hist(source_galaxy_redshift_lensed[index_mag_less_27_and_sigma_v_gtr_100], density=True, bins=100)
ax[1].set_xlim(0, 5)
ax[1].set_xlabel('$z$')
ax[1].grid()

ax[2].hist(Sigma_v_set[index_mag_less_27_and_sigma_v_gtr_100], density=True, bins=50)
ax[2].set_xlim(100, 400)
ax[2].set_xlabel('$\\sigma_v$')
ax[2].grid()

plt.savefig('./IntermediaPlot/Lens_Galaxy_Sample_collett_fig1.png', dpi=450)
plt.close()

#找出来这个里面星等大于26，然后爱因斯坦半径大于seeing和R_e的，然后成四像的
index_mag_less_26 = np.where((apparent_mag_set < 26))[0]

hdu=fits.open('JWST_catalogue/JADES_SF_mock_r1_v1.2.fits')
name = hdu[1].data.names
R_e_set = hdu[1].data['Re_circ'][acceptIndex]
angular_diameter_distance_to_source = Planck18.angular_diameter_distance(source_galaxy_redshift_lensed).to(u.kpc).value
R_e_arc_sec = R_e_set / angular_diameter_distance_to_source * u.rad.to(u.arcsec)
seeing = 0.135
Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images = np.where((theta_E_set**2 - R_e_arc_sec ** 2 - (seeing / 2) ** 2 > 0)&(apparent_mag_set < 26)&(image_num==4))[0]
JWST_Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_28p8_and_4_images = np.where((theta_E_set**2 - R_e_arc_sec ** 2 - (seeing / 2) ** 2 > 0)&(JWST_apparent_mag_set < 28.8)&(image_num==4))[0]

R_e_set_Lens = np.loadtxt('./SampleResult_Galaxy/Lens_R_e.csv', delimiter=',')
angular_diameter_distance_to_lens = Planck18.angular_diameter_distance(Lens_z_set).to(u.kpc).value
R_e_Lens_arc_sec = R_e_set_Lens / angular_diameter_distance_to_lens * u.rad.to(u.arcsec)


np.savetxt('SampleResult_Galaxy/R_e_arc_sec.csv', R_e_arc_sec, delimiter=',')

np.savetxt('SampleResult_Galaxy/R_e_Lens_arc_sec.csv', R_e_Lens_arc_sec, delimiter=',')

np.savetxt('SampleResult_Galaxy/Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images.csv', Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_26_and_4_images, delimiter=',')

np.savetxt('SampleResult_Galaxy/JWST_Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_28p8_and_4_images.csv', JWST_Index_in_accepted_Selected_size_less_theta_E_and_seeing_and_mag_less_28p8_and_4_images, delimiter=',')