import numpy as np
import matplotlib.pyplot as pp
from astropy.cosmology import Planck18
import astropy.units as u
import astropy.cosmology.units as cu #//FIXME
from astropy.constants import c
from MySample import sample
import random
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import special as sp
from scipy.optimize import fsolve, root
from tqdm import tqdm
import time
#================================
#for lenstronomy
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
# import the lens model class 
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions
# import the lens equation solver class (finding image plane positions of a source position)
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
# import lens model solver with 4 image positions constrains
from lenstronomy.LensModel.Solver.solver4point import Solver4Point
import lenstronomy.Util.param_util as param_util
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
#=================================
from Sersic_profile import kappa_s4

t0 = time.time()
#=============================================
#找一下所有的事例里有哪几个是被强透镜了的

SampleParam = np.loadtxt('./SampleResult/SampleParameter.csv', delimiter=',')

#Redshift, m1, m2, spin1, spin2, inclination, polarization, ra, dec

tau = []
acceptIndex = []
randomset = []
year = 10 
for i in range(len(SampleParam[0])):
    #按照徐菲的文章，把SIS的多像optical depth转换成了SIE的 
    tau.append((1 - np.exp(-4.17*10**(-6)*(Planck18.comoving_distance(SampleParam[0][i]).to(u.Gpc).value)**3)) * 2.69)

np.random.seed(1000000)

for year_i in range(year):
   
    
    for i in range(len(tau)):
        tmp = np.random.random()
        if tau[i] > tmp and i not in acceptIndex: #按道理说这里应该是没有i not in acceptIndex吧
            acceptIndex.append(i)
            randomset.append(tmp)
        else:
            pass
    

np.savetxt('./SampleResult/AcceptLensIndex.csv', acceptIndex, delimiter=',')
np.savetxt('./SampleResult/tau.csv', tau, delimiter=',')
np.savetxt('./SampleResult/AcceptRandomValue.csv', randomset, delimiter=',')






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
    kappa_s = []
    for i_tmp in range(len(x_image)):
        kappa_s.append(kappa_s4(x_image[i_tmp], y_image[i_tmp], q, theta_E, theta_E / 0.576, sigma_v))
    # print('kappa: ', kappa)
    # print('gamma: ', gamma)
    # print('time delays: ', t_days)
    return x_image, y_image, mag_, t_days, kappa, gamma, kappa_s
    



Lens_z_set = []
Q_set = []
Sigma_v_set = []
y1_set = []
y2_set = []
mag_set = []
t_days_set = []
kappa_set = []
gamma_set = []
kappa_s_set = []
x_image_set = []
y_image_set = []

for i in tqdm(range(len(acceptIndex))):
    z_s = SampleParam[0][int(acceptIndex[i])]
    
    comoving_dis_zs = Planck18.comoving_distance(z_s).value
    lens_z = GetLens_z(comoving_dis_zs)
    sigma_v, q = Get_axis_ratio_and_sigma_v()
    y1 = Parameter_y1(q)
    y2 = Parameter_y2(q)
    
    x_image, y_image, mag_, t_days, kappa, gamma, kappa_s = Solve_Lens_equation_lenstronomy(y1, y2, q, z_s, lens_z, sigma_v)
    while len(x_image) == 0 or len(x_image) == 1 or len(x_image) == 3:
        lens_z = GetLens_z(comoving_dis_zs)
        sigma_v, q = Get_axis_ratio_and_sigma_v()
        y1 = Parameter_y1(q)
        y2 = Parameter_y2(q)
        x_image, y_image, mag_, t_days, kappa, gamma, kappa_s = Solve_Lens_equation_lenstronomy(y1, y2, q, z_s, lens_z, sigma_v)
    
    
    Lens_z_set.append(lens_z)
    Q_set.append(q)
    Sigma_v_set.append(sigma_v)
    y1_set.append(y1)
    y2_set.append(y2)
    
    mag_set.append(list(mag_))
    t_days_set.append(list(t_days))
    kappa_set.append(list(kappa))
    gamma_set.append(list(gamma))
    kappa_s_set.append(list(kappa_s))
    x_image_set.append(list(x_image))
    y_image_set.append(list(y_image))
    

image_num = []
for i in range(len(mag_set)):
    image_num.append(len(mag_set[i]))
    
mag_set = [y for x in mag_set for y in x]
t_days_set = [y for x in t_days_set for y in x]
kappa_set = [y for x in kappa_set for y in x]
gamma_set = [y for x in gamma_set for y in x]
kappa_s_set = [y for x in kappa_s_set for y in x]
x_image_set = [y for x in x_image_set for y in x]
y_image_set = [y for x in y_image_set for y in x]
z_s_set = []
for i in tqdm(range(len(acceptIndex))):
    z_s_set.append(SampleParam[0][int(acceptIndex[i])])


np.savetxt('./SampleResult/magnification.csv', mag_set, delimiter=',')
np.savetxt('./SampleResult/timedelayDays.csv', t_days_set, delimiter=',')
np.savetxt('./SampleResult/kappa.csv', kappa_set, delimiter=',')
np.savetxt('./SampleResult/gamma.csv', gamma_set, delimiter=',')
np.savetxt('./SampleResult/kappa_s.csv', kappa_s_set, delimiter=',')
np.savetxt('./SampleResult/imagenumber.csv', image_num, delimiter=',')
np.savetxt('./SampleResult/lens_z.csv', Lens_z_set, delimiter=',')
np.savetxt('./SampleResult/axis_ratio.csv', Q_set, delimiter=',')
np.savetxt('./SampleResult/sigma_v.csv', Sigma_v_set, delimiter=',')
np.savetxt('./SampleResult/y1.csv', y1_set, delimiter=',')
np.savetxt('./SampleResult/y2.csv', y2_set, delimiter=',')
np.savetxt('./SampleResult/x1_image.csv', x_image_set, delimiter=',')
np.savetxt('./SampleResult/x2_image.csv', y_image_set, delimiter=',')
np.savetxt('./SampleResult/z_s.csv', z_s_set, delimiter=',')


import matplotlib.pyplot as plt
mag_set = np.loadtxt('./SampleResult/magnification.csv', delimiter=',')
t_days_set = np.loadtxt('./SampleResult/timedelayDays.csv', delimiter=',')
kappa_set = np.loadtxt('./SampleResult/kappa.csv', delimiter=',')
gamma_set = np.loadtxt('./SampleResult/gamma.csv', delimiter=',')
kappa_s_set = np.loadtxt('./SampleResult/kappa_s.csv', delimiter=',')
image_num = np.loadtxt('./SampleResult/imagenumber.csv', delimiter=',')
Lens_z_set = np.loadtxt('./SampleResult/lens_z.csv', delimiter=',')
x_image_set= np.loadtxt('./SampleResult/x1_image.csv', delimiter=',')
y_image_set= np.loadtxt('./SampleResult/x2_image.csv', delimiter=',')
y1_set = np.loadtxt('./SampleResult/y1.csv', delimiter=',')
y2_set = np.loadtxt('./SampleResult/y2.csv', delimiter=',')

#PDF
fig, ax = plt.subplots(2, 3, figsize=(12, 8))

ax[0,0].hist(mag_set, density=True, bins=1000)
ax[0,0].set_xlim(-10, 10)
ax[0,0].set_xlabel('$\mu$')
ax[0,0].grid()


ax[0,1].hist(t_days_set * 24, density=True, bins=1000)
ax[0,1].set_xlim(-4500, 1500)
ax[0,1].set_xlabel('$\Delta t$ [hr]')
ax[0,1].grid()

ax[0,2].hist(kappa_set, density=True, bins=50000)
ax[0,2].set_xlim(0, 2)
ax[0,2].set_xlabel('$\kappa_\mathrm{total}$')
ax[0,2].grid()

ax[1,0].hist(kappa_s_set, density=True, bins=10000)
ax[1,0].set_xlim(0.001, 1000)
ax[1,0].set_xlabel('$\kappa_*$')
ax[1,0].grid()
ax[1,0].semilogx()


ax[1,1].hist(kappa_s_set / kappa_set, density=True, bins=1000)
ax[1,1].set_xlim(0, 1.5)
ax[1,1].set_xlabel('$f_*$')
ax[1,1].grid()

ax[1,2].hist(Lens_z_set, density=True, bins=1000)
ax[1,2].set_xlabel('$z_l$')
ax[1,2].grid()


plt.savefig('./IntermediaPlot/Lens_sample_result.png', dpi=450)
plt.close()


plt.scatter(np.abs(mag_set), kappa_s_set / kappa_set)
plt.plot([0,1000],[0.5, 0.5], 'k')
plt.plot([1, 1], [0, 10], 'k')
plt.xlim(10**(-3), 100)
plt.ylim(10**(-1), 2)
plt.xlabel('$\mu$')
plt.ylabel('$f_*$')
plt.loglog()
plt.grid()
plt.savefig('./IntermediaPlot/mag_vs_f_s.png.png', dpi=450)
plt.close()


plt.scatter(np.abs(mag_set), kappa_s_set)
# plt.plot([0,1000],[0.5, 0.5], 'k')
plt.plot([1, 1], [0, 1000], 'k')
plt.loglog()
plt.xlim(10**(-4), 1000)
plt.ylim(10**(-2), 150)
plt.xlabel('$\mu$')
plt.ylabel('$\kappa_*$')
plt.grid()
plt.savefig('./IntermediaPlot/mag_vs_kappa_s.png.png', dpi=450)
plt.close()


plt.scatter(x_image_set, y_image_set)
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.grid()
plt.savefig('./IntermediaPlot/Lens_image_position.png', dpi=450)
plt.close()

plt.scatter(y1_set, y2_set)
plt.xlabel('$y_1$')
plt.ylabel('$y_2$')
plt.grid()
plt.savefig('./IntermediaPlot/Lens_source_position.png', dpi=450)
plt.close()
print('Run time: ', time.time() - t0)


