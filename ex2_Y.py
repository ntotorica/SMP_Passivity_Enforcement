""" ex2_Y.py

This is based on ex2_Y.m from the Matrix Fitting Toolbox [1].

Author: Jennifer Houle
Date: 3/18/2020

[1] B. Gustavsen, Matrix Fitting Toolbox, The Vector Fitting Website. March 20, 2013. Accessed on:
    Feb. 25, 2020. [Online]. Available:
    https://www.sintef.no/projectweb/vectorfitting/downloads/matrix-fitting-toolbox/.

"""


import numpy as np
import warnings
warnings.filterwarnings(action="error", category=np.ComplexWarning)
import time

from VFdriver import VFdriver
from RPdriver import RPdriver
from plots import plot_figure_11
from create_netlist import create_netlist_file
from SMP import SMP
from eig_plot import plot


bigY = np.load('ex2_Y_bigY.npy')
s = np.load('ex2_Y_s.npy')
s = np.ravel(s)

# Pole-Residue Fitting
vf_driver = VFdriver(N=50,
                     poletype='linlogcmplx',
                     weightparam='common_sqrt',
                     Niter1=7,
                     Niter2=4,
                     asymp='D',
                     logx=False,
                     plot=False
                     )
poles=None
SER, rmserr, bigYfit = vf_driver.vfdriver(bigY, s, poles)


prev_SER = SER.copy()
ylim1 = -2e-3
ylim2 = -1*ylim1

plot(plot_name='Original', s_pass=2*np.pi*1j*np.linspace(0, 2e5, 1001).T, ylim=np.array((ylim1, ylim2)), labels=['Original'], SER1=prev_SER)
################################################
# Initalize SMP object
################################################
start_time = time.process_time()
smp = SMP(plot=True)
SMP_SER = smp.SMP_driver(SER, Niter=10, s_pass=2*np.pi*1j*np.linspace(0, 2e5, 1001).T)
print("{} s".format(time.process_time() - start_time))


plot(plot_name='SMP Enforecment', s_pass=2*np.pi*1j*np.linspace(0, 2e5, 1001).T, ylim=np.array((ylim1, ylim2)), labels=['SMP Perturbed'], SER1=SMP_SER)
plot(plot_name='Orig Vs. SMP', s_pass=2*np.pi*1j*np.linspace(0, 2e5, 1001).T, ylim=np.array((ylim1, ylim2)), labels=['Original', 'SMP Perturbed'], SER1=prev_SER, SER2=SMP_SER)
################################################


# plot_figure_11(s, bigY, bigYfit, SER)
# # Passivity Enforcement
# start_time = time.process_time()
# rp_driver = RPdriver(parametertype='y',
                     # s_pass=2*np.pi*1j*np.linspace(0, 2e5, 1001).T,
                     # ylim=np.array((-2e-3, 2e-3)))
# RP_SER, bigYfit_passive, opts3 = rp_driver.rpdriver(SER, s)
# print("{} s".format(time.process_time() - start_time))
# plot_figure_11(s, bigY, bigYfit_passive, SER)
# plot(plot_name='Orig Vs. SMP', s_pass=2*np.pi*1j*np.linspace(0, 2e5, 1001).T, ylim=np.array((ylim1, ylim2)), labels=['Original', 'SMP Perturbed'], SER1=prev_SER, SER2=RP_SER)


poles = SMP_SER['poles']
residues = SMP_SER['C']
Ns = poles.shape[0]
Nc = int(residues.shape[1] / Ns)
poles = poles.reshape((1, -1))
residues = residues.reshape((Nc ** 2, Ns))
create_netlist_file(poles, residues)

