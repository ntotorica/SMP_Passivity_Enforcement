"""
Author: Nathan Totorica
Date: 4/29/2021

References:
[1] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions on Power Delivery,
    vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

[2] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: February 3, 2021, [Online]. Available: https://github.com/jenniferEhoule/circuit_synthesis

"""
import numpy as np
import warnings
warnings.filterwarnings(action="error", category=np.ComplexWarning)
import time
import skrf as rf

from VFdriver import VFdriver
from RPdriver import RPdriver
from plots import plot_figure_11
from create_netlist import create_netlist_file
from SMP import SMP
from eig_plot import plot
from eig_plot import plot_s_parameters
from convert_y_to_y_prime import convert_y_to_y_prime_matrix_values_3d

##########################
# Import Y-parameters
##########################
#SP = rf.Network('s_param_16.s16p')
SP = rf.Network('s_param_4.s4p')
#SP = rf.Network('s_param_8.s8p')
#SP = rf.Network('s_param_4_infini.s4p')
#SP = rf.Network('s_param_test.s2p') # Inifinity case
#SP = rf.Network('s_param_3.s3p')
s = 2*np.pi*1j*SP.frequency.f
bigY = SP.y
bigS = SP.s

################################################
# Need to reshape - skrf y_param -> (f,N,N)
################################################
bigY_new = np.zeros((bigY.shape[1], bigY.shape[2], s.shape[0]),dtype=complex)

for freq in range(len(s)):
    bigY_new[:,:,freq] = bigY[freq,:,:]

bigY=bigY_new

################################################
# Vector Fitting
################################################
bigY_prime = convert_y_to_y_prime_matrix_values_3d(bigY)

# Pole-Residue Fitting
vf_driver = VFdriver(N=30,
                     poletype='linlogcmplx',
                     weightparam='common_sqrt',
                     Niter1=7,
                     Niter2=4,
                     asymp='D',
                     logx=False,
                     plot=False
                     )
poles=None
SER, rmserr, bigYfit = vf_driver.vfdriver(bigY_prime, s, poles)

prev_SER = SER.copy()

ylim1 = None
ylim2 = None
plot(plot_name='Original', s_pass=s.T, ylim=np.array((ylim1, ylim2)), labels=['Original'], SER1=prev_SER)
################################################
# Passivity Enforcement
################################################
start_time = time.process_time()

# Initalize SMP object
smp = SMP(plot=False)
SMP_SER = smp.SMP_driver(SER, Niter=10, s_pass=s.T)
print("{} s".format(time.process_time() - start_time))

plot(plot_name='SMP Enforecment', s_pass=s.T, ylim=np.array((ylim1, ylim2)), labels=['SMP Perturbed'], SER1=SMP_SER)
plot(plot_name='Orig Vs. SMP', s_pass=s.T, ylim=np.array((ylim1, ylim2)), labels=['Original', 'SMP Perturbed'], SER1=prev_SER, SER2=SMP_SER)

################################################
# Generate Netlist
################################################
poles = SMP_SER['poles']
residues = SMP_SER['C']
Ns = poles.shape[0]
Nc = int(residues.shape[1] / Ns)
poles = poles.reshape((1, -1))
residues = residues.reshape((Nc ** 2, Ns))
create_netlist_file(poles, residues)
