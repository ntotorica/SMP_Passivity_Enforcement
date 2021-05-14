#!/user/bin/python
'''
 Used to compare S-parameters. Will plot magnitude differences and output RMS error
'''

import numpy as np
import skrf as rf
from eig_plot import plot
from eig_plot import plot_s_parameters
import argparse


if __name__ == "__main__":
    # User Input Arguments
    parser = argparse.ArgumentParser(description = 'Provide S-parameter files')
    parser.add_argument('-s1', help='')
    parser.add_argument('-s2', help='')
    
    
    args = parser.parse_args()
    
    ############################
    # S-parameter comparison
    ############################
    SP1 = rf.Network(args.s1)
    SP2 = rf.Network(args.s2)
    
    
    s_param1 = SP1.s
    s_param2 = SP2.s
    freq = SP2.f
    
    s_param1_new = np.zeros((s_param1.shape[1], s_param1.shape[2], freq.shape[0]),dtype=complex)
    s_param2_new = np.zeros((s_param1.shape[1], s_param1.shape[2], freq.shape[0]),dtype=complex)
    Nc = s_param2_new.shape[0]
    
    for f in range(len(freq)):
        s_param1_new[:,:,f] = s_param1[f,:,:]
        s_param2_new[:,:,f] = s_param2[f,:,:]
    plot_s_parameters(s_param1_new, s_param2_new, Nc, freq)