#!/usr/bin/python

from math import pi, sqrt
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from intercheig import intercheig
from rot import rot
import skrf as rf

'''
References:
[1] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions on Power Delivery,
    vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

[2] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: February 3, 2021, [Online]. Available: https://github.com/jenniferEhoule/circuit_synthesis


    Methods calculate_y_from_SER() and calculate_eigenvalues_for_EE() implemented based on [2]
'''


def plot(plot_name, s_pass, ylim, SER1, labels, SER2=None):
    SER_list = [SER1]
    if(SER2):
        SER_list.append(SER2)
    s = s_pass
    graphs = []
    
    for graph in SER_list:
        A = graph['A']
        B = graph['B']
        C = graph['C']
        D = graph['D']
        E = graph['E']
        Nc = np.shape(D)[0]
        oldT0 = np.array([])
        I = np.ones((A.shape[0], 1))
        EE0 = np.zeros((Nc, s_pass.shape[0]), dtype=complex)
        for k in range(s_pass.shape[0]):
            Y = calculate_y_from_SER(I, graph, k, s_pass) #
            # Fills in the eigenvalues into EE0
            EE0, oldT0 = calculate_eigenvalues_for_EE(EE0, Nc, Y, k, oldT0)
            
        graphs.append(EE0)
    plot_eigenvalues_of_gs(plot_name, s_pass, graphs, ylim, labels)
   


def calculate_y_from_SER(I, SER1, k, s_pass):
    """ Implements Eq. (15) in [6] """
    Y = SER1['C'] @ np.diagflat((s_pass[k] * I - np.diag(SER1['A']).reshape(-1, 1)) ** (-1)) @ SER1['B'] + SER1['D'] + s_pass[k] * SER1['E']
    return Y
    
def calculate_eigenvalues_for_EE(EE0, Nc, Y, k, oldT0):
    """
    Calculate the eiginevalues of G(s) and fill into EE0 matrix for the given value of k
    :param EE0: Matrix of the eigenvalues of G(s)
    :param Nc: Number of ports on which data is being fit
    :param Y: The result of Eq. (15) in [6]
    :param k: Index for the iteration number
    :param oldT0: Eigenvectors from pervious iteration
    :return: EE0 with additional data for given k; oldT0 with eigenvectors from the current index k
    """
    G = Y.real
    D, T0 = LA.eig(G) # Calculate eigenvalues / eigenvectors to evaluate
    T0 = rot(T0.astype(complex))  # Minimizing phase angle of eigenvectors in least squares sense
    T0, D = intercheig(T0, oldT0, np.diag(D).copy(), Nc, k) # Rearrange the eigenvalues / vectors to smooth them out over frequency
    oldT0 = T0
    EE0[:, k] = np.diag(D)
    return EE0, oldT0

def plot_eigenvalues_of_gs(plot_name, s_pass, graphs, ylim, labels):
    """
    Plots eigenvalues of G(s)
    :param s_pass: Frequencies being sampled
    :param EE0: Eigenvalues of the current model
    :param xlimflag: Indicates xlim should be used
    :param ylimflag: Indicates ylim should be used
    :return:
    """
    plt.rcParams['font.size'] = 12
    plt.rcParams['grid.color'] = 'gray'
    plt.rcParams['grid.linestyle'] = 'dotted'
    
    fig = plt.figure(plot_name, figsize=(8, 7))
    ax = fig.add_subplot(1, 1, 1)
    freq = (s_pass / (2 * np.pi * 1j)).real
    
    ax.plot(freq, (graphs[0].T).real, color='b', linewidth=1, label=labels[0])
    
    if len(graphs) > 1:
        ax.plot(freq, (graphs[1].T).real, color='r', linewidth=1, label=labels[1])
        
        rms_error(graphs[1].T, graphs[0].T, freq)
            
    plt.axhline(y=0.0, color='g', linestyle='-', label="zero boundary")
    if(np.any(ylim != None)):
        plt.ylim(ylim[0], ylim[1])
        plt.xlim(freq[0], freq[-1])
    
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Eigenvalues of G')
    plt.legend()
    plt.title("Eigenvalues of G(s)")
    plt.tight_layout()
    plt.savefig('eigenvalues_of_Gs')
    plt.show()

    return
    
def rms_error(prev_EE0, EE0,freq):
     diff = EE0 - prev_EE0
     print("RMSE 1: {}".format(np.sqrt(np.sum(np.abs(diff ** 2))) / np.sqrt(freq.shape[0])))
     
     mse = (np.square(np.subtract(prev_EE0.real,EE0.real))).mean(axis=None)
     print("MSE {}".format(mse))
     print("RMSE {}".format(np.sqrt(mse)))
     return
     
def plot_s_parameters(s1, s2, Nc, freq):
    Ns = freq.shape[0]
    plt.rcParams['font.size'] = 12
    plt.rcParams['grid.color'] = 'gray'
    plt.rcParams['grid.linestyle'] = 'dotted'
            
    fig = plt.figure(11, figsize=(8, 7))
    ax = fig.add_subplot(1, 1, 1)
    for row in range(Nc):
        for col in range(Nc):
            dum1 = s1[row, col, :]
            dum2 = s2[row, col, :]
            ax.plot(freq, np.abs(dum1), color='b', linewidth=1, label='Original')
            ax.plot(freq, np.abs(dum2), color='r', linewidth=1, label='Perturbed')
            diff = dum2 - dum1
            ax.plot(freq, np.abs(diff), color='g', linewidth=1, label='Difference')

            plt.xlabel('Frequency [Hz]')
            plt.ylabel('S-Parameter')
            plt.title("Model Comparison")
            plt.legend()
            plt.yscale('log')
            print(f"RMS error S{row+1}{col+1}: {np.sqrt(np.sum(np.abs(diff ** 2))) / np.sqrt(Ns)}")


    plt.show()
    
   