#!/usr/bin/python
"""
Author: Nathan Totorica
Date: 4/15/2021

References:
[1] E. Medina, A. Ramirez, J. Morales and K. Sheshyekani, "Passivity Enforcement of FDNEs via Perturbation of Singularity Test Matrix," in IEEE Transactions on Power Delivery,
    vol. 35, no. 4, pp. 1648-1655, Aug. 2020, doi: 10.1109/TPWRD.2019.2949216.

[2] Houle, Jennifer, GitHub. May 10, 2020. Accessed on: February 3, 2021, [Online]. Available: https://github.com/jenniferEhoule/circuit_synthesis

"""
import numpy as np
import sys
import scipy.linalg as LA
import argparse
import math

from matplotlib import pyplot as plt
from intercheig import intercheig
from rot import rot
from pr2ss import pr2ss
from pass_check import separate_real_imag_in_state_space 
from numpy.linalg import multi_dot

class SMP(object):

    def __init__(self, plot):
        self.plot = plot


    def SMP_driver(self, SER, Niter, s_pass):
        ''' 
            Implementation of algorithm described in [1]
            
            :param SER: State space model (usually an input from running VFdriver.py)
                SER dictionary consists of ABCDE matrices
            :param Niter: Number of perturbation iterations to attempt
            :param s_pass: Frequency range for passivity enforcement
            
            :return perturbued SER dictionary
        '''
        print("\n\nStarting SMP Perturbation")
        
        # Assign SER inputs
        self.SER = pr2ss(SER)
        
        self.A = A = np.array(self.SER['A'].copy())
        self.B = B = np.array(self.SER['B'].copy())
        self.C = C = np.array(self.SER['C'].copy())
        self.D = D = np.array(self.SER['D'].copy())
        self.E = E = np.array(self.SER['E'].copy())
        self.s_pass = s_pass
        
        
        for n in range(Niter):
            print("Round {}".format(n))      
            C = self.C.copy()
            viol_idx, viol_w, delta_w = self.pass_check_Y(A, B, C, D) # Assess passivity by looking at postive real eigenvalues. 
            
            if(self.plot): # Plot eigenvalue
                self.plot_w(n_violations = len(viol_idx), iteration = n)
            
            if not np.any(viol_idx): # If singularity matrix is passive, then break from loop
                if(n > 0):
                    print("Passivity enforced in {} iteration(s)\n\n".format(n))
                else:
                    print("No cross over frequencies identified - Singular Matrix is Passive!\n\n")
                break
            
            # Calc and return array of delta w's (n,)
            for n in range(viol_idx.shape[0]):
                j = viol_idx[n]
                delta_lambda = np.zeros((self.C.shape[1]), dtype=float)
                delta_lambda[j] = 2*viol_w[n]*delta_w[n] + (delta_w[n]**2)        # i.e. Δλj=2ωjΔωj+Δωj^2 eqt. (11) from [1]
                self.prelim_perturbation(delta_lambda)                              # Preliminary perturbation step
                self.refine_perturbation(delta_lambda, j)                           # Refinement perturbation step

        # Return SER matrix
        self.SER['C'] = self.C
        return self.SER

    def pass_check_Y(self, A, B, C, D):
        """
        Calculate Singularity Matrix and assess passivity - based on pass_check.py in [2] from [1]
        
        
        A violation extending to infinity is denoted by s_end=j*1e16
        
        :param A: poles in a vector ((ports x num_poles) x (ports x num_poles))
        :param B: 
        :param C: residues in a 3D matrix (ports x (ports x num_poles))
        :param D: D values in a 2D matrix (ports x ports)
        
        :return: 
            delta_w : List of delta's for each cross over frequency
            viol_idx : List of indexs for positively real eigenvalues
            viol_wval : List of crossover frequencies
        """
        N = A.shape[0]
        Nc = D.shape[0]
        E = np.zeros((Nc, Nc))
        
        A = self.chomp(A) # Get rid of tiny imaginary numbers
        
        # Check if D is singular
        if np.sum(LA.eig(D) == 0):
            Ahat = LA.solve(A, np.eye(N))
            Bhat = - Ahat @ B
            Chat = C * Ahat
            Dhat = D - C * Ahat * B
            A = Ahat.copy()
            B = Bhat.copy()
            C = Chat.copy()
            self.D_inv = D = Dhat.copy()
        else:
            self.D_inv = np.linalg.matrix_power(D, -1)
        
        # Calculate S matrix using ABCD inputs
        S1 = A @ (B @ (np.linalg.matrix_power(D, -1)) @ C - A) # Eq. (23a) in [7]
        
        wS1, wS2 = LA.eig(S1)
        
        # Save eigenvalues and eigenvectors
        self.eig = wS1.copy() 
        self.vec = wS2.copy() 
        
        
        wS1 = self.chomp(np.sqrt(wS1))

        
        if np.sum(LA.eig(D) == 0) > 0:
            ws1 = 1 / wS1
        
        ind = np.nonzero(np.abs(wS1.imag) < 1e-6)
        viol_idx = np.array(ind[0])[:]
        
        wS1_copy = wS1.copy()
        wS1 = (wS1[ind]).real
        
        
        sing_w = np.sort(wS1)
        keep = np.nonzero(sing_w < (self.s_pass[-1] / 1j).real)
        sing_w = sing_w[keep]
        
        viol_wval = sing_w.copy()
        
        # Reorder violation index to match sorted w values
        idx = viol_idx.copy()
        for j in range(viol_idx.shape[0]):
            eig  = wS1_copy[viol_idx[j]].real
            try:
                index = np.where(viol_wval == eig)[0][0]
                idx[index] = viol_idx[j]
            except(IndexError):
                idx = np.delete(idx, np.where(idx == viol_idx[j]))
        viol_idx = idx

        if sing_w.shape[0] == 0:
            intervals = viol_idx = viol_wval = delta_w = np.array([])
            return viol_idx, viol_wval, delta_w
        
        # Establishing frequency list at midpoint of all bands defined by sing_w
        midw = np.zeros((1 + sing_w.shape[0], 1), dtype=float)
        midw[0] = sing_w[0] / 2
        midw[-1] = 2 * sing_w[-1]
        for k in range(sing_w.shape[0] - 1):
            midw[k + 1] = (sing_w[k] + sing_w[k + 1]) / 2
        
        EE = np.zeros((Nc, midw.shape[0]), dtype=complex)
        viol = np.zeros(midw.shape[0], dtype=float)
        
        # Checking passivity at all midpoints
        for k in range(midw.shape[0]):
            sk = 1j * midw[k]
            G = (self.fitcalcABCDE(sk[0], np.diag(A), B, C, D, E)).real
            EE[:, k], EE_temp = LA.eig(G)
            
            if np.any(EE[:, k] < 0, axis=0):
                viol[k] = 1
            else:
                viol[k] = 0
        
        # Establishing intervals for passivity violations:
        intervals = np.empty((2,0))
        delta_w = []
        for k in range(midw.shape[0]):
            if viol[k] == 1:
                if k == 0:
                    intervals = (np.vstack((0, sing_w[0])))
                    delta_w.append(-1*sing_w[0])

                ## Infinity case
                elif k == midw.shape[0] - 1:
                    intervals = np.hstack((intervals, (np.vstack((sing_w[k - 1], (self.s_pass[-1] / 1j).real)))))
                    delta_w.append((2* sing_w[k - 1] - sing_w[k - 1])/2)

                    
                ## Overlap case
                elif viol[k] == viol[k-1]: # If two non passive regions in a row, then there is proabably overlap for different ports
                    intervals = (np.vstack((0, sing_w[k])))
                    delta_w.append(-1*sing_w[k])
                else:
                    intervals = np.hstack((intervals, (np.vstack((sing_w[k - 1], sing_w[k])))))
                    delta_w.append((sing_w[k] - sing_w[k - 1])/2)
                    delta_w.append(-1*(sing_w[k] - sing_w[k - 1])/2)
                    
        return viol_idx, viol_wval, delta_w
    
    def prelim_perturbation(self, delta_lambda):
        """
            Perform preliminary perturbation process as described in [1]:
            
            Eigenvalue decomposition
            ΔS=PΔΛQT          (1)
               ΔS = ΔλjpjqTj  (jth position)
            where,
               Δλj=2ωjΔωj+Δω2j
            ΔS = ABD−1ΔC.     (2)
            Sub (2) into (1) and solve for ΔC via least-squares
        """
        P = self.vec.copy()                #size [(Nc*N),(Nc*N)] 
        Q = LA.inv(self.vec.copy())        #size [(Nc*N),(Nc*N)]
        
        # ΔS = PΔλQ
        delta_s = P @ np.diagflat(delta_lambda) @ Q #Eqt. (22) from [1] -- size [(Nc*N),(Nc*N)]
        
        # Least-Squares to solve ΔS =ABD−1ΔC for ΔC
        delta_c = LA.lstsq(self.A @ self.B @ self.D_inv, delta_s)[0]  # Eqt. (16) from [1] --  !! small values !!

        # Guarantee symmetry by averaging-out Δcijm = Δcjim
        self.delta_c_sym = self.force_symmetry(delta_c) # size [Nc,(Nc*N)]
        
        # ΔSp = ABD^−1ΔC
        self.delta_s_prelim = self.A @ self.B @ self.D_inv @ self.delta_c_sym # Eqt. (25) from [1] -- size [(Nc*N),(Nc*N)]

    def refine_perturbation(self, delta_lambda, j):
        '''
            Perform refinement perturbation to increase convergence, as described in [1]
            
        '''
        P = self.vec.copy()
        Q = LA.inv(self.vec.copy())
        
        epsilon = delta_lambda[j] / (Q[j,:] @ self.delta_s_prelim @ P[:,j])# Eqt. (28) from [1]

        # ΔC array - Added to original C matrix after each cross over frequency has been calculated
        delta_c_refine = epsilon*(self.delta_c_sym) # Eqt. (29) from [1]
        
        self.C += delta_c_refine     # Eqt. (15) from [1] -- Get ΔC = ΔCs + C
    
    def fitcalcABCDE(self, sk, A, B, C, D, E):
        """
        Calculate Yfit from the state space model as in Eq. (5) in [7]
        :param sk: frequency
        :param A: poles
        :param B: vector of 1's
        :param C: residues
        :param D: D values
        :param E: E values
        :return: Yfit calculation using the state space model
        """
        sk = sk
        Nc = D.shape[0]

        dum = np.tile(1 / (sk - A), (Nc, 1)).swapaxes(0, 1)
        C = C * dum.T
        Yfit = (C @ B) + D + (sk * E)
        return Yfit
    
    def plot_w(self, n_violations, iteration):
        """
            Create scatter plot of real and imaginary values of w (square root of eig(S)
        """
        w = np.sqrt(self.eig)
        w.flatten()
        x = [value.real for value in w]
        y = [value.imag for value in w]
        plt.scatter(x,y, label="Round %d"%iteration)
        plt.xlim(0, (self.s_pass[-1] / 1j).real)
        plt.ylim(-1*(self.s_pass[-1] / 1j).real, (self.s_pass[-1] / 1j).real)
        plt.xlabel('Real')
        plt.ylabel('Imag')
        plt.title("Square root eigenvalues (w) of S")
        plt.legend()
        
        if(n_violations == 0):
            plt.show()
        
    def force_symmetry(self, matrix):
        '''
            Force symmetry so that Δcijm = Δcjim
                where:
                    i = Nc
                    j = Nc
                    m = N
            
            :param matrix
                    size [Nc, (Nc*N)]
            :return sym_matrix
                    symmetric [Nc, (Nc*N)]
        '''
        x,y = np.shape(matrix)
        num_poles = int(y/x)
        sym_matrix = np.zeros(np.shape(matrix), dtype=complex)
        
        regroup = np.zeros((x,x,num_poles), dtype=complex)
        for i in range(x):
            offset = 0
            for j in range(x):
                for m in range(num_poles):
                    regroup[i][j][m] = matrix[i][offset+m]
                offset+=num_poles
        
        for i in range(x):
            offset = 0
            for j in range(x):
                for m in range(len(regroup[i][j])):
                    if(i != j):
                        avg = (regroup[i][j][m] + regroup[j][i][m]) / 2
                        regroup[i][j][m] = avg
                        regroup[j][i][m] = avg
                    sym_matrix[i][m+offset] = regroup[i][j][m]
                offset+=num_poles
        return sym_matrix
    
    def chomp(self, eig):
        """
            Remove any noisy close to zero imaginary values
        """
        check = np.isclose(eig.imag, 0, atol=1e-1)      # Remove noisy imaginary values < 0.1
        eig[check] = eig[check].real
        return eig
        
if __name__ == '__main__':
    pass
