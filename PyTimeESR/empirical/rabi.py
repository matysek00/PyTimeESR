# -*- coding: utf-8 -*-

import numpy as np

class Rabi_2Level():
    """ Implementing Equatons from [1]. I didn't get them to work though.

    [1], Reina-Galvez, M Nachtigall, N Lorente, J Martinek, C Wolf
            arXiv preprint arXiv:2503.24046, 2025•arxiv.org
    """
    def __init__(self, ham_dict, floquet_dict):
        self.ham_dict = ham_dict
        self.floquet_dict = floquet_dict
       
    def T2(self):
        """Eq 26 
            Reina-Galvez, M Nachtigall, N Lorente, J Martinek, C Wolf
            arXiv preprint arXiv:2503.24046, 2025•arxiv.org
        """
        
        gamma0 = self.floquet_dict['gamma0']
        A = self.floquet_dict['A']
        efermi = self.floquet_dict['bias']
        epsilon = self.ham_dict['eps_QD']
        U = self.ham_dict['Hubbard']
        T2inv = 0
        
        for i in range(2):
            T2inv += gamma0[i] * np.pi * (1+A[i]**2/2) * (1 - np.sign(efermi[i]-epsilon)/2 + np.sign(efermi[i]-epsilon-U)/2)
        
        return 1/T2inv

    def QualitySeq(self, Electrode):
        """Eq 28 
            Reina-Galvez, M Nachtigall, N Lorente, J Martinek, C Wolf
            arXiv preprint arXiv:2503.24046, 2025•arxiv.org
        """
        print('WARNING: This function assumes that gammaC >> hbar ')

        efermi = self.floquet_dict['bias']
        vdc = efermi[Electrode] - efermi[1-Electrode]
        
        A = self.floquet_dict['A'][Electrode]
        P = self.floquet_dict['Spin_polarization'][Electrode]
        
        sin = 1. # TODO: Fix this.
        print(A,P,efermi, vdc, sin)
        Q = -np.sign(vdc)*A*P/(2*(1+A**2/2))
        return Q 
        
    def QualityCot(self, Electrode):
        """Eq 29 
            Reina-Galvez, M Nachtigall, N Lorente, J Martinek, C Wolf
            arXiv preprint arXiv:2503.24046, 2025•arxiv.org
        """
        print('WARNING: This function assumes that gammaC >> hbar ')

        efermi = self.floquet_dict['bias']
        vdc = efermi[Electrode] - efermi[1-Electrode]
        sin = 1. # TODO: Fix this.
        A = self.floquet_dict['A'][Electrode]
        P = self.floquet_dict['Spin_polarization'][Electrode]
        g0 = self.floquet_dict['gamma0']
        epsilon = self.ham_dict['eps_QD']
        U = self.ham_dict['Hubbard']

        sin = 1. # TODO: Fix this.

        log = np.log((vdc+epsilon+U)**2/(vdc+epsilon)**2)
        nomin = -g0[Electrode]*A*P*sin*log

        denom1 = g0[Electrode]   * (2 - np.sign(vdc-epsilon) + np.sign(vdc-epsilon-U))
        denom2 = g0[1-Electrode] * (2 - np.sign(epsilon) + np.sign(epsilon+U))
        denom  = np.pi*((1+A**2/2)*denom1 + denom2)
        Q = nomin/denom
        
        return Q 
        