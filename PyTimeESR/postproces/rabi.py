import numpy as np


def Rabi_oscilation(fn_rates: str, n_max: int, hbar: float= 2.418883402334483e-8,):
    G = np.loadtxt(fn_rates, skiprows=1)
    Nlevels = 4

    tauinv = np.zeros(2*n_max+1)
    omega_SST = np.zeros(2*n_max+1)
    omega_FLT = np.zeros(2*n_max+1)

    for n in range(-n_max, n_max+1):
        nidx = n+n_max

        Gn = G[np.where(G[:, 4]==n)]
        RGn = (Gn[:,5]+ Gn[:,7]).reshape((Nlevels, Nlevels, Nlevels, Nlevels))
        IGn = (Gn[:,6]+ Gn[:,8]).reshape((Nlevels, Nlevels, Nlevels, Nlevels))
        
        for v in range(2,3):
            
            omega_SST[nidx] += RGn[0,v,v,1]
            omega_FLT[nidx] += IGn[0,v,v,1]

            for j in range(0,1):
                tauinv[nidx] += RGn[j,v,v,j]
    
    tauinv *= 1/hbar
    omega_FLT *= 1/hbar
    omega_SST *= 1/hbar

    return tauinv, omega_SST, omega_FLT

