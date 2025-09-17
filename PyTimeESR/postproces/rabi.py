import numpy as np
import itertools

     
def expand_rates(G_red, N):
    """Turn a sparse version of Gamma into a full tensor. 
    In the sparse version only nonzero components are pressent,
    each is of the form [u,v,i,j,n, Re(GL), Im(GL), Re(GR), Im(GR)]. 
    NOTE: it is assumed that only one Fourier component is being passed. 
    In the full tensor, has N*N*N*N componens of the form [Re(GL), Im(GL), Re(GR), Im(GR)] 

    Parameters: 
        G_red (np.array, (:, 9)): 
            Gamma in reduced form
        N (int):
            Number of levels in the system
    Returns:
        G (np.array, (N, N, N, N,4)):
            Gamma as a tensor 
    """
    G = np.zeros((N, N, N, N,4))

    for comp in G_red:
        G[tuple((comp[:4]-1).astype(int))] = comp[5:]
            
    return G 

def Rabi_oscilation(fn_rates: str, n_max: int, hbar: float= 2.418883402334483e-8,):
    """Calculates the folowing 
    Omega_SST = Re Gamma(g,v,v,e)
    Omega_FLT = Im Gamma(g,v,v,e)
    tau^-1 = Re Gamma(j,v,v,j)
    where v runs over 0,2 and j runs over e,g  sum over elctrodes is implied

    Parameters:
    -----------
    fn_rates (str):
        file with stored rates
    n_max (int):
        maximum fourier component to consider
    hbar (float, optional): 
        hbar, with apropriate time energy conversion. 

    Returns:
    --------
    tauinv (np.array, (2*nmax+1)): 
        inverse decay time 
    omegaFLT (np.array, (2*nmax+1)): 
        Rabi rate in CB regime
    omegaSST (np.array, (2*nmax+1)): 
        Rabi rate in sequential regime. 
    """
    G = np.loadtxt(fn_rates, skiprows=1)
    Nlevels = 4

    tauinv = np.zeros(2*n_max+1)
    omega_SST = np.zeros(2*n_max+1)
    omega_FLT = np.zeros(2*n_max+1)

    for n in range(-n_max, n_max+1):
        nidx = n+n_max
        
        # select correct fourier components
        Gn = G[np.where(G[:, 4]==n)]
        if Gn.size == 0: 
            continue
        Gn = expand_rates(Gn, Nlevels)   

        RGn = (Gn[:,:,:,:,0]+ Gn[:,:,:,:,2])
        IGn = (Gn[:,:,:,:,1]+ Gn[:,:,:,:,3])
        
        for v in range(2,4): # v = 0,2
            
            omega_SST[nidx] += RGn[0,v,v,1]
            omega_FLT[nidx] += IGn[0,v,v,1]

            for j in range(0,2): # j  = g,e
                tauinv[nidx] += RGn[j,v,v,j]
    
    tauinv *= 1/hbar
    omega_FLT *= 1/hbar
    omega_SST *= 1/hbar

    return tauinv, omega_SST, omega_FLT

