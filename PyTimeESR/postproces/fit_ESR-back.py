import numpy as np
from scipy.optimize import curve_fit



def create_fit_esr(x, y, Zeeman, lb=None, ub=None, maxfev=None, size = 1.,tol=1e-5, weight_scale=0):
    """Probably redundant funcition
   
    p0 = (res, gamma, Isym, Iasym, p0, p1, p2, pinv)
    
    Ips = Isym/(1 + 4 * (x - res)**2 / gamma**2)
    Ipa = 2 * Iasym * (x - res) / gamma / (1 + 4 * (x - res)**2/ gamma**2)
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    I = I0 + Ips + Ipa
    weighscale: if >0, apply a weight that penalizes points away from guessed resonance
        to put more emphasis on fitting the peak region.
    """
    ### NEED To FIT ALL Peaks at SIMULTANEOUSLY
    harmmax = Zeeman/np.min(x)
    harmmin = Zeeman/np.max(x)

    harms = np.arange(int(harmmin)+1, int(harmmax)+1, 1)

    if len(harms) < 1:
        print("No harmonics in the fitting range")
        opt = np.zeros(8)
        opt[0] = Zeeman
        return opt, np.zeros_like(opt)
    
    error_return_opt = (np.zeros(4), np.zeros((len(harms),4)))
    error_return_cov = (np.zeros(4), np.zeros((len(harms),4)))
    
    b = (y[-1]-y[0])/(x[-1]-x[0])
    a = y[0] - b*x[0]

    p0sback = [a,b,0.,0.]  # initial background guess
    lbsback = [-np.inf, -np.inf, -np.inf, -np.inf]  # background lower bound
    ubsback = [np.inf, np.inf, np.inf, np.inf]      # background upper bound
    
    lbs = lbsback.copy()
    ubs = ubsback.copy()
    p0s = p0sback.copy()

    new_harms = []
    errors = []

    # Guess initial parameters for each harmonic
    max_var = 0.
    
    for ih, harm in enumerate(harms):

        Zeeman_h = Zeeman / harm
        use = np.abs(x - Zeeman_h) < size
        var = np.max(y[use]) - np.min(y[use])
    
        max_var = max(var, max_var)
        
        oscilations = (np.sum(y[use] == np.max(y[use])) > 2
                       ) and (np.sum(y[use] == np.min(y[use])) > 2)
       
        if var/max_var < tol:
            # too small to fit
            continue
        if oscilations:
            # oscillatory data, skip
            continue

        p0 = [Zeeman_h, 0., 0., 0.]#guess_p0(x, y, Zeeman_h, a, use)
        ub = [Zeeman_h + size/2, np.inf, 1.1*var, 1.1*var]
        lb = [Zeeman_h - size/2, -np.inf, -1.1*var, -1.1*var]
        
        popt, pcov = curve_fit(ESR_fit_fun_single, x[use], y[use], p0=p0sback+p0, bounds=(lbsback+lb,ubsback+ub), maxfev=maxfev)
        fit = ESR_fit_fun_single(x[use],*popt)
 
        errors.append(get_error(x[use], y[use], size, fit, popt[4]))
        print(errors)

        p0s += popt[4:].tolist()
        lbs += lb
        ubs += ub

        new_harms.append(harm)
        error_return_opt[1][ih,:2] = p0[:2]
    
    harms = new_harms
    p0s = np.array(p0s)
    bounds = (np.array(lbs), np.array(ubs))
    bad_bounds = ~np.logical_and(p0s > bounds[0], p0s < bounds[1])

    if bad_bounds.any():
        exit

    if len(harms) < 1:
        return (p0s, []), (np.zeros_like(p0s), []
                           )
    # Initial parameter vector
    def fitting_fun(x, *params):
        I0 = background_current(x, *params[:4])
        Ipeak = 0

        for ih, harm in enumerate(harms):
            res = params[4 + ih*4 + 0]
            gamma = params[4 + ih*4 + 1]
            Isym = params[4 + ih*4 + 2]
            Iasym = params[4 + ih*4 + 3]
            Ipeak += peak_curent(x, res, gamma, Isym, Iasym)
        return I0 + Ipeak
    
    try:
        popt, pcov = curve_fit(fitting_fun, x, y, 
            p0=p0s, bounds=bounds, maxfev=maxfev)       
    except RuntimeError:
        return (p0s[:4], []), (np.zeros_like(p0s[:4]), [])
       
    back_par = popt[:4]
    peak_parmas = np.reshape(popt[4:], (len(harms), 4))
    peak_params_single = np.reshape(p0s[4:], (len(harms), 4))

    total_fit = ESR_fit_fun(x, back_par, peak_parmas)

    for i,p in enumerate(peak_parmas):
        err = get_error(x, y, size, total_fit, p[0])
        if errors[i] < err: 
            
            # the single peak fit is preferable
            peak_parmas[i] = peak_params_single[i]

        if peak_parmas[i][1] < 0:
        # negative gamma 
            peak_parmas[i][1] = -peak_parmas[i][1]
            peak_parmas[i][3] = -peak_parmas[i][3]

    pvar = np.diag(pcov)
    back_cov = pvar[:4]
    peak_cov = np.reshape(pvar[4:], (len(harms), 4))

    return (back_par, peak_parmas), (back_cov, peak_cov)

def get_error(x, y, size, total_fit, Zeeman):
    use = use = np.abs(x - Zeeman) < size
    var = np.max(y[use]) - np.min(y[use])
    err = np.sqrt(np.sum((total_fit[use] - y[use])**2))/np.sum(use)/var
    return err

def ESR_fit_fun(x, back_par, peak_parmas):
    total_fit = background_current(x, *back_par)
    for p in peak_parmas:
        total_fit += peak_curent(x, *p)
    return total_fit

def ESR_fit_fun_single(x, res, gamma, Isym, Iasym, p0, p1, p2, pinv):
    """
    ESR fitting function.
    
    Ips = Isym/(1 + (x - res)**2 / gamma**2)
    Ipa = 2 * Iasym * (x - res) / (2 * (x - res)**2 + gamma**2)
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    I = I0 + Ips + Ipa

    Parameters:
    x : array_like
        The independent variable (frequency or field).
    res : float
        Resonance position.
    gamma : float
        Lorentzian width (damping).
    Isym : float
        Symmetric peak current amplitude.
    Iasym : float
        Asymmetric peak current amplitude (Fano).
    p0, p1, p2 : float
        Coefficients for the polynomial background current.
    pinv : float
        Inverse background current coefficient.
    Returns:
    y : array_like
        The fitted values for the given x.
    """

    Background = background_current(x, p0, p1, p2, pinv) 
    peak = peak_curent(x, res, gamma, Isym, Iasym)
              
    return peak + Background

def peak_curent(x, res, gamma, Isym, Iasym):
    """
    Peak current function for ESR fitting.
    Coppied from Jose's Matlab code.

    peak_sym = peak0/(1 + (x - res)**2 / gamma**2)
    peak_asym = 2 * peak1 * (x - res) / (2 * (x - res)**2 + gamma**2)
    return peak_sym + peak_asym
    """
    
    y = 2*(x - res)/gamma 
    # don't know why the factor of 2 is needed, but it is in the original code

    nominator = Isym + Iasym * y
    denominator = y**2 + 1
    peak = nominator / denominator
    return peak

def background_current(x, p0, p1, p2, pinv,):
    """
    Polynomial background current function for ESR fitting.
    Coppied from Jose's Matlab code.

    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    """
    I0 = pinv/(x+1e-10) + p0 + p1 * x + p2 * x**2
    return I0

def guess_p0(x, y, x0, a, use):

    y = y[use]
    x = x[use]

    argmin = np.argmin(y)
    argmax = np.argmax(y)

    
    ymax= y[argmax]
    ymin= y[argmin]

    maxmin = (ymax - a)/(ymin - a)

    rat0 = 100
    if maxmin > rat0:
        Ia = 0
        Is = ymax - a
        
        half_width1 = np.argmin(ymax/2 - a - y[:argmax])
        half_width2 = np.argmin(ymax/2 - a - y[argmax:])
        gamma = (x[argmax - half_width1] - x[argmax + half_width2])/2

    elif 1/maxmin > rat0:
        Ia = 0
        Is = -ymin - a

        half_width1 = np.argmin(-ymin/2 - a - y[:argmin])
        half_width2 = np.argmin(-ymin/2 - a - y[argmin:])
        gamma = (x[argmin - half_width1] - x[argmin + half_width2])/2

    else:
        Is = -(ymax + ymin)
        Ia = (ymax- ymin)
        gamma = np.abs(x[argmax] - x[argmin])/2

    p0 = [x0, gamma, Is, Ia]
    return p0