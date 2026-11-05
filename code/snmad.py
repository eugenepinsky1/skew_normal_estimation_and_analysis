"""Quartile-MAD estimation for the skew-normal distribution SN(m, s, alpha)."""
import numpy as np
from scipy.special import owens_t, ndtr
from scipy.stats import norm, skewnorm
from scipy.optimize import least_squares, minimize

B = np.sqrt(2/np.pi)

def F0(a, al):
    """CDF of standard SN(0,1,al)."""
    return ndtr(a) - 2*owens_t(a, al)

def H(a, al):
    """E|X-a| for X ~ SN(0,1,al)."""
    c = np.sqrt(1+al**2); d = al/c
    return a*(2*F0(a, al)-1) + 4*norm.pdf(a)*ndtr(al*a) - B*d*(2*ndtr(a*c)-1)

def dH_da(a, al):
    return 2*F0(a, al)-1

def dH_dal(a, al):
    c2 = 1+al**2
    return B*(1-2*ndtr(a*np.sqrt(c2)))/c2**1.5

def sample_stats(y):
    q = np.quantile(y, [0.25, 0.5, 0.75])   # type-7 sample quartiles
    d = np.abs(y[:, None] - q[None, :]).mean(0)
    return q, d

def resid(theta, q, d):
    m, ls, al = theta
    s = np.exp(ls)
    return (s*H((q-m)/s, al) - d)/d

def jac(theta, q, d):
    m, ls, al = theta
    s = np.exp(ls); z = (q-m)/s
    J = np.empty((3, 3))
    J[:, 0] = -dH_da(z, al)
    J[:, 1] = s*(H(z, al) - z*dH_da(z, al))   # d/d log s
    J[:, 2] = s*dH_dal(z, al)
    return J/d[:, None]

def start_values(y, q):
    # quick moment-type start: sign of skew from quartile skewness (Bowley)
    bow = (q[2]+q[0]-2*q[1])/max(q[2]-q[0], 1e-12)
    al0 = np.clip(10*bow, -8, 8)
    return al0

def qmad_fit(y, alpha_max=50.0, starts=(None,)):
    q, d = sample_stats(y)
    best = None
    al_guesses = [start_values(y, q), 1.0, -1.0, 4.0, -4.0]
    for al0 in al_guesses:
        dl = al0/np.sqrt(1+al0**2)
        s0 = (q[2]-q[0])/1.349/np.sqrt(1-2*dl**2/np.pi)
        m0 = q[1] - s0*B*dl
        try:
            r = least_squares(resid, [m0, np.log(max(s0,1e-8)), al0], jac=jac, args=(q, d),
                              bounds=([-np.inf, -np.inf, -alpha_max], [np.inf, np.inf, alpha_max]),
                              xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=200)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
        if best.cost < 1e-20:
            break
    m, ls, al = best.x
    return np.array([m, np.exp(ls), al]), best.cost

def nll(theta, y):
    m, ls, al = theta
    s = np.exp(ls); z = (y-m)/s
    return -(np.log(2) - ls + norm.logpdf(z) + norm.logcdf(al*z)).sum()

def nll_grad(theta, y):
    m, ls, al = theta
    s = np.exp(ls); z = (y-m)/s
    w = np.exp(norm.logpdf(al*z) - norm.logcdf(al*z))   # Mills-type ratio
    g_m = -(z/s - al*w/s).sum()
    g_ls = -(-1 + z**2 - al*z*w).sum()
    g_al = -(z*w).sum()
    return np.array([g_m, g_ls, g_al])

def ml_fit(y, alpha_max=50.0, init=None):
    n = len(y)
    best = None
    inits = []
    if init is not None: inits.append(init)
    mu, sd = y.mean(), y.std()
    g1 = ((y-mu)**3).mean()/sd**3
    g1c = np.clip(g1, -0.99, 0.99)
    # method-of-moments start
    r = np.sign(g1c)*(2*abs(g1c)/(4-np.pi))**(1/3)
    dl = np.clip(r/np.sqrt(2/np.pi*(1+r**2)), -0.99, 0.99)
    al0 = dl/np.sqrt(1-dl**2)
    s0 = sd/np.sqrt(1-2*dl**2/np.pi); m0 = mu - s0*B*dl
    inits += [[m0, np.log(s0), al0], [mu, np.log(sd), 0.5], [mu, np.log(sd), -0.5]]
    for th0 in inits:
        r = minimize(nll, th0, args=(y,), jac=nll_grad, method='L-BFGS-B',
                     bounds=[(None, None), (None, None), (-alpha_max, alpha_max)],
                     options=dict(maxiter=500, gtol=1e-9, ftol=1e-14))
        if best is None or r.fun < best.fun:
            best = r
    m, ls, al = best.x
    return np.array([m, np.exp(ls), al]), best.fun
