"""Corrected quartile-MAD (QMAD) estimator for SN(m, s, alpha) + naive plug-in variant."""
import numpy as np
from scipy.stats import skewnorm
from scipy.optimize import brentq, least_squares
from snmad import H, F0, dH_da, dH_dal, jac as plug_jac, resid as plug_resid, sample_stats, ml_fit, B

P = np.array([.25, .5, .75])

def xi(al):
    """Model quartiles of standard SN(0,1,al) (numerical root of the Owen-T CDF)."""
    return np.array([brentq(lambda z: F0(z, al) - p, -12, 12, xtol=1e-14) for p in P])

def Dmod(al):
    return H(xi(al), al)

def kappa_of(al):
    D = Dmod(al); return (D[2]-D[0])/D[1]

# ---------- lookup tables on delta in [0,1] ----------
_dgrid = np.concatenate([np.linspace(0, 0.9, 1801)[:-1], 1-np.logspace(np.log10(0.1), -9, 2200)])
_agrid = _dgrid/np.sqrt(1-_dgrid**2)
_tab = np.array([[*Dmod(a), *xi(a)] for a in _agrid])
_kap = (_tab[:, 2]-_tab[:, 0])/_tab[:, 1]
KAPPA_MAX = _kap[-1]
assert np.all(np.diff(_kap) > 0)

def _interp_al(k):
    """Invert kappa -> alpha (alpha>=0), with Newton polish."""
    d = np.interp(k, _kap, _dgrid)
    if d >= _dgrid[-1]:
        return np.inf
    al = d/np.sqrt(1-d**2)
    for _ in range(3):  # secant/Newton polish on kappa(al)=k
        h = 1e-6*max(1, al)
        f0 = kappa_of(al)-k; f1 = (kappa_of(al+h)-kappa_of(al-h))/(2*h)
        if f1 <= 0: break
        step = f0/f1
        al = max(al-step, 0.0)
        if abs(step) < 1e-12*max(1, al): break
    return al

def qmad_from_stats(q, d):
    k = (d[2]-d[0])/d[1]
    sgn = 1.0 if k >= 0 else -1.0
    if abs(k) >= KAPPA_MAX:
        al = np.inf; D2 = _tab[-1, 1]; xi2 = _tab[-1, 4]; dl = 1.0
    else:
        al = _interp_al(abs(k)); D = Dmod(al); D2 = D[1]; xi2 = xi(al)[1]; dl = al/np.sqrt(1+al**2)
    s = d[1]/D2
    m = q[1] - s*sgn*xi2
    return np.array([m, s, sgn*al, sgn*dl])

def qmad_fit(y):
    q, d = sample_stats(y)
    return qmad_from_stats(q, d)

# ---------- naive plug-in variant (empirical quartiles inside theoretical MAD) ----------
def plugin_fit(y, alpha_max=100.0):
    q, d = sample_stats(y); best = None
    for al0 in [-8, -4, -2, -1, -.5, .5, 1, 2, 4, 8]:
        dl = al0/np.sqrt(1+al0**2); s0 = (q[2]-q[0])/1.349/np.sqrt(1-2*dl**2/np.pi); m0 = q[1]-s0*B*dl
        r = least_squares(plug_resid, [m0, np.log(s0), al0], jac=plug_jac, args=(q, d),
                          bounds=([-np.inf, -np.inf, -alpha_max], [np.inf, np.inf, alpha_max]),
                          xtol=1e-12, ftol=1e-12, gtol=1e-12)
        th = np.array([r.x[0], np.exp(r.x[1]), r.x[2]])
        mis = np.sum((F0((q-th[0])/th[1], th[2])-P)**2)
        crit = (r.cost > 1e-12, mis)    # exact roots first, then quartile consistency
        if best is None or crit < best[0]: best = (crit, th)
    th = best[1]
    return np.array([*th, th[2]/np.sqrt(1+th[2]**2)])

def ml_fit4(y, alpha_max=100.0):
    st = qmad_fit(y)
    init = None
    if np.isfinite(st[2]):
        init = [st[0], np.log(st[1]), float(np.clip(st[2], -alpha_max, alpha_max))]
    th, _ = ml_fit(y, alpha_max=alpha_max, init=init)
    return np.array([*th, th[2]/np.sqrt(1+th[2]**2)])

# ---------- comparators ----------
_bow = (_tab[:, 5]+_tab[:, 3]-2*_tab[:, 4])/(_tab[:, 5]-_tab[:, 3])   # Bowley skewness of model
_kb = np.flatnonzero(np.diff(_bow) <= 0)[0]   # Bowley flat to machine precision beyond here
_bow = _bow[:_kb+1]; _dgb = _dgrid[:_kb+1]
BOWLEY_MAX = _bow[-1]

def _al_from_delta(d):
    return np.inf if d >= 1 else d/np.sqrt(1-d**2)

def qm_fit(y):
    """Quartile matching: model quartiles = sample quartiles."""
    q = np.quantile(y, [0.25, 0.5, 0.75])
    b = (q[2]+q[0]-2*q[1])/(q[2]-q[0]); sg = 1.0 if b >= 0 else -1.0
    if abs(b) >= BOWLEY_MAX:
        dl = 1.0; xi_ = _tab[-1, 3:6]; al = np.inf
    else:
        dl = np.interp(abs(b), _bow, _dgb); al = _al_from_delta(dl)
        for _ in range(3):
            h = 1e-6*max(1, al); bb = lambda a: (lambda x: (x[2]+x[0]-2*x[1])/(x[2]-x[0]))(xi(a))
            f0 = bb(al)-abs(b); f1 = (bb(al+h)-bb(al-h))/(2*h)
            if f1 <= 0: break
            st = f0/f1; al = max(al-st, 0.0)
            if abs(st) < 1e-12*max(1, al): break
        xi_ = xi(al); dl = al/np.sqrt(1+al**2)
    s = (q[2]-q[0])/(xi_[2]-xi_[0])
    xi2 = xi_[1] if sg > 0 else -xi_[1]
    m = q[1]-s*xi2
    return np.array([m, s, sg*al, sg*dl])

G1_MAX = (4-np.pi)/2*(2/np.pi)**1.5/(1-2/np.pi)**1.5   # skewness of half-normal

def mom_fit(y):
    mu = y.mean(); sd = y.std(ddof=0); g = np.mean((y-mu)**3)/sd**3
    sg = 1.0 if g >= 0 else -1.0
    if abs(g) >= G1_MAX:
        dl = 1.0
    else:
        r = (2*abs(g)/(4-np.pi))**(1/3)
        dl = np.sqrt(np.pi/2)*r/np.sqrt(1+r**2)
    s = sd/np.sqrt(1-2*dl**2/np.pi); m = mu-sg*s*B*dl
    return np.array([m, s, sg*_al_from_delta(dl), sg*dl])
