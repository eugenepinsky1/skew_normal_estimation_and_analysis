"""L-moment estimator for SN(m, s, alpha): match l1, l2, t3 = l3/l2."""
import numpy as np
from scipy.stats import skewnorm
from scipy.integrate import quad
from qmad import _dgrid, B

from snmad import F0
_XG = np.linspace(-12, 12, 24001)
def lmom_std(al):
    """lambda_1..3 of SN(0,1,al): lambda2 = int F(1-F), lambda3 = int F(1-F)(2F-1)."""
    F = F0(_XG, al); g = F*(1-F)
    l1 = B*al/np.sqrt(1+al**2)
    return np.array([l1, np.trapezoid(g, _XG), np.trapezoid(g*(2*F-1), _XG)])

_ag = _dgrid/np.sqrt(1-_dgrid**2)
_LT = np.array([lmom_std(a) for a in np.minimum(_ag, 1e4)])
_t3 = _LT[:, 2]/_LT[:, 1]
_kb = np.flatnonzero(np.diff(_t3) <= 0)
_kb = _kb[0] if len(_kb) else len(_t3)-1
T3_MAX = _t3[_kb]

def sample_lmom(y):
    x = np.sort(y); n = len(x); i = np.arange(n)
    b0 = x.mean(); b1 = np.sum(i*x)/(n*(n-1)); b2 = np.sum(i*(i-1)*x)/(n*(n-1)*(n-2))
    return np.array([b0, 2*b1-b0, 6*b2-6*b1+b0])

def from_lmom(L):
    t3 = L[2]/L[1]; sg = 1.0 if t3 >= 0 else -1.0
    if abs(t3) >= T3_MAX:
        l = _LT[_kb]; al = np.inf; dl = 1.0
    else:
        dl = np.interp(abs(t3), _t3[:_kb+1], _dgrid[:_kb+1])
        l = np.array([np.interp(dl, _dgrid, _LT[:, k]) for k in range(3)])
        al = dl/np.sqrt(1-dl**2)
    s = L[1]/l[1]; m = L[0]-s*sg*l[0]
    return np.array([m, s, sg*al, sg*dl])

def lm_fit(y):
    return from_lmom(sample_lmom(np.asarray(y, float)))
