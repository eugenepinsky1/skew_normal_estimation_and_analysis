"""Winsorized quartile-MAD (WQMAD): MADs about quartiles with |y-Q_j| truncated at t = c*IQR."""
import numpy as np
from scipy.stats import skewnorm
from snmad import H, F0
from qmad import xi, P, _dgrid, _agrid

def W(a, t, al):
    """E min(|X-a|, t) for X ~ SN(0,1,al)."""
    return t + H(a, al) - 0.5*(H(a+t, al) + H(a-t, al))

def Wmod(al, c):
    x = xi(al); t = c*(x[2]-x[0]); return W(x, t, al), x
class WQMAD:
    def __init__(self, c=2.0):
        self.c = c
        tab = []
        for a in _agrid:
            Wv, x = Wmod(a, c); tab.append([*Wv, *x])
        self.tab = np.array(tab)
        self.kap = (self.tab[:, 2]-self.tab[:, 0])/self.tab[:, 1]
        bad = np.flatnonzero(np.diff(self.kap) <= 0)
        self.kb = bad[0] if len(bad) else len(self.kap)-1
        self.kmax = self.kap[self.kb]
    def kappa(self, al):
        Wv, _ = Wmod(al, self.c); return (Wv[2]-Wv[0])/Wv[1]
    def stats(self, y):
        q = np.quantile(y, P); t = self.c*(q[2]-q[0])
        d = np.minimum(np.abs(y[:, None]-q[None, :]), t).mean(0)
        return q, d
    def from_stats(self, q, d):
        k = (d[2]-d[0])/d[1]; sg = 1.0 if k >= 0 else -1.0
        if abs(k) >= self.kmax:
            row = self.tab[self.kb]; dl = _dgrid[self.kb]; al = dl/np.sqrt(1-dl**2)
            W2, xi2 = row[1], row[4]
        else:
            dl = np.interp(abs(k), self.kap[:self.kb+1], _dgrid[:self.kb+1]); al = dl/np.sqrt(1-dl**2)
            for _ in range(3):
                h = 1e-6*max(1, al); f0 = self.kappa(al)-abs(k); f1 = (self.kappa(al+h)-self.kappa(al-h))/(2*h)
                if f1 <= 0: break
                st = f0/f1; al = max(al-st, 0.0)
                if abs(st) < 1e-12*max(1, al): break
            Wv, x = Wmod(al, self.c); W2, xi2 = Wv[1], x[1]; dl = al/np.sqrt(1+al**2)
        s = d[1]/W2; m = q[1]-sg*s*xi2
        return np.array([m, s, sg*al, sg*dl])
    def fit(self, y):
        return self.from_stats(*self.stats(y))
