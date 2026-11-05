import numpy as np, json
from scipy.stats import skewnorm, norm
from scipy.integrate import quad
from scipy.optimize import brentq, minimize
from snmad import H, F0, dH_dal, jac as plug_jac
from qmad import xi, Dmod, kappa_of, qmad_from_stats, P
from asym import fisher, V_mad

def derivs(al, h=1e-5):
    D = lambda a: Dmod(a)
    kp = (kappa_of(al+h)-kappa_of(al-h))/(2*h)
    D2p = (D(al+h)[1]-D(al-h)[1])/(2*h)
    xi2p = (xi(al+h)[1]-xi(al-h)[1])/(2*h)
    return kp, D2p, xi2p

def if_qmad(y, al):
    """Influence function of (m, s, alpha) for QMAD at SN(0,1,al). y array."""
    y = np.asarray(y, float)
    q = xi(al); D = Dmod(al); f = skewnorm.pdf(q, al); kap = (D[2]-D[0])/D[1]
    kp, D2p, xi2p = derivs(al)
    IFd = [np.abs(y-q[j]) - D[j] + (2*P[j]-1)*(P[j]-(y <= q[j]))/f[j] for j in range(3)]
    IFq2 = (0.5-(y <= q[1]))/f[1]
    IFk = (IFd[2]-IFd[0]-kap*IFd[1])/D[1]
    IFa = IFk/kp
    IFs = (IFd[1]-D2p*IFa)/D[1]
    IFm = IFq2 - q[1]*IFs - xi2p*IFa
    return np.vstack([IFm, IFs, IFa])

def score(y, al):
    w = np.exp(norm.logpdf(al*y)-norm.logcdf(al*y))
    return np.vstack([y-al*w, -1+y**2-al*y*w, y*w])   # (m, log s = s at s=1, alpha)

def if_ml(y, al):
    return np.linalg.solve(fisher(al), score(np.asarray(y, float), al))

def avar_qmad(al):
    q = xi(al); pts = [-np.inf, *q, np.inf]
    A = np.empty((3, 3))
    for i in range(3):
        for j in range(i, 3):
            A[i, j] = A[j, i] = sum(quad(lambda y: (if_qmad(y, al)[i]*if_qmad(y, al)[j])[0]*skewnorm.pdf(y, al),
                                          a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
    return A

def avar_plug(al):
    q = xi(al); J = plug_jac([0, 0, al], q, np.ones(3)); Ji = np.linalg.inv(J)
    return Ji@V_mad(al, q)@Ji.T

# ---- asymptotic bias under point-mass contamination ----
def qmad_functional(al, eps, x):
    Fm = lambda t: (1-eps)*F0(t, al) + eps*(t >= x)
    q = np.empty(3)
    for j, p in enumerate(P):
        # generalized inverse of mixture cdf
        if (1-eps)*F0(x, al) >= p:  # quantile below x
            q[j] = brentq(lambda t: (1-eps)*F0(t, al)-p, -40, 40, xtol=1e-13)
        elif (1-eps)*F0(x, al)+eps >= p:
            q[j] = x
        else:
            q[j] = brentq(lambda t: (1-eps)*F0(t, al)+eps-p, -40, 40, xtol=1e-13)
    d = (1-eps)*H(q, al) + eps*np.abs(x-q)
    return qmad_from_stats(q, d)

_yg = np.linspace(-9, 9, 6001); _dy = _yg[1]-_yg[0]
def ml_functional(al, eps, x, init):
    w = skewnorm.pdf(_yg, al)*_dy; w /= w.sum()
    def obj(th):
        m, ls, a = th; s = np.exp(ls)
        lf = lambda y: np.log(2)-ls+norm.logpdf((y-m)/s)+norm.logcdf(a*(y-m)/s)
        return -((1-eps)*np.dot(w, lf(_yg)) + eps*lf(x))
    best = None
    for th0 in [init, [0, 0, al], [0, 0, -al], [0, 0.5, 1]]:
        r = minimize(obj, th0, method='L-BFGS-B', bounds=[(None, None), (-3, 3), (-100, 100)],
                     options=dict(ftol=1e-15, gtol=1e-10, maxiter=2000))
        if best is None or r.fun < best.fun: best = r
    m, ls, a = best.x
    return np.array([m, np.exp(ls), a, a/np.sqrt(1+a**2)]), best.x

if __name__ == "__main__":
    out = {}
    print("alpha* =", brentq(lambda a: np.linalg.det(plug_jac([0, 0, a], xi(a), np.ones(3))), 1.9, 2.2, xtol=1e-12))
    rows = []
    for al in [1, 2, 3, 5, 8]:
        Aq = avar_qmad(al); Am = np.linalg.inv(fisher(al)); Ap = avar_plug(al)
        rows.append(dict(al=al, qmad=np.diag(Aq).tolist(), ml=np.diag(Am).tolist(), plug=np.diag(Ap).tolist()))
        print(al, "QMAD", np.round(np.diag(Aq), 3), "ML", np.round(np.diag(Am), 3), "PLUG", np.round(np.diag(Ap), 1),
              "ARE", np.round(np.diag(Am)/np.diag(Aq), 3))
    json.dump(rows, open("results/avar.json", "w"), indent=1)
