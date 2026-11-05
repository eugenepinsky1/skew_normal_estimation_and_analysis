import numpy as np, json
from scipy.stats import skewnorm
from scipy.integrate import quad
from snmad import F0
from qmad import xi, P
from wqmad import W, WQMAD
from asym import fisher
_cache = {}
def est(c):
    if c not in _cache: _cache[c] = WQMAD(c)
    return _cache[c]
def if_parts(al, c):
    q = xi(al); t = c*(q[2]-q[0]); f = skewnorm.pdf(q, al); Wv = W(q, t, al)
    a = 2*F0(q, al)-F0(q+t, al)-F0(q-t, al); b = 1-F0(q+t, al)+F0(q-t, al)
    def IFT(y):
        y = np.atleast_1d(np.asarray(y, float))
        IFq = [(P[j]-(y <= q[j]))/f[j] for j in range(3)]
        IFt = c*(IFq[2]-IFq[0])
        IFd = [np.minimum(np.abs(y-q[j]), t)-Wv[j]+a[j]*IFq[j]+b[j]*IFt for j in range(3)]
        return np.vstack(IFd+[IFq[1]])
    w = est(c); T0 = np.array([*Wv, q[1]]); h = 1e-6
    G = lambda T: w.from_stats(np.array([0, T[3], 0]), T[:3])[:3]
    J = np.column_stack([(G(T0+h*e)-G(T0-h*e))/(2*h) for e in np.eye(4)])
    pts = sorted(set(np.concatenate([q, q+t, q-t])))
    return IFT, J, pts
def if_wqmad(y, al, c):
    IFT, J, _ = if_parts(al, c); return J@IFT(y)
def avar_wqmad(al, c):
    IFT, J, pts = if_parts(al, c); edges = [-40]+pts+[40]
    V = np.zeros((4, 4))
    for i in range(4):
        for k in range(i, 4):
            V[i, k] = V[k, i] = sum(quad(lambda y: (IFT(y)[i]*IFT(y)[k])[0]*skewnorm.pdf(y, al), lo, hi, limit=200)[0]
                                     for lo, hi in zip(edges[:-1], edges[1:]))
    return J@V@J.T
if __name__ == "__main__":
    rows = json.load(open("results/avar.json"))
    for c in [1.0, 1.5, 2.0, 3.0]:
        out = []
        for r in rows:
            A = np.diag(avar_wqmad(r["al"], c)); r[f"wq{c}"] = A.tolist()
            out.append(f"a={r['al']}: " + " ".join(f"{x:.2f}" for x in A) + f" (QMAD {r['qmad'][1]:.2f} QM {r['qm'][1]:.2f})")
        print("c=", c, " | ".join(out))
    json.dump(rows, open("results/avar.json", "w"), indent=1)
