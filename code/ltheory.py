import numpy as np, json
from scipy.stats import skewnorm
from snmad import F0
import lmom
X = np.linspace(-12, 12, 6001); DX = X[1]-X[0]
def if_lambda(y, al):
    F = F0(X, al); mu = lmom.lmom_std(al)[0]
    d2 = 1-2*F; d3 = -(2*F-1)**2+2*F*(1-F)
    y = np.atleast_1d(y)
    ind = (y[:, None] <= X[None, :]).astype(float)
    I2 = ((ind-F[None, :])*d2[None, :]).sum(1)*DX
    I3 = ((ind-F[None, :])*d3[None, :]).sum(1)*DX
    return np.vstack([y-mu, I2, I3])
def G(L):
    return lmom.from_lmom(np.array(L))[:3]
def jac(al, h=1e-6):
    L0 = lmom.lmom_std(al)
    # use exact (non-interpolated) inverse for derivative: finite differences on the table-based map are fine at this h? use larger h
    h = 1e-5
    return np.column_stack([(G(L0+h*e)-G(L0-h*e))/(2*h) for e in np.eye(3)])
def avar_lm(al):
    yg = np.linspace(-8, 8, 2001); w = skewnorm.pdf(yg, al)*(yg[1]-yg[0])
    IF = if_lambda(yg, al); S = (IF*w)@IF.T
    J = jac(al); return J@S@J.T
def if_lm(y, al):
    return jac(al)@if_lambda(y, al)
if __name__ == "__main__":
    rows = json.load(open("results/avar.json"))
    for r in rows:
        A = np.diag(avar_lm(r["al"])); r["lm"] = A.tolist()
        print(r["al"], np.round(A, 3), "ML", np.round(r["ml"], 3), "WQ", np.round(r["wq1.0"], 3))
    json.dump(rows, open("results/avar.json", "w"), indent=1)
    rng = np.random.default_rng(1); n = 2000
    E = np.array([lmom.lm_fit(skewnorm.rvs(3, size=n, random_state=rng)) for _ in range(2000)])
    print("MC n*var alpha=3:", n*E[:, :3].var(0))
