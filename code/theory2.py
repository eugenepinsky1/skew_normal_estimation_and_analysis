import numpy as np, json
from scipy.stats import skewnorm
from scipy.integrate import quad
import qmad
from qmad import xi, P, B
from theory import if_qmad, if_ml, avar_qmad
from asym import fisher

def G_qm(q):
    b = (q[2]+q[0]-2*q[1])/(q[2]-q[0])
    dl = np.interp(abs(b), qmad._bow, qmad._dgb); al = dl/np.sqrt(1-dl**2)
    from scipy.optimize import brentq
    bb = lambda a: (lambda x: (x[2]+x[0]-2*x[1])/(x[2]-x[0]))(xi(a))
    al = brentq(lambda a: bb(a)-b, max(al*0.8, 1e-4), al*1.2+0.1, xtol=1e-14)
    x = xi(al); s = (q[2]-q[0])/(x[2]-x[0]); return np.array([q[1]-s*x[1], s, al])
def G_mom(r):   # raw moments m1,m2,m3
    m1, m2, m3 = r; v = m2-m1**2; mu3 = m3-3*m1*m2+2*m1**3; g = mu3/v**1.5
    rr = (2*g/(4-np.pi))**(1/3); dl = np.sqrt(np.pi/2)*rr/np.sqrt(1+rr**2)
    s = np.sqrt(v)/np.sqrt(1-2*dl**2/np.pi); return np.array([m1-s*B*dl, s, dl/np.sqrt(1-dl**2)])
def jacnum(G, x, h=1e-6):
    return np.column_stack([(G(x+h*e)-G(x-h*e))/(2*h) for e in np.eye(len(x))])
def avar_qm(al):
    q = xi(al); f = skewnorm.pdf(q, al)
    V = np.array([[(min(P[i], P[j])-P[i]*P[j])/(f[i]*f[j]) for j in range(3)] for i in range(3)])
    J = jacnum(G_qm, q); return J@V@J.T
def raw(k, al): return skewnorm.moment(k, al)
def avar_mom(al):
    r = np.array([raw(k, al) for k in (1, 2, 3)])
    V = np.array([[raw(i+j, al)-r[i-1]*r[j-1] for j in (1, 2, 3)] for i in (1, 2, 3)])
    J = jacnum(G_mom, r); return J@V@J.T
def if_qm(y, al):
    q = xi(al); f = skewnorm.pdf(q, al)
    IFq = np.vstack([(P[j]-(y <= q[j]))/f[j] for j in range(3)])
    return jacnum(G_qm, q)@IFq
def if_mom(y, al):
    r = np.array([raw(k, al) for k in (1, 2, 3)])
    IFr = np.vstack([y**k-r[k-1] for k in (1, 2, 3)])
    return jacnum(G_mom, r)@IFr
if __name__ == "__main__":
    rows = json.load(open("results/avar.json"))
    for r in rows:
        al = r["al"]; r["qm"] = np.diag(avar_qm(al)).tolist(); r["mom"] = np.diag(avar_mom(al)).tolist()
        print(al, "QMAD", np.round(r["qmad"], 2), "QM", np.round(r["qm"], 2), "MOM", np.round(r["mom"], 2), "ML", np.round(r["ml"], 2))
    json.dump(rows, open("results/avar.json", "w"), indent=1)
    c = json.load(open("results/curves.json")); y = np.array(c["y"])
    c["IFqm"] = if_qm(y, c["al"]).tolist(); c["IFmom"] = if_mom(y, c["al"]).tolist()
    json.dump(c, open("results/curves.json", "w"))
