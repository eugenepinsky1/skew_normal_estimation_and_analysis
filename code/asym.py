import numpy as np
from scipy.stats import skewnorm, norm
from scipy.integrate import quad
from snmad import *
p=np.array([.25,.5,.75])
def E(f, al):
    return quad(lambda x: f(x)*skewnorm.pdf(x,al), -np.inf, np.inf, limit=200, epsabs=1e-11)[0]
def V_mad(al, q):
    D = H(q, al)
    V=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            V[i,j]=E(lambda x: abs(x-q[i])*abs(x-q[j]), al)-D[i]*D[j]
    return V
def fisher(al):
    def sc(x):
        w = np.exp(norm.logpdf(al*x)-norm.logcdf(al*x))
        return np.array([x-al*w, -1+x**2-al*x*w, x*w])  # m, log s, alpha at m=0,s=1
    I=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            I[i,j]=E(lambda x: sc(x)[i]*sc(x)[j], al)
    return I
def asym_qmad(al):
    q=skewnorm.ppf(p,al)
    J=jac([0,0,al],q,np.ones(3))
    Ji=np.linalg.inv(J)
    return Ji@V_mad(al,q)@Ji.T
if __name__=="__main__":
    for al in [0.5,1,1.5,2,2.03,2.5,3,4,5,8]:
        A=asym_qmad(al); Iinv=np.linalg.inv(fisher(al))
        # log s -> s same at s=1
        print(al, "qmad var(m,s,a)=", np.round(np.diag(A),2), " ML=", np.round(np.diag(Iinv),2))
