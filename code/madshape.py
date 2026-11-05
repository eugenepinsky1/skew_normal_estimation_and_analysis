"""MAD skewness G=(mu-M)/H and MAD kurtosis K=(H_L+H_R)/(2H) for SN(0,1,alpha)."""
import numpy as np
from scipy.stats import norm
from snmad import F0, H
from qmad import xi
B = np.sqrt(2/np.pi)
def Lpart(x, al):
    c = np.sqrt(1+al**2); d = al/c
    return -2*norm.pdf(x)*norm.cdf(al*x) + B*d*norm.cdf(c*x)
def GK(al):
    q1, M, q3 = xi(al); mu = B*al/np.sqrt(1+al**2); Hm = H(M, al)
    HL = 2*(q1*(2*F0(q1, al)-F0(M, al)) - 2*Lpart(q1, al) + Lpart(M, al))
    HR = 2*(q3*(2*F0(q3, al)-F0(M, al)-1) - 2*Lpart(q3, al) + Lpart(M, al) + mu)
    return (mu-M)/Hm, (HL+HR)/(2*Hm)
def sample_GK(y):
    y = np.asarray(y, float); q1, M, q3 = np.quantile(y, [.25, .5, .75]); Hm = np.mean(np.abs(y-M))
    HL = np.mean(np.abs(y[y <= M]-q1)); HR = np.mean(np.abs(y[y > M]-q3))
    return (y.mean()-M)/Hm, (HL+HR)/(2*Hm)
if __name__ == "__main__":
    from scipy.stats import skewnorm
    for al in [0, 1, 2, 3, 5, 8, 1e4]:
        g, k = GK(al); y = skewnorm.rvs(al, size=2_000_000, random_state=1); gs, ks = sample_GK(y)
        print(al, round(g, 5), round(k, 5), "| MC", round(gs, 5), round(ks, 5))
