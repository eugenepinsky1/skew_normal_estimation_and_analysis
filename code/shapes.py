import numpy as np, json
from scipy.stats import skewnorm, skew, kurtosis
from sim import generate
from lmom import sample_lmom
from madshape import sample_GK
def lm4(y):
    x = np.sort(np.asarray(y, float)); n = len(x); i = np.arange(n, dtype=float)
    b0 = x.mean()
    b1 = np.sum(i*x)/(n*(n-1))
    b2 = np.sum(i*(i-1)*x)/(n*(n-1)*(n-2))
    b3 = np.sum(i*(i-1)*(i-2)*x)/(n*(n-1)*(n-2)*(n-3))
    l2 = 2*b1-b0; l3 = 6*b2-6*b1+b0; l4 = 20*b3-30*b2+12*b1-b0
    return l3/l2, l4/l2
def measures(y):
    q1, M, q3 = np.quantile(y, [.25, .5, .75]); d = np.abs(y[:, None]-np.array([q1, M, q3])).mean(0)
    G, K = sample_GK(y); t3, t4 = lm4(y)
    return dict(g1=skew(y), b2=kurtosis(y, fisher=False), t3=t3, t4=t4, G=G, K=K,
                kappa=(d[2]-d[0])/d[1], bowley=(q3+q1-2*M)/(q3-q1))
if __name__ == "__main__":
    al, n, R = 3.0, 500, 1000
    pop = measures(skewnorm.rvs(al, size=2_000_000, random_state=0))
    out = {"pop": pop}
    for tag, eps, kind in [("clean", 0, "none"), ("L1", .01, "L"), ("L5", .05, "L"), ("R1", .01, "R"), ("R5", .05, "R")]:
        rng = np.random.default_rng(777); M = []
        for r in range(R): M.append(measures(generate(n, al, eps, kind, 5, rng)))
        out[tag] = {k: [float(np.mean([m[k] for m in M])), float(np.std([m[k] for m in M]))] for k in pop}
    json.dump(out, open("results/shapes.json", "w"), indent=1, default=float)
    for k in pop: print(k, round(pop[k], 4), {t: [round(v, 3) for v in out[t][k]] for t in out if t != "pop"})
