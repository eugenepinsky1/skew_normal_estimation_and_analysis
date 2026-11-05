import numpy as np, pandas as pd, json
from scipy.stats import skewnorm
from ais import df, F
from sim import XG
y0 = df["BMI"].values.astype(float); n = len(y0)
KEYS = ["QMAD", "WQMAD", "QM", "MoM", "LM", "ML"]
def fit(y): return {k: F[k](y) for k in KEYS}
def w1(a, b):
    lo, hi = 10, 40; x = np.linspace(lo, hi, 6001)
    cdf = lambda th: skewnorm.cdf(x, np.clip(th[2], -1e4, 1e4), loc=th[0], scale=th[1])
    return np.trapezoid(np.abs(cdf(a)-cdf(b)), x)
out = {"n": n, "summary": dict(mean=y0.mean(), sd=y0.std(ddof=1), g1=float(((y0-y0.mean())**3).mean()/y0.std()**3),
        q=np.quantile(y0, [.25, .5, .75]).tolist(), min=y0.min(), max=y0.max())}
clean = fit(y0); out["clean"] = {k: v.tolist() for k, v in clean.items()}
rng = np.random.default_rng(2026); B = 500; boot = {k: [] for k in KEYS}
for b in range(B):
    yb = y0[rng.integers(0, n, n)]
    for k, v in fit(yb).items(): boot[k].append(v)
out["boot_se"] = {}
for k in KEYS:
    a = np.array(boot[k]); out["boot_se"][k] = [float(np.std(a[:, 0])), float(np.std(a[:, 1])), float(np.std(a[:, 3])),
                                                 float(np.median(np.abs(a[:, 2]-np.median(a[:, 2]))))]
    out["boot_se"][k].append(float(np.mean(np.sign(a[:, 2]) != np.sign(clean[k][2]))*100))
# corrupted data: decimal-point errors
rng = np.random.default_rng(7); idx = rng.choice(n, 4, replace=False)
scen = {}
for lab, k_, op in [("2 errors /10", 2, "div"), ("4 errors /10", 4, "div"), ("2 errors x10", 2, "mul")]:
    y = y0.copy(); j = idx[:k_]; y[j] = y[j]/10 if op == "div" else y[j]*10
    f = fit(y); scen[lab] = {k: dict(est=f[k].tolist(), w1=w1(f[k], clean[k]), w1_ml=w1(f[k], clean["ML"])) for k in KEYS}
    scen[lab]["values"] = y[j].tolist()
out["scen"] = scen; out["idx"] = idx.tolist(); out["orig_vals"] = y0[idx].tolist()
from madshape import GK, sample_GK
g, k = sample_GK(y0); rng = np.random.default_rng(11)
bs = np.array([sample_GK(y0[rng.integers(0, n, n)]) for _ in range(2000)])
out["GK_sample"] = [g, k]; out["GK_ci"] = [np.percentile(bs[:, 0], [2.5, 97.5]).tolist(), np.percentile(bs[:, 1], [2.5, 97.5]).tolist()]
out["GK_fit"] = {kk: list(GK(min(abs(v[2]), 1e4)*np.sign(v[2]))) for kk, v in out["clean"].items()}
json.dump(out, open("results/ais.json", "w"), indent=1, default=float)
print(json.dumps(out["summary"], default=float)); 
for k in KEYS: print(k, np.round(clean[k], 3), np.round(out["boot_se"][k], 3))
for lab, d in scen.items():
    print(lab, d["values"]); 
    for k in KEYS: print("   ", k, np.round(d[k]["est"], 3), round(d[k]["w1"], 3))
