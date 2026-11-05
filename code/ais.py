import numpy as np, pandas as pd, json
from qmad import qmad_fit, qm_fit, mom_fit, ml_fit4
from wqmad import WQMAD
from lmom import lm_fit
WQ = WQMAD(1.0)
df = pd.read_csv("/home/claude/data/ais_sn.csv")
F = dict(QMAD=qmad_fit, WQMAD=WQ.fit, QM=qm_fit, MoM=mom_fit, ML=ml_fit4, LM=lm_fit)
def fit_all(y): return {k: f(np.asarray(y, float)) for k, f in F.items()}
if __name__ == "__main__":
    for var in ["BMI", "Bfat", "SSF", "LBM", "Fe", "Wt"]:
        for grp, d in [("all", df), ("female", df[df.sex == "female"]), ("male", df[df.sex == "male"])]:
            y = d[var].values; r = fit_all(y)
            g1 = ((y-y.mean())**3).mean()/y.std()**3
            print(f"{var:5s} {grp:6s} n={len(y)} g1={g1:.2f} ", " | ".join(f"{k}: {v[0]:.2f},{v[1]:.2f},{v[2]:.2f}" for k, v in r.items()))
