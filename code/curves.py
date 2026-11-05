import numpy as np, json
from theory import if_qmad, if_ml, qmad_functional, ml_functional
al = 3.0
y = np.linspace(-6, 10, 401)
IFq = if_qmad(y, al); IFm = if_ml(y, al)
xs = np.linspace(-10, 15, 101)
bq = []; bm = []; init = [0, 0, al]
for x in xs:
    bq.append(qmad_functional(al, 0.05, x).tolist())
    th, init = ml_functional(al, 0.05, x, init); bm.append(th.tolist()); init = list(init)
json.dump(dict(al=al, y=y.tolist(), IFq=IFq.tolist(), IFm=IFm.tolist(), xs=xs.tolist(), bq=bq, bm=bm),
          open("results/curves.json", "w"))
for i in range(0, 101, 10): print(xs[i], np.round(bq[i], 3), np.round(bm[i], 3))
