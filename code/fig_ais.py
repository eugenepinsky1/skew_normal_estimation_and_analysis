import json, numpy as np
from scipy.stats import skewnorm
from figstyle import *
from ais import df
o = json.load(open("results/ais.json")); y = df["BMI"].values
x = np.linspace(14, 36, 600); K = dict(qmad="QMAD", wqmad="WQMAD", ml="ML", qm="QM", mom="MoM", lm="LM")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
for a, (fits, title) in zip(ax, [(o["clean"], "original data ($n=202$)"),
                                 ({k: v["est"] for k, v in o["scen"]["2 errors /10"].items() if k != "values"},
                                  "two decimal-point errors (1%)")]):
    a.hist(y, bins=np.arange(16, 35.5, 1), density=True, color="#d9d9d9", edgecolor="white")
    for k in ORDER:
        th = fits[K[k]]; st = dict(STY[k]); st.pop("marker")
        a.plot(x, skewnorm.pdf(x, np.clip(th[2], -1e4, 1e4), loc=th[0], scale=th[1]), **st)
    a.set_title(title); a.set_xlabel("BMI (kg/m$^2$)")
ax[0].set_ylabel("density"); ax[0].legend(loc="upper right", fontsize=7)
fig.tight_layout(); fig.savefig("/home/claude/paper/figures/fig_ais.pdf"); fig.savefig("fig_ais.png", dpi=140)
