import json, numpy as np
from figstyle import *
c = json.load(open("results/curves.json"))
y = np.array(c["y"]); IFq = np.array(c["IFq"]); IFm = np.array(c["IFm"])
names = [r"location $m$", r"scale $s$", r"shape $\alpha$"]
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
IFs = {"qmad": IFq, "wqmad": np.array(c["IFwq"]), "ml": IFm, "qm": np.array(c["IFqm"]), "mom": np.array(c["IFmom"]), "lm": np.array(c["IFlm"])}
for i in range(3):
    for k in ORDER:
        st = dict(STY[k]); st.pop("marker")
        ax[i].plot(y, IFs[k][i], **st)
    ax[i].set_title(names[i]); ax[i].set_xlabel(r"contaminating point $y$")
    ax[i].axhline(0, color="#999", lw=0.6)
lim = [(-12, 12), (-10, 60), (-400, 400)]
for i in range(3): ax[i].set_ylim(*lim[i])
ax[0].set_ylabel("influence function"); ax[0].legend(loc="upper left", fontsize=7)
fig.savefig("/home/claude/paper/figures/fig_if.pdf"); fig.savefig("fig_if.png", dpi=150)

xs = np.array(c["xs"]); bq = np.array(c["bq"]); bm = np.array(c["bm"]); al = c["al"]
dt = al/np.sqrt(1+al**2); tru = [0, 1, None, dt]
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
B = {"qmad": bq, "wqmad": np.array(c["bwq"]), "ml": bm, "qm": np.array(c["bqm"]), "mom": np.array(c["bmom"]), "lm": np.array(c["blm"])}
for a, j, t in zip(ax, [0, 1, 3], [r"$m$", r"$s$", r"$\delta=\alpha/\sqrt{1+\alpha^2}$"]):
    for k in ORDER:
        st = dict(STY[k]); st.pop("marker")
        a.plot(xs, B[k][:, j]-tru[j], **st)
    a.axhline(0, color="#999", lw=0.6); a.set_title("bias of "+t); a.set_xlabel(r"outlier position $x$")
ax[0].set_ylim(-3, 4); ax[1].set_ylim(-0.5, 4.5)
ax[0].legend(loc="upper right", fontsize=7)
fig.savefig("/home/claude/paper/figures/fig_bias.pdf"); fig.savefig("fig_bias.png", dpi=150)
