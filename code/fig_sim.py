import numpy as np, os
from figstyle import *
from analyze import cell
FD = "/home/claude/paper/figures"
def get(n, al, eps, kind):
    tag = f"main_n{n}_a{al}_e{eps if eps else 0}_{'none' if eps==0 else kind}"
    return cell(tag)[1]
def fig_eps(n=500):
    eps = [0, .01, .025, .05, .10]
    mets = [("rmse_s", r"RMSE of $\hat s$"), ("rmse_d", r"RMSE of $\hat\delta$"), ("w1", r"mean $W_1$")]
    kinds = [("R", "right tail"), ("L", "left tail"), ("B", "both tails")]
    for al in [1, 3, 5]:
        fig, ax = plt.subplots(3, 3, figsize=(7.2, 6.0), sharex=True)
        for j, (k, kl) in enumerate(kinds):
            res = [get(n, al, e, k) for e in eps]
            for i, (mk, ml) in enumerate(mets):
                a = ax[i, j]
                for meth in ORDER:
                    st = dict(STY[meth]); st["ms"] = 3.5
                    a.plot(eps, [r[meth][mk] for r in res], **st)
                if i == 0: a.set_title(kl)
                if j == 0: a.set_ylabel(ml)
                if i == 2: a.set_xlabel(r"contamination fraction $\varepsilon$")
                a.set_ylim(bottom=0)
        ax[0, 0].legend(loc="upper left", fontsize=7)
        fig.tight_layout(); fig.savefig(f"{FD}/fig_eps_a{al}.pdf"); fig.savefig(f"fig_eps_a{al}.png", dpi=130); plt.close(fig)
def fig_dist():
    ks = [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20]
    mets = [("rmse_s", r"RMSE of $\hat s$"), ("rmse_d", r"RMSE of $\hat\delta$"), ("w1", r"mean $W_1$")]
    fig, ax = plt.subplots(2, 3, figsize=(7.2, 4.6), sharex=True)
    for row, (kind, lab) in enumerate([("R", "right-tail outliers"), ("L", "left-tail outliers")]):
        res = [cell(f"dist_{kind}_k{k}")[1] for k in ks]
        for a, (mk, ml) in zip(ax[row], mets):
            for meth in ORDER:
                st = dict(STY[meth]); st["ms"] = 3.5
                a.plot(ks, [r[meth][mk] for r in res], **st)
            a.set_ylim(bottom=0)
            if row == 0: a.set_title(ml)
            if row == 1: a.set_xlabel(r"outlier distance $k$ (IQRs)")
        ax[row, 0].set_ylabel(lab)
    ax[0, 0].legend(loc="upper left", fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FD}/fig_dist.pdf"); fig.savefig("fig_dist.png", dpi=130); plt.close(fig)
if __name__ == "__main__":
    import sys
    if "dist" in sys.argv: fig_dist()
    if "eps" in sys.argv: fig_eps()
