import numpy as np, json, os, sys, glob, time
from sim import generate, w1, cells_main, cells_dist, cells_plug, OUT
from qmad import qm_fit, mom_fit
def run(cell):
    tag, n, al, eps, kind, k, R, methods, seed = cell
    fn = f"{OUT}/extra_{tag}.npz"
    if os.path.exists(fn) or not os.path.exists(f"{OUT}/{tag}.npz"): return
    rng = np.random.default_rng(seed)
    est = {m: np.empty((R, 4)) for m in ("qm", "mom")}; wd = {m: np.empty(R) for m in est}
    for r in range(R):
        y = generate(n, al, eps, kind, k, rng)
        for m, fn_ in (("qm", qm_fit), ("mom", mom_fit)):
            th = fn_(y); est[m][r] = th; wd[m][r] = w1(th, al)
    np.savez(fn, **{f"est_{m}": est[m] for m in est}, **{f"w1_{m}": wd[m] for m in est})
if __name__ == "__main__":
    cells = dict(main=cells_main, dist=cells_dist, plug=cells_plug)[sys.argv[1]]()
    t = time.time()
    for c in cells: run(c)
    print("done", sys.argv[1], time.time()-t)
