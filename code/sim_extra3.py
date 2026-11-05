import numpy as np, os, sys, time
from sim import generate, w1, cells_main, cells_dist, cells_plug, OUT
from lmom import lm_fit

def run(cell):
    tag, n, al, eps, kind, k, R, methods, seed = cell
    fn = f"{OUT}/extra3_{tag}.npz"
    if os.path.exists(fn) or not os.path.exists(f"{OUT}/{tag}.npz"): return
    rng = np.random.default_rng(seed); est = np.empty((R, 4)); wd = np.empty(R)
    for r in range(R):
        y = generate(n, al, eps, kind, k, rng); th = lm_fit(y); est[r] = th; wd[r] = w1(th, al)
    np.savez(fn, est_lm=est, w1_lm=wd)
if __name__ == "__main__":
    cells = dict(main=cells_main, dist=cells_dist, plug=cells_plug)[sys.argv[1]]()
    t = time.time()
    for c in cells: run(c)
    print("done", sys.argv[1], time.time()-t)
