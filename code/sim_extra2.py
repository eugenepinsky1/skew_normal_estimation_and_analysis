import numpy as np, os, sys, time
from sim import generate, w1, cells_main, cells_dist, cells_plug, OUT
from wqmad import WQMAD
WQ = WQMAD(1.0)
def run(cell):
    tag, n, al, eps, kind, k, R, methods, seed = cell
    fn = f"{OUT}/extra2_{tag}.npz"
    if os.path.exists(fn) or not os.path.exists(f"{OUT}/{tag}.npz"): return
    rng = np.random.default_rng(seed); est = np.empty((R, 4)); wd = np.empty(R)
    for r in range(R):
        y = generate(n, al, eps, kind, k, rng); th = WQ.fit(y); est[r] = th; wd[r] = w1(th, al)
    np.savez(fn, est_wqmad=est, w1_wqmad=wd)
if __name__ == "__main__":
    cells = dict(main=cells_main, dist=cells_dist, plug=cells_plug)[sys.argv[1]]()
    t = time.time()
    for c in cells: run(c)
    print("done", sys.argv[1], time.time()-t)
