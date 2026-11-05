import numpy as np, itertools, sys, os, json, time
from multiprocessing import Pool
from scipy.stats import skewnorm
from scipy.special import owens_t, ndtr
from qmad import qmad_fit, ml_fit4, plugin_fit

OUT = "results"; os.makedirs(OUT, exist_ok=True)
XG = np.linspace(-8, 10, 3601); DX = XG[1]-XG[0]

def cdf(x, m, s, al):
    if not np.isfinite(s): return np.full_like(x, np.nan)
    al = np.clip(al, -1e4, 1e4); z = (x-m)/s
    return ndtr(z) - 2*owens_t(z, al)

def w1(th, al_true):
    return np.sum(np.abs(cdf(XG, *th[:3]) - cdf(XG, 0, 1, al_true)))*DX

def generate(n, al, eps, kind, k, rng):
    y = skewnorm.rvs(al, size=n, random_state=rng)
    if eps == 0 or kind == 'none': return y
    q1, q3 = skewnorm.ppf([.25, .75], al); iqr = q3-q1
    idx = np.flatnonzero(rng.random(n) < eps); nk = len(idx)
    jit = 0.25*iqr*rng.standard_normal(nk)
    if kind == 'R': y[idx] = q3 + k*iqr + np.abs(jit)
    elif kind == 'L': y[idx] = q1 - k*iqr - np.abs(jit)
    elif kind == 'B':
        side = rng.random(nk) < 0.5
        y[idx] = np.where(side, q3 + k*iqr + np.abs(jit), q1 - k*iqr - np.abs(jit))
    return y

def run_cell(cell):
    tag, n, al, eps, kind, k, R, methods, seed = cell
    fn = f"{OUT}/{tag}.npz"
    if os.path.exists(fn): return tag
    rng = np.random.default_rng(seed)
    est = {mth: np.empty((R, 4)) for mth in methods}
    wd = {mth: np.empty(R) for mth in methods}
    fns = dict(qmad=qmad_fit, ml=ml_fit4, plug=plugin_fit)
    for r in range(R):
        y = generate(n, al, eps, kind, k, rng)
        for mth in methods:
            th = fns[mth](y); est[mth][r] = th; wd[mth][r] = w1(th, al)
    np.savez(fn, **{f"est_{m}": est[m] for m in methods}, **{f"w1_{m}": wd[m] for m in methods},
             meta=json.dumps(dict(n=n, al=al, eps=eps, kind=kind, k=k, R=R)))
    return tag

def cells_main(R=500):
    out = []; s = 1000
    for n, al in itertools.product([100, 250, 500, 1000], [1, 3, 5]):
        s += 1
        out.append((f"main_n{n}_a{al}_e0_none", n, al, 0.0, 'none', 5, R, ['qmad', 'ml'], s))
        for eps, kind in itertools.product([.01, .025, .05, .10], ['R', 'L', 'B']):
            s += 1
            out.append((f"main_n{n}_a{al}_e{eps}_{kind}", n, al, eps, kind, 5, R, ['qmad', 'ml'], s))
    return out

def cells_dist(R=1000):
    out = []; s = 5000
    for kind, k in itertools.product(['R', 'L'], [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20]):
        s += 1
        out.append((f"dist_{kind}_k{k}", 500, 3, 0.05, kind, k, R, ['qmad', 'ml'], s))
    return out

def cells_plug(R=500):
    out = []; s = 9000
    for n, al in itertools.product([100, 500, 2000], [1, 2, 3, 5]):
        s += 1
        out.append((f"plug_n{n}_a{al}", n, al, 0.0, 'none', 5, R, ['plug', 'qmad', 'ml'], s))
    return out

if __name__ == "__main__":
    which = sys.argv[1]
    cells = dict(main=cells_main, dist=cells_dist, plug=cells_plug)[which]()
    cells.sort(key=lambda c: -c[1])  # big n first
    t = time.time()
    with Pool(2) as pool:
        for i, tag in enumerate(pool.imap_unordered(run_cell, cells)):
            print(f"{i+1}/{len(cells)} {tag} {time.time()-t:.0f}s", flush=True)
