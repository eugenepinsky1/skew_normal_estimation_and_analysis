import numpy as np, json, glob, os, csv
R = "results"
def load(tag):
    d = np.load(f"{R}/{tag}.npz"); meta = json.loads(str(d["meta"]))
    return d, meta
def metrics(est, w1, al):
    dt = al/np.sqrt(1+al**2)
    e = est[:, :2]-np.array([0, 1])
    out = dict(
        bias_m=np.nanmean(e[:, 0]), bias_s=np.nanmean(e[:, 1]),
        rmse_m=np.sqrt(np.nanmean(e[:, 0]**2)), rmse_s=np.sqrt(np.nanmean(e[:, 1]**2)),
        bias_d=np.nanmean(est[:, 3]-dt), rmse_d=np.sqrt(np.nanmean((est[:, 3]-dt)**2)),
        wrong=100*np.mean(np.sign(est[:, 2]) != np.sign(al)),
        mae_a=np.nanmedian(np.abs(est[:, 2]-al)), w1=np.nanmean(w1),
        inf=100*np.mean(np.abs(est[:, 2]) >= 99.9))
    return out
def cell(tag, methods=("qmad", "wqmad", "ml", "lm", "qm", "mom")):
    d, meta = load(tag); d = dict(d)
    for pre in ("extra_", "extra2_", "extra3_"):
        if os.path.exists(f"{R}/{pre}{tag}.npz"): d.update(dict(np.load(f"{R}/{pre}{tag}.npz")))
    return meta, {m: metrics(d[f"est_{m}"], d[f"w1_{m}"], meta["al"]) for m in methods if f"est_{m}" in d}

def f(x, k=3): return f"{x:.{k}f}"
def bold_pair(a, b, k=3, lower=True):
    sa, sb = f(a, k), f(b, k)
    if (a < b) == lower and sa != sb: sa = r"\textbf{"+sa+"}"
    elif sa != sb: sb = r"\textbf{"+sb+"}"
    return sa, sb

def main_table(n=500, al=3):
    rows = []
    for eps, kind in [(0, 'none')]+[(e, k) for k in ['R', 'L', 'B'] for e in [.01, .025, .05, .10]]:
        tag = f"main_n{n}_a{al}_e{eps if eps else 0}_{kind}"
        if not os.path.exists(f"{R}/{tag}.npz"): return None
        meta, r = cell(tag); q, m = r["qmad"], r["ml"]
        cells = []
        for key, k in [("rmse_m", 3), ("rmse_s", 3), ("rmse_d", 3), ("wrong", 1), ("w1", 3)]:
            cells += bold_pair(q[key], m[key], k)
        lab = "---" if eps == 0 else kind
        rows.append(f"{lab} & {eps:g} & " + " & ".join(cells) + r"\\")
    return rows

def all_rows():
    out = []
    for p in sorted(glob.glob(f"{R}/main_*.npz")):
        tag = os.path.basename(p)[:-4]; meta, r = cell(tag)
        for m in r: out.append(dict(**meta, method=m, **r[m]))
    return out

if __name__ == "__main__":
    rows = all_rows()
    if rows:
        with open(f"{R}/main_summary.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    t = main_table()
    print("\n".join(t) if t else "main not ready")
    for p in sorted(glob.glob(f"{R}/plug_*.npz")):
        tag = os.path.basename(p)[:-4]; meta, r = cell(tag, ("plug", "qmad", "ml"))
        print(tag, {m: (round(v['rmse_m'], 3), round(v['rmse_s'], 3), round(v['rmse_d'], 3), round(v['wrong'], 1), round(v['w1'], 3)) for m, v in r.items()})
