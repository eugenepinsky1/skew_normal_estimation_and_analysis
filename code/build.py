import json, numpy as np, re, os
from analyze import cell, f
R_MAIN = 500
s = open("main_src.tex").read()
def sci(x):
    if x < 1000: return f"{x:.2f}" if x < 100 else f"{x:.0f}"
    e = int(np.floor(np.log10(x))); return f"${x/10**e:.1f}\\times10^{{{e}}}$"
rows = json.load(open("results/avar.json"))
L = [r"\begin{table}[H]",
r"\caption{Asymptotic variances ($n\times$variance) at $\mathrm{SN}(0,1,\alpha)$ and asymptotic relative efficiency with respect to ML (ARE $=$ ML variance / variance). QMAD: quartile-MAD estimator; WQMAD: winsorized QMAD with $c=1$ (Section~\ref{sec:wqmad}); QM: quartile matching; MoM: method of moments; LM: L-moment estimator; ML: maximum likelihood estimator.\label{tab:avar}}",
r"\newcolumntype{C}{>{\centering\arraybackslash}X}",
r"\begin{tabularx}{\textwidth}{ccCCCCCC}", r"\toprule",
r" & & \multicolumn{3}{c}{\textbf{Asymptotic variance}} & \multicolumn{3}{c}{\textbf{ARE}}\\",
r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
r"$\alpha$ & Method & $m$ & $s$ & $\alpha$ & $m$ & $s$ & $\alpha$\\", r"\midrule"]
for r in rows:
    for k, lab in [("qmad", "QMAD"), ("wq1.0", "WQMAD"), ("qm", "QM"), ("mom", "MoM"), ("lm", "LM"), ("ml", "ML")]:
        v = r[k]; are = [r["ml"][i]/v[i] for i in range(3)]
        first = str(r["al"]) if k == "qmad" else ""
        L.append(f"{first} & {lab} & " + " & ".join(sci(x) for x in v) + " & " + " & ".join(f"{a:.3f}" for a in are) + r"\\")
    if r is not rows[-1]: L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
P2 = [r"\begin{table}[H]",
r"\caption{Asymptotic variances ($n\times$variance) of the plug-in estimator~\eqref{eq:plugin} and of QMAD at $\mathrm{SN}(0,1,\alpha)$. The plug-in Jacobian is singular at $\alpha^\ast=2.021$.\label{tab:plugavar}}",
r"\newcolumntype{C}{>{\centering\arraybackslash}X}",
r"\begin{tabularx}{\textwidth}{cCCCCCC}", r"\toprule",
r" & \multicolumn{3}{c}{\textbf{Plug-in}} & \multicolumn{3}{c}{\textbf{QMAD}}\\",
r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
r"$\alpha$ & $m$ & $s$ & $\alpha$ & $m$ & $s$ & $\alpha$\\", r"\midrule"]
for r in rows:
    P2.append(f"{r['al']} & " + " & ".join(sci(x) for x in r["plug"]+r["qmad"]) + r"\\")
P2 += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
s = s.replace("%%TABLE_PLUGAVAR%%", "\n".join(P2))
s = s.replace("%%TABLE_AVAR%%", "\n".join(L))
# ---- plug-in table ----
def plug_table():
    L = [r"\begin{table}[H]",
    r"\caption{Clean data ($\varepsilon=0$, %d replications): plug-in estimator~\eqref{eq:plugin} with quartile-consistent root selection, QMAD and ML. Wrong sign is the percentage of samples with $\sgn\hat\alpha\neq\sgn\alpha$. Plug: plug-in estimator; QMAD: quartile-MAD estimator; ML: maximum likelihood estimator.\label{tab:plug}}" % 500,
    r"\newcolumntype{C}{>{\centering\arraybackslash}X}",
    r"\small", r"\begin{tabularx}{\textwidth}{ccCCCCCCCCC}", r"\toprule",
    r" & & \multicolumn{3}{c}{\textbf{RMSE of $\hat\delta$}} & \multicolumn{3}{c}{\textbf{Wrong sign (\%)}} & \multicolumn{3}{c}{\textbf{Mean $W_1$}}\\",
    r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
    r"$n$ & $\alpha$ & Plug & QMAD & ML & Plug & QMAD & ML & Plug & QMAD & ML\\", r"\midrule"]
    for n in [100, 500, 2000]:
        for al in [1, 2, 3, 5]:
            _, r = cell(f"plug_n{n}_a{al}", ("plug", "qmad", "ml"))
            v = [f(r[m]["rmse_d"]) for m in ("plug", "qmad", "ml")] + [f(r[m]["wrong"], 1) for m in ("plug", "qmad", "ml")] + [f(r[m]["w1"]) for m in ("plug", "qmad", "ml")]
            L.append(f"{n} & {al} & " + " & ".join(v) + r"\\")
        if n != 2000: L.append(r"\midrule")
    L += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
    return "\n".join(L)



M4 = [("qmad", "QMAD"), ("wqmad", "WQMAD"), ("ml", "ML"), ("lm", "LM"), ("qm", "QM"), ("mom", "MoM")]
def best_fmt(vals, k):
    strs = [f"{v:.{k}f}" for v in vals]; b = min(vals)
    return [r"\textbf{"+s_+"}" if f"{v:.{k}f}" == f"{b:.{k}f}" else s_ for s_, v in zip(strs, vals)]
def main_tables(n=500, al=3):
    confs = [(0, 'none')]+[(e, k) for k in ['R', 'L', 'B'] for e in [.01, .025, .05, .10]]
    res = []
    for eps, kind in confs:
        res.append(((eps, kind), cell(f"main_n{n}_a{al}_e{eps if eps else 0}_{kind}")[1]))
    out = []
    for tabid, mets, cap in [
        ("tab:mainA", [("rmse_s", 3, r"RMSE of $\hat s$"), ("w1", 3, r"Mean $W_1$")],
         r"RMSE of the scale estimate and mean Wasserstein-1 distance between the fitted and the true distribution"),
        ("tab:mainB", [("rmse_d", 3, r"RMSE of $\hat\delta$"), ("wrong", 1, r"Wrong sign of $\hat\alpha$ (\%)")],
         r"RMSE of $\hat\delta=\hat\alpha/\sqrt{1+\hat\alpha^2}$ and percentage of samples with the wrong sign of $\hat\alpha$")]:
        L = [r"\begin{table}[H]",
             r"\caption{%s; $n=%d$, $\alpha=%d$, outliers at $k=5$ IQRs, %d replications. Contamination: R right tail, L left tail, B both tails. The smallest value in each row and block is in bold. Abbreviations: QMAD, quartile-MAD estimator; WQMAD, winsorized QMAD ($c=1$); MLE, maximum likelihood estimator; LM, L-moment estimator; QM, quartile matching; MoM, method of moments.\label{%s}}" % (cap, n, al, R_MAIN, tabid),
             r"\begin{adjustwidth}{-\extralength}{0cm}", r"\newcolumntype{C}{>{\centering\arraybackslash}X}",
             r"\footnotesize\setlength{\tabcolsep}{3pt}", r"\begin{tabularx}{\fulllength}{cc" + "C"*12 + "}", r"\toprule",
             r" & & " + " & ".join(r"\multicolumn{6}{c}{\textbf{%s}}" % m[2] for m in mets) + r"\\",
             r"\cmidrule(lr){3-8}\cmidrule(lr){9-14}",
             r"Type & $\varepsilon$ & " + " & ".join([lab for _, lab in M4]*2) + r"\\", r"\midrule"]
        prev = None
        for (eps, kind), r in res:
            if prev is not None and kind != prev: L.append(r"\midrule")
            prev = kind
            cells_ = []
            for key, k, _ in mets: cells_ += best_fmt([r[m][key] for m, _ in M4], k)
            L.append(("---" if eps == 0 else kind) + f" & {eps:g} & " + " & ".join(cells_) + r"\\")
        L += [r"\bottomrule", r"\end{tabularx}", r"\end{adjustwidth}", r"\end{table}"]
        out.append("\n".join(L))
    return out, dict(res)
def appendix_tables(eps=0.05):
    out = []
    for al in [1, 3, 5]:
        L = [r"\begin{table}[H]",
             r"\caption{Mean $W_1$ distance and RMSE of $\hat\delta$ for $\alpha=%d$, all sample sizes; clean data and $\varepsilon=%g$ with outliers at $k=5$ IQRs (R right, L left, B both tails); %d replications. The smallest value in each row and block is in bold. Abbreviations: QMAD, quartile-MAD estimator; WQMAD, winsorized QMAD ($c=1$); MLE, maximum likelihood estimator; LM, L-moment estimator; QM, quartile matching; MoM, method of moments.\label{tab:app%d}}" % (al, eps, R_MAIN, al),
             r"\begin{adjustwidth}{-\extralength}{0cm}", r"\newcolumntype{C}{>{\centering\arraybackslash}X}",
             r"\footnotesize\setlength{\tabcolsep}{3pt}", r"\begin{tabularx}{\fulllength}{cc" + "C"*12 + "}", r"\toprule",
             r" & & \multicolumn{6}{c}{\textbf{Mean $W_1$}} & \multicolumn{6}{c}{\textbf{RMSE of $\hat\delta$}}\\",
             r"\cmidrule(lr){3-8}\cmidrule(lr){9-14}",
             r"Data & $n$ & " + " & ".join([lab for _, lab in M4]*2) + r"\\", r"\midrule"]
        for e, kd, lab in [(0, "none", "Clean"), (eps, "R", "R"), (eps, "L", "L"), (eps, "B", "B")]:
            for n in [100, 250, 500, 1000]:
                r = cell(f"main_n{n}_a{al}_e{e if e else 0}_{kd}")[1]
                row = best_fmt([r[m]["w1"] for m, _ in M4], 3) + best_fmt([r[m]["rmse_d"] for m, _ in M4], 3)
                L.append((lab if n == 100 else "") + f" & {n} & " + " & ".join(row) + r"\\")
            if kd != "B": L.append(r"\midrule")
        L += [r"\bottomrule", r"\end{tabularx}", r"\end{adjustwidth}", r"\end{table}"]
        out.append("\n".join(L))
    return "\n\n".join(out)




if __name__ == "__main__":
    import sys
    fills = json.load(open("results/fill_static.json"))
    for extra in ("results/fill_dyn.json", "results/fill_text.json", "results/fill_gk.json", "results/fill_shapes.json"):
        if os.path.exists(extra): fills.update(json.load(open(extra)))
    if "full" in sys.argv:
        fills["%%TABLE_PLUG%%"] = plug_table()
        A, _ = main_tables(); fills["%%TABLE_MAINA%%"] = A[0]; fills["%%TABLE_MAINB%%"] = A[1]
        fills["%%TABLE_APP%%"] = appendix_tables()
    for k, v in fills.items(): s = s.replace(k, v)
    s = s.replace("maximum likelihood (ML)", "maximum likelihood estimation (MLE)")
    s = s.replace("ML estimates of", "MLEs of").replace("ML estimate of", "MLE of").replace("ML estimates", "MLEs").replace("ML estimate", "MLE")
    s = s.replace(r"\mathrm{ML}", r"\mathrm{MLE}")
    s = re.sub(r"(?<![A-Za-z\\])ML(?![A-Za-z])", "MLE", s)
    s = s.replace("MLE & Maximum likelihood\\\\", "MLE & Maximum likelihood estimation\\\\")
    s = s.replace(r"\mathrm{QMAD}", r"\mathrm{Q\text{-}MAD}").replace("WQMAD", "\\mbox{WQ-MAD}").replace("QMAD", "\\mbox{Q-MAD}")
    open("/home/claude/paper/main.tex", "w").write(s); open("/home/claude/mdpi_qmad/main.tex", "w").write(s)
    print("unfilled:", re.findall(r"%%[A-Z_0-9]+%%", s))
