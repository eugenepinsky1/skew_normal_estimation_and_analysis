import json, numpy as np
o = json.load(open("results/ais.json")); F = {}
S = o["summary"]
F["%%AIS_MEAN%%"] = f"{S['mean']:.2f}"; F["%%AIS_SD%%"] = f"{S['sd']:.2f}"; F["%%AIS_G1%%"] = f"{S['g1']:.2f}"
F["%%AIS_Q%%"] = "$%.2f$, $%.2f$ and $%.2f$" % tuple(S["q"]); F["%%AIS_RANGE%%"] = f"${S['min']:.2f}$--${S['max']:.2f}$"
C = o["clean"]; BS = o["boot_se"]; SC = o["scen"]
F["%%AIS_MOM_A%%"] = f"{C['MoM'][2]:.2f}"
for k, key in [("WQMAD", "WQ"), ("QMAD", "QMAD"), ("QM", "QM")]:
    F[f"%%AIS_{key}_FLIP_C%%"] = f"{100-BS[k][4]:.0f}"
F["%%AIS_IDX%%"] = ", ".join(f"${v:.2f}$" for v in o["orig_vals"])
s1, s2, s3 = SC["2 errors /10"], SC["4 errors /10"], SC["2 errors x10"]
F["%%AIS_S1_ML%%"] = f"{s1['ML']['est'][3]:.2f}"; F["%%AIS_S1_MOM%%"] = f"{s1['MoM']['est'][3]:.2f}"
F["%%AIS_S1_ML_W%%"] = f"{s1['ML']['w1']:.2f}"; F["%%AIS_S2_WQ%%"] = f"{s2['WQMAD']['est'][3]:.2f}"
F["%%AIS_ML_S%%"] = f"{C['ML'][1]:.2f}"; F["%%AIS_S3_ML_S%%"] = f"{s3['ML']['est'][1]:.2f}"
F["%%AIS_S3_QMAD_S%%"] = f"{s3['QMAD']['est'][1]:.2f}"
F["%%AIS_S3_MAXW%%"] = f"{max(s3['WQMAD']['w1'], s3['QM']['w1']):.2f}"
KEYS = ["QMAD", "WQMAD", "QM", "MoM", "LM", "ML"]
def fa(x): return r"$\infty$" if not np.isfinite(x) else (f"${x:.2f}$" if abs(x) < 99.9 else r"$100^{\dagger}$")
T = [r"\begin{table}[H]",
 r"\caption{Skew-normal fits to the BMI of $202$ AIS athletes. Bootstrap standard errors (500 resamples) are in parentheses. The last column is the percentage of bootstrap resamples in which $\sgn\hat\alpha$ agrees with the estimate from the full sample. Abbreviations: QMAD, quartile-MAD estimator; WQMAD, winsorized QMAD ($c=1$); MLE, maximum likelihood estimator; LM, L-moment estimator; QM, quartile matching; MoM, method of moments.\label{tab:ais}}",
 r"\newcolumntype{C}{>{\centering\arraybackslash}X}", r"\begin{tabularx}{\textwidth}{lCCCCC}", r"\toprule",
 r"Method & $\hat m$ & $\hat s$ & $\hat\alpha$ & $\hat\delta$ & Sign agreement (\%)\\", r"\midrule"]
for k in KEYS:
    c = C[k]; b = BS[k]
    T.append(f"{k} & {c[0]:.2f} ({b[0]:.2f}) & {c[1]:.2f} ({b[1]:.2f}) & {fa(c[2])} & {c[3]:.3f} ({b[2]:.3f}) & {100-b[4]:.1f}" + r"\\")
T += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
F["%%TABLE_AIS%%"] = "\n".join(T)
T = [r"\begin{table}[H]",
 r"\caption{Sensitivity of the BMI fits to decimal-point errors. For each scenario: the estimated $\delta$ and the $W_1$ distance between the fit on the corrupted data and the same method's fit on the original data (Table~\ref{tab:ais}). The original fits all have $\hat\delta>0$. Abbreviations: QMAD, quartile-MAD estimator; WQMAD, winsorized QMAD ($c=1$); MLE, maximum likelihood estimator; LM, L-moment estimator; QM, quartile matching; MoM, method of moments.\label{tab:aissens}}",
 r"\newcolumntype{C}{>{\centering\arraybackslash}X}", r"\begin{tabularx}{\textwidth}{lCCCCCC}", r"\toprule",
 r" & \multicolumn{2}{c}{\textbf{2 values $\div10$}} & \multicolumn{2}{c}{\textbf{4 values $\div10$}} & \multicolumn{2}{c}{\textbf{2 values $\times10$}}\\",
 r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
 r"Method & $\hat\delta$ & $W_1$ & $\hat\delta$ & $W_1$ & $\hat\delta$ & $W_1$\\", r"\midrule"]
for k in KEYS:
    T.append(k + " & " + " & ".join(f"{d[k]['est'][3]:.3f} & {d[k]['w1']:.3f}" for d in (s1, s2, s3)) + r"\\")
T += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
F["%%TABLE_AISSENS%%"] = "\n".join(T)
F["%%AIS_LM_S1%%"] = f"{s1['LM']['est'][3]:.2f}"; F["%%AIS_LM_S2%%"] = f"{s2['LM']['est'][3]:.2f}"
F["%%AIS_LM_S3_S%%"] = f"{s3['LM']['est'][1]:.2f}"; F["%%AIS_LM_FLIP_C%%"] = f"{100-BS['LM'][4]:.1f}"
F["%%AIS_G%%"] = f"{o['GK_sample'][0]:.3f}"; F["%%AIS_K%%"] = f"{o['GK_sample'][1]:.3f}"
F["%%AIS_GCI%%"] = "$[%.3f,\\,%.3f]$" % tuple(o["GK_ci"][0]); F["%%AIS_KCI%%"] = "$[%.3f,\\,%.3f]$" % tuple(o["GK_ci"][1])
gf = o["GK_fit"]; F["%%AIS_KFIT%%"] = "$%.3f$--$%.3f$" % (min(v[1] for v in gf.values()), max(v[1] for v in gf.values()))
json.dump(F, open("results/fill_dyn.json", "w"), indent=1)
print(len(F))
