import json
sh = json.load(open("results/shapes.json")); ai = json.load(open("results/ais.json"))["shapes"]
NAMES = [("g1", r"Pearson skewness $\gamma_1$"), ("t3", r"L-skewness $\tau_3$"), ("G", r"MAD skewness $G$"),
         ("kappa", r"Quartile-MAD skewness $\kappa$"), ("bowley", r"Bowley skewness $b$"),
         ("b2", r"Pearson kurtosis $\beta_2$"), ("t4", r"L-kurtosis $\tau_4$"), ("K", r"MAD kurtosis $K$")]
L = [r"\begin{table}[H]",
 r"\caption{Effect of tail contamination on sample shape measures: mean (standard deviation) over $1000$ samples of size $n=500$ from $\mathrm{SN}(0,1,3)$, with a fraction $\varepsilon$ of the observations replaced by outliers $5$ IQRs beyond the lower (L) or upper (R) quartile. The first column gives the population value.\label{tab:shapes}}",
 r"\begin{adjustwidth}{-\extralength}{0cm}", r"\newcolumntype{C}{>{\centering\arraybackslash}X}", r"\small",
 r"\begin{tabularx}{\fulllength}{lCCCCCC}", r"\toprule",
 r"Measure & Population & Clean & L, $\varepsilon=0.01$ & L, $\varepsilon=0.05$ & R, $\varepsilon=0.01$ & R, $\varepsilon=0.05$\\", r"\midrule"]
for k, lab in NAMES:
    if k == "b2": L.append(r"\midrule")
    k_ = 2 if k in ("g1", "b2") else 3
    L.append(lab + f" & {sh['pop'][k]:.{k_}f} & " + " & ".join(f"{sh[t][k][0]:.{k_}f} ({sh[t][k][1]:.{k_}f})" for t in ["clean", "L1", "L5", "R1", "R5"]) + r"\\")
L += [r"\bottomrule", r"\end{tabularx}", r"\end{adjustwidth}", r"\end{table}"]
T = [r"\begin{table}[H]",
 r"\caption{Sample shape measures of the AIS BMI data before and after the decimal-point errors of Table~\ref{tab:aissens}.\label{tab:aisshapes}}",
 r"\newcolumntype{C}{>{\centering\arraybackslash}X}", r"\begin{tabularx}{\textwidth}{lCCCC}", r"\toprule",
 r"Measure & Original & 2 values $\div10$ & 4 values $\div10$ & 2 values $\times10$\\", r"\midrule"]
for k, lab in NAMES:
    if k == "b2": T.append(r"\midrule")
    k_ = 2 if k in ("g1", "b2") else 3
    T.append(lab + " & " + " & ".join(f"{ai[s][k]:.{k_}f}" for s in ["orig", "2 errors /10", "4 errors /10", "2 errors x10"]) + r"\\")
T += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
json.dump({"%%TABLE_SHAPES%%": "\n".join(L), "%%TABLE_AISSHAPES%%": "\n".join(T)}, open("results/fill_shapes.json", "w"))
print("ok")
