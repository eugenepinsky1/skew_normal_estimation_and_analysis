import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
QC, MC = "#2a78d6", "#eb6834"
STY = {  # colour, linestyle, marker  (line style + marker = secondary encoding)
 "qmad":  dict(color="#2a78d6", ls="-",  marker="o", label="Q-MAD"),
 "wqmad": dict(color="#14946a", ls=(0, (6, 1.5)), marker="v", label="WQ-MAD"),
 "ml":    dict(color="#eb6834", ls="--", marker="s", label="MLE"),
 "qm":    dict(color="#a0522d", ls="-.", marker="^", label="QM"),
 "mom":   dict(color="#4a4a4a", ls=":",  marker="D", label="MoM"),
 "lm":    dict(color="#8e24aa", ls=(0, (3, 1, 1, 1, 1, 1)), marker="P", label="LM"),
}
ORDER = ["qmad", "wqmad", "ml", "lm", "qm", "mom"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6, "axes.edgecolor": "#888",
    "lines.linewidth": 1.6, "legend.frameon": False, "savefig.bbox": "tight", "font.family": "serif"})
