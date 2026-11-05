Code for "Robust Estimation of the Skew-Normal Distribution from Mean Absolute
Deviations about the Quartiles". Python 3, numpy, scipy, pandas, matplotlib.

snmad.py       closed-form E|X-a| (H), its derivatives, ML fit, plug-in system
qmad.py        QMAD estimator (model quartiles, lookup tables), plug-in fit,
               quartile matching (qm_fit), method of moments (mom_fit), ML wrapper
lmom.py        L-moment estimator (lm_fit); ltheory.py its influence function
madshape.py    MAD skewness G and MAD kurtosis K (closed form and sample versions)
wqmad.py       winsorized QMAD (class WQMAD, c = 1 by default in the paper)
asym.py, theory.py, theory2.py, wtheory.py
               influence functions, asymptotic variances, bias functionals
sim.py         Monte Carlo driver (QMAD, ML; plug-in study); seeds fixed per cell
sim_extra.py, sim_extra2.py
               QM, MoM and WQMAD on the identical regenerated samples
analyze.py     performance measures; fig_*.py figures; ais*.py AIS example
ais_sn.csv     AIS data (from the R package sn)
Paths in the scripts point to the authors' working folders; adjust as needed.
