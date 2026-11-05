import json
BIAS = r"""The left half of the figure shows the failure mode that motivates this paper.
A $5\%$ point mass at any $x\le-1.25$ makes the ML shape functional negative. The
true distribution is clearly right-skewed with $\delta=0.949$, yet for every
$x\le-2.5$ ML fits a left-skewed distribution with $\delta$ between $-0.99$ and $-0.90$,
and its location is displaced by $1.5$ to $2.1$ scale units. The point $x=-1.25$ lies only
$1.8$ interquartile ranges below $q_1$, so the outliers involved are not
exceptionally remote. The QMAD shape functional stays within $0.10$ of the truth for every
point mass outside the interval $[0.85,1.18]$, and its scale bias at $x=-10$ is
$0.75$, against $2.10$ for ML. For outliers in the right tail, both
estimators keep the correct sign. QMAD has the smaller bias for scale ($1.18$
against $2.71$ at $x=15$) and for shape. ML has the smaller bias for location
($-0.47$ against $-0.71$ at $x=15$).

The interval $[0.85,1.18]$ lies just below the upper quartile $q_3=1.150$. A point
mass of $5\%$ there distorts the midmean--midhinge contrast enough to reverse the
sign of $\hat\kappa$. This is the price of identifying the shape from a small
contrast: $\kappa(3)=0.060$ only. This case is an \emph{inlier} cluster, not
tail contamination.

The two benchmarks bracket QMAD. MoM behaves like ML, only worse: it reverses
the skewness for the same left-tail outliers, and its scale bias reaches $4.26$ at
$x=15$. QM has bounded bias. For outliers with $x<-1.5$ or $x>2$, its biases never
exceed $0.05$ (location), $0.13$ (scale) and $0.03$ ($\delta$). For tail
contamination it is therefore the most robust of the four estimators. Like QMAD,
it reverses the sign for inlier clusters below $q_3$, and its interval is about
twice as wide: $[0.73,1.10]$ against $[0.98,1.15]$ for QMAD."""
json.dump({"%%BIAS_TEXT%%": BIAS}, open("results/fill_static.json", "w"))

DISC = r"""\textbf{What the method buys and what it costs.} The QMAD estimator replaces
the quadratic influence of the ML score with a linear one. The skewness
information comes from a bounded contrast, $\mathrm{IF}_{d_3}-\mathrm{IF}_{d_1}$,
normalised by the MAD about the median. This is why it does not reverse the sign
of the skewness when outliers fall in the short tail, and why its scale estimate
is less inflated by remote points (Remark~\ref{rem:rates}). The cost is efficiency
on clean data (Table~\ref{tab:avar}). This cost is moderate for $1\le\alpha\le5$ and
severe for strongly skewed populations. There $\kappa(\alpha)$ is close to
its ceiling, and the quartiles carry little information about $\alpha$.
Practitioners who are confident that their data are uncontaminated should use ML,
possibly with the bias corrections of~\cite{Sartori2006,AzzaliniArellanoValle2013}.
QMAD is the better choice when a few per cent of gross errors in the tails cannot
be excluded. It is also useful as a fast, deterministic starting value for ML, and as
a diagnostic: a large discrepancy between the QMAD and ML fits signals
observations that drive the likelihood.

\textbf{Limitations.} First, the QMAD influence function is unbounded, so the
estimator has breakdown point zero in the formal sense. A single observation
moved to infinity carries $\hat s$ with it, only more slowly than it carries
$\hat s_{\mathrm{ML}}$. Second, because $\kappa$ is small, a cluster of
inliers just below the upper quartile can reverse the sign of $\hat\alpha$
(Figure~\ref{fig:bias}). Third, near $\alpha=0$ the shape is identified only to
order $\alpha^3$, so small shapes are hard to estimate, as they are for ML. Large
samples are needed to estimate $\alpha\approx1$ reliably by either method.

\textbf{Extensions.} The construction is not specific to the quartiles or to the
normal kernel. With MADs about $k>3$ quantiles, the equations
$\hat d_j=s\,D_j(\alpha)$ are over-identified and can be combined by generalized
method of moments or by the method of simulated quantiles~\cite{DominicyVeredas2013}.
This should recover part of the efficiency lost on clean data. Replacing
$\hat d_j$ by a trimmed or winsorized MAD would bound the influence function at
the cost of a modified $D_j$. The skew-$t$ and other skew-symmetric families have
closed-form MAD functions as well~\cite{CompanionMAD}, so the same three equations
apply to them with $H$ replaced accordingly."""

CONC = r"""Equating empirical and theoretical mean absolute deviations about the quartiles
yields a simple estimator of the skew-normal parameters. The equations must be
formulated with the \emph{model} quartiles, not the sample quartiles. The model
quartiles are cheap to compute, and once they are used the location parameter
drops out of the MAD equations. The resulting system is triangular and has a unique
solution. The shape is determined by the robust skewness ratio
$\kappa=(\text{midhinge}-\text{midmean})/\text{MAD}_{\text{median}}$, the scale by the
MAD about the median, and the location by the median. The estimator is
asymptotically normal with a linearly growing influence function. Under
contamination in the tails, it estimates scale and shape more accurately than
maximum likelihood, and it never reverses the sign of the skewness in our
experiments with tail outliers. Maximum likelihood remains preferable for
clean data."""
d = json.load(open("results/fill_static.json")); d["%%DISCUSSION%%"] = DISC; d["%%CONCLUSIONS%%"] = CONC
json.dump(d, open("results/fill_static.json", "w"))
