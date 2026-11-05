import json
T = {}
T["%%RMAIN%%"] = "$500$"; T["%%RDIST%%"] = "$1000$"
T["%%ABSTRACT_SIM%%"] = r"""In a simulation study with outliers in the left, right or both tails, the
proposed estimator fits the true distribution better than ML in all $144$
contaminated configurations, measured by the Wasserstein distance. ML reverses
the sign of the skewness when as little as $1\%$ of the data are outliers in the
tail opposite to the skewness. The winsorized version is the most accurate
estimator in most contaminated configurations. It is matched by quartile matching
for the scale, and it is clearly better than quartile matching for the shape. An
analysis of the body mass index of $202$ Australian athletes, with and without
artificial decimal-point errors, illustrates these findings. On clean data ML
remains the most efficient method, and the efficiency costs of the
quantile-based methods are quantified."""
T["%%SIM_RESULTS%%"] = r"""\subsection{Results under Tail Contamination}

Tables~\ref{tab:mainA} and~\ref{tab:mainB} give the results for $n=500$ and
$\alpha=3$, and Figure~\ref{fig:eps} shows the same configurations graphically.
Appendix~\ref{app:tables} reports all sample sizes and shapes. The main findings
are as follows.

\emph{QMAD against ML.} With any contamination, QMAD has a smaller mean $W_1$
distance than ML in all $144$ contaminated configurations ($4$ sample sizes $\times$
$3$ shapes $\times$ $4$ fractions $\times$ $3$ schemes). It has a smaller RMSE of
$\hat s$ in $137$ of them. At $\varepsilon=0.05$ the RMSE of $\hat s$ is $0.44$, $0.33$ and
$0.39$ for QMAD under schemes R, L and B, against $0.76$, $0.77$ and $0.56$ for
ML. For the shape the picture depends on the side of the outliers. With outliers in
the left tail, ML gets the wrong sign of $\alpha$ in $97\%$ of samples already
at $\varepsilon=0.01$, and in all samples from $\varepsilon=0.025$ on. QMAD's
wrong-sign rate rises only from $2\%$ (clean) to $7\%$ ($\varepsilon=0.05$). QMAD has
the smaller RMSE of $\hat\delta$ in all $48$ left-tail and $44$ of the $48$ two-tail
configurations. With outliers only in the right tail, the side of the skewness,
ML's $\hat\delta$ is pushed towards $1$. Since the true $\delta=0.949$ is close to $1$,
this produces a small RMSE of $\hat\delta$, and ML is better than QMAD for $\delta$ in
$45$ of the $48$ right-tail configurations. Even there, ML's fitted
\emph{distribution} is worse, because its scale is inflated.

\emph{L-moments.} LM is nearly as accurate as ML on clean data. At $n=500$ and
$\alpha=3$ its RMSE of $\hat s$ is $0.048$, against $0.046$ for ML. Under
contamination it behaves like a damped version of ML. Its mean $W_1$ is smaller than
ML's in all $144$ contaminated configurations, but larger than QMAD's in $140$ of
them, and its RMSE of $\hat s$ is larger than QMAD's in $137$. Like ML, it reverses
the skewness when the outliers lie in the short tail. At $n=500$ and $\alpha=3$, $1\%$
of left-tail outliers produce the wrong sign in $39\%$ of samples, and $5\%$ in all
samples. This agrees with Proposition~\ref{prop:if}(e). The L-skewness
numerator is pulled linearly by the outliers, whereas the QMAD skewness contrast
is bounded.

\emph{The benchmarks.} MoM behaves like ML. Its mean $W_1$ is worse than ML's in every
contaminated configuration. QM is the most robust method for the scale under heavy
contamination, but its shape estimate is the least precise of the quartile
methods (RMSE of $\hat\delta$ $0.25$--$0.48$ at $\varepsilon=0.05$).

\emph{WQMAD.} The winsorized estimator has the smallest mean $W_1$ in $98$ of the
$144$ contaminated configurations; QM has it in the other $46$. WQMAD also has the
smallest RMSE of $\hat s$ in $106$ configurations and of $\hat\delta$ in $90$. At
$\varepsilon=0.05$ and $\alpha=3$ its RMSE of $\hat\delta$ is $0.035$, $0.26$ and $0.035$
under R, L and B, against $0.25$, $0.48$ and $0.29$ for QM. Its weak spot is heavy
contamination on the side opposite to the skewness. With $10\%$ left-tail outliers,
WQMAD reverses the sign of $\alpha$ in $47\%$ of samples at $\alpha=3$, more often than
QMAD ($15\%$) and QM ($18\%$). For $\alpha=1$, whose skewness is weak, this happens
from $\varepsilon=0.025$ on. Every truncated outlier contributes the maximal amount
$\hat t$ to $\tilde d_1$ but little to $\tilde d_3$. Once the outliers make up a
sizeable fraction, they outweigh the small skewness contrast.

\emph{Clean data.} Without contamination, ML, LM and MoM are the most accurate, and the
ordering follows Table~\ref{tab:avar}. For example, at $n=500$ and $\alpha=3$ the RMSE
of $\hat s$ is $0.046$ (ML), $0.048$ (LM), $0.055$ (MoM), $0.064$ (WQMAD), $0.104$
(QMAD) and $0.108$ (QM). The corresponding mean $W_1$ values are $0.028$, $0.029$,
$0.031$, $0.035$, $0.044$ and $0.049$. In terms of the fitted distribution, the loss from using WQMAD
instead of ML on clean data is therefore small. For $\alpha=1$, all methods find the
sign of $\alpha$ hard to estimate. At $n=500$ the wrong-sign rate is $13\%$ for ML,
LM and MoM, $26\%$ for WQMAD and $37\%$ for QMAD and QM.

%%TABLE_MAINA%%

%%TABLE_MAINB%%

\begin{figure}[H]
\includegraphics[width=13.5cm]{figures/fig_eps_a3.pdf}
\caption{RMSE of $\hat s$, RMSE of $\hat\delta$ and mean $W_1$ as functions of the
contamination fraction, for $n=500$, $\alpha=3$ and outliers at $k=5$ IQRs beyond the
quartiles.\label{fig:eps}}
\end{figure}

\subsection{Effect of the Outlier Distance}

Figure~\ref{fig:dist} varies the distance $k$ of the outliers at $\varepsilon=0.05$,
$n=500$ and $\alpha=3$. The errors of ML and MoM grow linearly in $k$, with slopes
that reflect their quadratic and cubic influence functions. The ML scale error
reaches $3.5$ at $k=20$, and for every left-tail distance $k\ge2$ ML reverses
the skewness in at least $99.6\%$ of samples. QMAD's errors also grow linearly in $k$, but with a slope about
$3$ times smaller than ML's (Remark~\ref{rem:rates}). The errors of WQMAD and QM
do not grow with $k$ at all. LM's errors grow linearly too. Its scale RMSE at $k=20$
is $2.61$, against $1.47$ for QMAD, and with left-tail outliers it reverses the
skewness in $89\%$ of samples at $k=2$ and in all samples for $k\ge5$. For $k\le1$ the outliers are not yet outlying, and all
methods perform similarly.

\begin{figure}[H]
\includegraphics[width=13.5cm]{figures/fig_dist.pdf}
\caption{Effect of the outlier distance $k$ (in IQRs beyond the quartiles) for
$\varepsilon=0.05$, $n=500$ and $\alpha=3$, with $1000$ replications. Top: right-tail
outliers; bottom: left-tail outliers.\label{fig:dist}}
\end{figure}

\subsection{Effect of Contamination on Shape Measures}\label{sec:sim:shapes}

Gross errors affect not only parameter estimates but also the descriptive shape
measures used to choose a model. Table~\ref{tab:shapes} compares the sample Pearson
skewness and kurtosis, the L-skewness and L-kurtosis, the MAD skewness and kurtosis
$G$ and $K$ of Section~\ref{sec:madshape}, and the two quartile-based skewness
measures $\kappa$ and $b$. Both are computed under the contamination schemes of
Section~\ref{sec:sim}, for $\alpha=3$ and $n=500$.

For kurtosis the picture is clear. Pearson's $\beta_2$ quadruples, from $3.5$ to $14.0$,
when only $1\%$ of the observations are outliers. The L-kurtosis rises from $0.124$
to between $0.29$ and $0.33$ at $\varepsilon=0.05$, an increase by a factor of $2.4$ to
$2.6$. The MAD kurtosis $K$ rises from $0.593$ to at most $0.704$, an increase of
$19\%$. It is bounded by $1$ by construction~\cite{PinskyKundu2026}, and it is
by far the most resilient of the three.

For skewness, the side of the outliers matters. With $5\%$ of left-tail outliers,
the Pearson skewness, the L-skewness and the MAD skewness $G$ all change sign. They
reach $-2.40$, $-0.222$ and $-0.164$, against population values of $0.667$, $0.119$
and $0.165$. Pearson's $\gamma_1$ is already negative at $\varepsilon=0.01$. $G$ is
affected because it uses the mean. The two quartile-anchored measures, $\kappa$ and
Bowley's $b$, keep the correct sign and move only moderately. This is the
mechanism behind the robustness of QMAD, and it is why the estimator is built on
the quartiles rather than on the mean or the sub-means
(Section~\ref{sec:madshape}).

%%TABLE_SHAPES%%

\subsection{The Plug-in Estimator}\label{sec:sim:plugin}

Table~\ref{tab:plug} shows the plug-in estimator~\eqref{eq:plugin} on clean data.
The root is chosen by the quartile-consistency rule of Section~\ref{sec:plugin}.
Even at $n=2000$, the plug-in estimator gets the sign of $\alpha$ wrong in $13$--$16\%$
of samples for $\alpha\in\{2,3\}$, where QMAD and ML practically never do. Its RMSE
of $\hat\delta$ does not decrease appreciably with $n$. This confirms that the
spurious roots of Proposition~\ref{prop:plugin} matter in practice.

%%TABLE_PLUG%%"""
json.dump(T, open("results/fill_text.json", "w"), indent=1)
