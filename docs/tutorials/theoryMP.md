
# Targeted Molecular Dynamics

In order to study the dynamics of biomolecular processes by computer simulation, biased simulation methods have been developed that can accelerate the sampling of otherwise rare events. One such method is Targeted Molecular Dynamics (TMD) developed by [Schlitter et al. 1994](https://doi.org/10.1016/0263-7855(94)80072-3).
In the specific implementation of TMD of this software package, a constraint is applied on the distance $s$ between two subsets of atoms, moving them with constant velocity $v$ from an initial state $s_0$ along the pulling coordinate towards a target conformation. This constraint,

$$ \Phi(s,t) =  s -(s_0 + v t)  = 0 \, , $$

modifies the equation of motion of the $K$ atoms $\mathbf{r}=(\mathbf{r}_1,\dots,\mathbf{r}_K)$ bound by the potential $U(\mathbf{r})$ by an additional constraint force term

$$  m_i \ddot{\symbf{r}}_i = - \frac{\partial U}{\partial \mathbf{r}} + f \frac{\partial \Phi}{\partial \mathbf{r}} \; , $$

that is, the constraint force $f$ on $s$ is given by the Lagrange multiplier which in practice is calculated via the SHAKE algorithm ([Ryckaert et al. 1977](https://doi.org/10.1016/0021-9991(77)90098-5)). The constraint forces are then used to compute the work

$$ W(s)=\int_{s_0}^{s}  \mathrm{d} s'  \; f(s') $$

done to move the atom, or subset of atoms, along the pulling coordinate. Note that in the follwing we use the notations $f(s) \equiv f(t(s))$ interchangeable. 



# Jarzynski's Equality

In order to obtain equilibrium estimates, such as the free energy profile $\Delta G$ along the pulling coordinate $s$, from the nonequilibrium simulations described above, we use the method proposed by [Jarzynski 2004](https://doi.org/10.1088/1742-5468/2004/09/P09005). While for finite pulling times, the average work $\left<W(s)\right>$ is always higher that $\Delta G = G(s)-G(s_0)$, it is also possible to give an exact relation between the work and free energy, namely

$$
\begin{align}
\mathrm{e}^{-\beta\Delta G(s)} = \left< \mathrm{e}^{-\beta W(s)} \right>. 
\end{align}
$$

$\left< A  \right>$ denotes an ensemble average of some function of the atoms positions $\mathbf{r}$ and momenta $\mathbf{p}$ over independent pulling trajectories starting from an initial Boltzmann distribution at fixed $s$, i.e.,

$$
\begin{align*}
	\left< A  \right> &= \frac{1}{\int  \mathrm{d} \mathbf{r}_0  \mathrm{d} \mathbf{p}_0 \,  \mathrm{e}^{-\beta H(\mathbf{r}_0,\mathbf{p}_0,0)}} \int  \mathrm{d} \mathbf{r}_0  \mathrm{d} \mathbf{p}_0 \, A(\mathbf{r}_0,\mathbf{p}_0,t) \; \mathrm{e}^{-\beta H(\mathbf{r}_0,\mathbf{p}_0,0)}\\
	&\approx \frac{1}{\mathtt{N}} \sum_{\mathtt{i}=1}^{\mathtt{N}} A_\mathtt{i} \; ,
\end{align*} 
$$

where $\beta^{-1} = \mathrm{k}_\mathrm{B} T$ is the inverse temperature and $H$ the time-dependent Hamiltonian. $\Delta G$ denotes the Gibbs free energy. While the original formulation of Eq. (1) was proved for an NVT ensemble predicting the Helmholtz free energy $\Delta F$, it was later generalized for NPT ensembles (see [Park et al. 2004](https://aip.scitation.org/doi/abs/10.1063/1.1651473)). 

For a finite sample of $\mathtt{N}$ pulling trajectories, the bottom equation, a sample average, can be used. The exponential average, however, highly relies on a sufficient sampling of the small values of the work distribution, which occur only rarely. Thus, in practice a large $\mathtt{N}$ of non-equilibrium trajectories is required due to the slow convergence of the otherwise biased exponential average. This problem can be approached by a cumulant expansion:

$$
\begin{align}
\Delta G(s) = \left< W(s) \right> -\frac{\beta}{2}\left<\delta W(s)^2\right> + \dots
\end{align} 
$$

with $\delta W = W -\left< W\right>$. The expansion can be truncated at second order if the work distribution is Gaussian. In this case, the average dissipated work $W_{\text{diss}}$ is given entirely by the variance, i.e.,

$$	\left<W_{\text{diss}}(s)\right> = \frac{\beta}{2}\left< \delta W(s)^2 \right> \; .
$$


# Dissipation-Corrected TMD

A method developed by [Wolf, Stock 2018](https://doi.org/10.1021/acs.jctc.8b00835) called dissipation-corrected TMD (dcTMD) combines Langevin dynamics (see, e.g., [Zwanzig 2001](https://global.oup.com/academic/product/nonequilibrium-statistical-mechanics-9780195140187?q=zwanzig&lang=en&cc=de)) with the second order cumulant expansion given in Eq. (2). It evaluates $\Delta G$ as well as a non-equilibrium friction coefficient ${\it{\Gamma}}$ directly from TMD simulations.

Here, the memory-free Langevin equation models the dynamics projected to a small set of coordinates, here the motion of the distance $s$, in terms of the mean force $\frac{\partial G}{\partial s}$, friction ${\it{\Gamma}}$ and fluctuations $\xi$, as well as the applied constraint force, as

$$	
\begin{align}
m \ddot{s}= - \frac{\partial G}{\partial s} - {\it{\Gamma}}(s)\,\dot{s} + \xi(t) + f(t) \quad ( = 0 ).
\end{align}
$$

$m$ is the reduced mass of the two atom subsets. This approach assumes that $f_c$ does not change the Langevin fields $G$, ${\it{\Gamma}}$, and the noise $\xi$ and makes the ansatz, that $f$ can be simply included as an additive term. 

When considering an ensemble of TMD trajectories instead of a single one, the stochastic term in Eq. (3) cancels out, $\left< \xi \right>=0$. Together with $\dot{s}=v$, integrating on both sides from the initial state $s_0$ to final state $s$ gives

$$
\begin{align}
	\Delta G(s) &= \int_{s_0}^{s}  \text{d}s' \left<f(s')\right> - v \int_{s_0}^{s}  \text{d}s' \; {\it{\Gamma}}(s')  
\end{align}
$$

where the first term corresponds to the average external work done on the system and the second term describes the dissipated work of the process in terms of the friction ${\it{\Gamma}}$. Evaluating the latter,

$$
\begin{align*}
	\left<W_\text{diss}\right> &= \frac{\beta}{2} \left<\left(\int_{s_0}^{s} \mathrm{d}s^\prime \delta f(s') \right)^2\right>\\
	&= \beta \int_{s_0}^{s}  \text{d}s^\prime \int_{s_0}^{s^\prime}  \text{d}s^{\prime\prime} \left<\delta f(s^{\prime\prime}) \delta f (s^{\prime})\right> \; ,\\
\end{align*}
$$

provides a connection between the dissipated work and the friction,

$$
\begin{align}
 {\it{\Gamma}}(s) &= \frac{1}{v}\frac{\mathrm{d} W_\text{diss}}{\mathrm{d} s} \\
&= \beta \int_{0}^{t(s)} \text{d}\tau \left<\delta f(t(s))\delta f(\tau)\right>  \; ,
\end{align}
$$

assuming  a Gaussian work distribution.




# Classes calculating the dcTMD quantities


### WorkEstimator

One way to analyze nonequilibrium force time traces from constraint pulling simulations is by calculating the work first and then estimating the free energy  $\Delta G$ via Eq. (2), as well as the friction ${\it{\Gamma}}$ via Eq. (5). This approach is implemented in the WorkEstimator class.


### ForceEstimator 

Another option to directly analyze the constraint force time traces, calculate $\Delta G$ on-the-fly by via Eq. (4) and the friction ${\it{\Gamma}}$ via Eq. (6).  This approach is implemented in the ForceEstimator class

We advise the use of the WorkEstimator class, which is computationally less demanding, since it does not require the full resolution of the force time traces.

