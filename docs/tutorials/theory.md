
# Targeted MD Simulations

For accelerated sampling of rare events, biased simulation methods have been developed.  
One particular method is targeted MD developed Schlitter et al. \cite{Schlitter94}. Here, a constraint is applied to a subset of atoms to move it towards a target conformation along a predescribed one-dimensional path in conformational space with constant velocity $v_c$ from $s_0 (t=0)$ to $s (t=t_f)$  along the pulling coordinate $s(t)=s_0+v_{c}t$. The corresponding constraint function is

$$	\Phi(s(t)) = \left( s(t) -(s_0 + v_0 t) \right) = 0. $$

The necessary constraint force $f_c$, 

$$	f_c = \lambda \frac{\mathrm{d} \Phi(s(t))}{\mathrm{d}s} = \lambda $$

added through a Lagrangian-multiplier λ to the equation of motion, is stored to compute the work
$W(s)=\int_{s_0}^{s} f_c(s')\; ds' \geq\Delta G$ needed to move the atom, or subset of atoms, along the pulling
coordinate.


# Jarzynski's Equality

An approach to estimate ∆G for finite switching times was proposed by Jarzynski [4]. It is possible
to calculate the free energy directly from a thermostated system in equilibrium that is driven away
from equilibrium (e.g. in non-equilibrium pulling trajectories) along a pulling coordinate $s(t)$ via

$$	e^{-\beta\Delta G(s)} = \left< e^{-\beta W(s)} \right>_\text{N}. $$

Here $\Delta G$ denotes the Gibbs free energy,
$\left< \cdot  \right>_\text{N}$ denotes an ensemble average over independent pulling trajectories starting from an initial Boltzmann distribution. The original formulation of Eq. (3) is for an NVT ensemble predicting the Helmholtz free energy $\Delta F$ instead of $\Delta G$ and was later generalized for NPT ensembles \cite{Park04}. 

Equation (3) allows to directly calculate the free energy profile from biased simulations. However, the exponential average in Eq. (3) highly relies on a sufficient sampling of the small values of the work distribution, which occur only rarely \cite{Park04}. Thus, in practice a large number of non-equilibrium trajectories is required, due to the slow and erratic convergence of the exponential average. This problem can be approached by a cumulant expansion of Eq. (3):

$$	\Delta G(s) = \left< W(s) \right>_\text{N} -\frac{\beta}{2}\left<\delta W(s)^2\right>_\text{N} + h.o.t. 
$$

with $\delta W(s) = W(s) -\left< W(s) \right>_\text{N}$. 
The expansion in Eq. (4) can be truncated at second order if the work distribution is Gaussian. This allows to formulate the dissipated work 

$$	W_{\text{diss}}(s) = \frac{\beta}{2}\left< \delta W(s)^2 \right>_\text{N}
$$

in terms of the variance of the work.


# Dissipation-Corrected TMD

A method developed by Wolf and Stock \cite{Wolf18} combines Langevin dynamics with the second order cumulant expansion given in Eq. (\ref{eq:Jary2ndcum}). It evaluates $\Delta G(s)$ as well as a non-equilibrium friction coefficient $\Gamma_\text{N}$ directly from TMD trajectories and is called dissipation-corrected TMD (dcTMD).

As stated before, in TMD simulations a constraint force $f_c$ is applied to the system. We assume that $f_c$ does not change the Langevin fields $f$, $\Gamma$, and $K$ and make the ansatz, that $f_c$ can be simply included as an additive term in the memory free Langevin equation \cite{Marquart83}. This yields:

$$	m \ddot{s}(t) = - \frac{\text{d}G}{\text{d}s} - \Gamma(s)\dot{s} + K(s)\xi(t) + f_c(t).
$$

To keep the velocity constant, $f_c(t)$ needs to counter al occurring forces, thus

$$	m \ddot{s}(t)=0. $$

Equation (7), however, only holds if the following requirements on the constraint force $f_c$ are fulfilled: first, $v_c$ is kept constant against the drag of the friction $\Gamma(s)$ and the acceleration due to the free energy gradient; second, $f_c$ counterbalances the effect of the stochastic force $K(s)\xi(t)$.

The first requirement is fulfilled if $v_c$ is slow compared to the bath fluctuations. Then the pulling can be considered as a slow adiabatic change \cite{Servantie03} and the system is virtually in equilibrium at every point along $s$. 

The second requirement is harder to satisfy, since $f_c$ and the bath fluctuations necessarily occur in the same time scales. When considering an ensemble of TMD trajectories instead of a single trajectory the stochastic term cancels out by definition. Thus, performing an ensemble average of Eq. (6) over a set of trajectories one obtains:

$$
	\frac{\text{d}G}{\text{d}s} = -\Gamma(s)v_c + \left<\delta f_c(s)\right>_\text{N}
$$   

by making use of $\left< \xi \right>_\text{N}=0$ and $\left< \dot{s} \right>_\text{N}=v_c$.
Integrating on both sides from the initial state $s_0$ to final state $s$ gives:

$$
\begin{align}
	\Delta G(s) &= - v_c\int_{s_0}^{s} \text{d}s' \Gamma(s') + \int_{s_0}^{s} \text{d}s' \left<f_c(s')\right>_\text{N} \\
				&=  - W_{diss}(s) + \left<W(s)\right>_\text{N} 
\end{align} 
$$

where the first term describes the dissipated work of the process in terms of the friction $\Gamma$
and the second term corresponds to the average external work carried out on the system.

With this, the dissipated work

$$
\begin{align}
	W_\text{diss} &= \frac{\beta}{2} \left<\left(\int_{s_0}^{s} \delta f_c (s) \text{d}s\right)^2\right>_\text{N}\\
	&= \beta \int_{s_0}^{s} \text{d}s_2 \int_{s_0}^{s_2} \text{d}s_1 \left<\delta f_c (s_2) \delta f_c (s_1)\right>_\text{N}\\
\end{align}
$$

with the force $\delta f_c = f_c - \left<f_c\right>$, can be 
compared with Eq. \ref{eq:dGGammaWdiss}, suggesting that the NEQ friction term $\Gamma$ can be expressed by

$$
\begin{align}
	\Gamma &= \frac{\beta}{v_c} \int_{s_0}^{s} \text{d}s_1 \left<\delta f_c (s) \delta f_c (s_1)\right>_\text{N} \\
	&= \beta \int_{0}^{t} \text{d}\tau \left<\delta f_c(t)\delta f_c(t-\tau)\right>_\text{N},
\end{align}
$$

where a substitution of variables was performed using $s_1 = s_0 + v_c t'$ and $\tau = t - t'$. Thus, a non-equilibrium friction correction is obtained 

1. under the assumption of a Gaussian work distribution, which determines the friction arising from a set of TMD simulations

2. pulled with constant mean velocity $v_c$ along a reaction coordinate $s$.

Equation \ref{eq:GammaNEQ} has a similar form to the well known friction expression in equilibrium \cite{Zwanzig01}

$$
	\Gamma(s)_{\text{EQ}} = \beta \int_{0}^{\infty} \left<\delta f_c(t)\delta f_c(0)\right>_{\text{EQ},s} \text{d}t,
$$


# Classes calculating the dcTMD quantities.

There are two ways to calculate the dcTMD quantities $\Delta G$ and $\Gamma$. 

## WorkEstimator

One is via Eq. (9) and (10). This is implemented in the WorkEstimator class.

$$ 
\begin{align}
\Delta G(s) &=  - W_{diss}(s) + \left<W(s)\right>_\text{N},\\
\Gamma &= \frac{1}{v_c}\frac{d~W_{diss}(s)}{ds}
\end{align}
$$

Another way is via the force autocorrelation function. This is implemented in the ForceEstimator class using equations 

## Force Estimator 

Class for performing dcTMD analysis via force $f_c(t)$ auto correlation:

$$
\begin{align}
\Gamma &= \beta \int_{0}^{t} \text{d}\tau \left<\delta f_c(t)\delta f_c(t-\tau)\right>_\text{N},\\
\Delta G(s) &= - v_c\int_{s_0}^{s} \text{d}s' \Gamma(s') + \int_{s_0}^{s} \text{d}s' \left<f_c(s')\right>_\text{N}
\end{align}
$$

