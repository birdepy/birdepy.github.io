========================
Transition Probabilities
========================
Many of the :ref:`general models implemented in Birdepy <Birth-and-death Processes>` have
no explicit expression for their transition probabilities
:math:`p_{i,j}(t) := \mathbb P(Z(t)=j\mid Z(0)=i)`, :math:`i,j\in \mathcal S`, :math:`t\geq 0`. 
However, it is often desirable to compute these transition probabilities since they allow for a deeper understanding of the future (random) evolution of the process. 
Practically speaking, this may underpin some performance analysis of a system which is being modelled by the process. 
Additionally, transition probabilities allow likelihood functions of discretely-observed sample paths to be evaluated, which opens up the possibility of parameter estimation (as discussed in the next section). 
This page outlines nine methods of approximating :math:`p_{i,j}(t)`, as implemented in :ref:`birdepy.probability()`.  
Three of these methods are matrix-based and rely on state space truncation, another method uses Laplace transforms, and the remaining four methods use approximation models. 

To obtain transition probabilities in BirDePy simply call :ref:`birdepy.probability()` (after importing BirDePy). For example, to obtain an array containing the probabilities of transitioning from :math:`Z(0)=20` to :math:`Z(1)=25` for the 'Verhulst' model with :math:`\boldsymbol \theta = (0.8, 0.4, 0.01, 0.001)`:: 
	 
   import birdepy as bd
   bd.probability(20, 25, 1.0, [0.8, 0.4, 0.01, 0.001], model='Verhulst', method='expm')

Which returns::

   array([[0.08189476]])

The method 'expm' can be replaced by one of the nine choices described on this page. 
Note that output is an array since it is possible to input initial and final population sizes as lists. 


Matrix exponential
^^^^^^^^^^^^^^^^^^
The Kolmogorov forward equation for CTMCs provides the foundational property that the collection of transition probability functions :math:`p_{i,j}(t)` satisfy, for :math:`t\ge0`, the system of first order ordinary differential equations

	.. math::

		\frac{ d p_{i,0}(t)}{d t} &= \mu_1 p_{i,1}(t) - \lambda_0p_{i,0}(t),\\
		\frac{ d p_{i,j}(t)}{d t} &= \lambda_{j-1}p_{i,j-1}(t) + \mu_{j+1}p_{i,j+1}(t)- (\lambda_j+\mu_j)p_{i,j}(t), \quad j\ge 1,\\

with :math:`p_{i,i}(0) = 1` and :math:`p_{i,j}(0)=0` for :math:`j\ne i`. 
Let :math:`Q` be the *generator* of :math:`Z`, that is, the square matrix with diagonal entries :math:`q_{i,i}=-(\lambda_i+\mu_i)`, upper diagonal entries :math:`q_{i,i+1}=\lambda_i`, and lower diagonal entries :math:`q_{i-1,i}=\mu_i` (and zeros elsewhere). 
Also collect the transition probabilities into a matrix :math:`P(t) = \big(p_{i,j}(t); i,j\in\mathcal S\big)`. 
Then the equation above can be written in matrix form as :math:`\frac{d}{d t} P(t) = P(t)Q`. 
When :math:`Q` is finite it immediately follows that 

	.. math::

		P(t) = \lim_{k\to\infty}\sum_{n=0}^k \frac{1}{n!}(Qt)^{n} =: \exp(Qt).

See :ref:`birdepy.probability(method='expm')` for the implementation of this method in BirDePy. 


Uniformization
^^^^^^^^^^^^^^
Due to the special nature of the generator :math:`Q` (specifically that it has non-negative off-diagonal entries and row-sums equal to 0), an alternative procedure known as *uniformization*, introduced by [Jensen1953]_, is available for computing :math:`P(t)`.  
Consider a discrete-time Markov chain (DTMC) with probability transition matrix :math:`A:=Q/a+I`, where :math:`a:=\max_z |\lambda_z+\mu_z|` when this exists (for example when :math:`|\mathcal S|<\infty`), and where :math:`I` denotes the identity matrix of appropriate size.
Suppose that this DTMC transitions at the event times of a Poisson process :math:`(N(t),~t\ge0)` with rate :math:`a`. 
Using this construction, we can show that 

	.. math::

		\mathbb P(Z(t)=j\mid Z(0)=i, N(t)=n)=(A^n)_{i,j}.

Therefore, conditioning on :math:`N(t)` provides

	.. math::

		P(t) = \lim_{k\to\infty}\sum_{n=0}^k \frac{(ta)^ne^{-ta}}{n!} A^n. 

In addition to potentially requiring truncation of the state space for many models of interest, another source of approximation error for this method is that a finite :math:`k` must be chosen in the above infinite series. 

This method is discussed at length in [Grassman1977]_ and [vanDijkEtAl2018]_. 

See :ref:`birdepy.probability(method='uniform')` for the implementation of this method in BirDePy. 

Erlangization
^^^^^^^^^^^^^
Using probabilistic arguments, *Erlangization* is yet another matrix-based method for computing :math:`P(t)`. 
Let :math:`T` be an exponentially distributed random variable with mean :math:`\eta^{-1}>0`. 
Define :math:`r^{(\eta)}_{i,j} := \mathbb P(Z(T) = j\mid Z(0)= i)` and collect these quantities into  :math:`R(\eta) = \big(r^{(\eta)}_{i,j},~i,j\in \mathcal S)`. 
The matrix :math:`R(\eta)` can be viewed as the one-step transition probability matrix of a DTMC embedded in :math:`Z` at the epochs of a Poisson process with rate :math:`\eta`. 
By conditioning on the first transition of :math:`Z` and the expiry of the time :math:`T`, we obtain the recursive expressions 

	.. math:: 

		r^{(\eta)}_{i,j}= \frac{\mu_ir^{(\eta)}_{i-1,j} + \lambda_ir^{(\eta)}_{i+1,j} + \eta 1_{\{i=j\}}}{\lambda_i + \mu_i + \eta},\quad i,j\in\mathcal S. 

This system of equations can be written in matrix form in terms of :math:`R(\eta)`, whose solution is given in terms of the generator :math:`Q` by

	.. math:: 

		R(\eta) = \eta\big(\eta I-Q)\big)^{-1}.

Therefore, if we let :math:`\eta:=k/t` and :math:`S_{k,t}` be an Erlang distributed random variable with rate parameter :math:`k/t` and shape parameter :math:`k`, then the :math:`(i,j)`th entry of the matrix :math:`R(k/t)^k` contains :math:`\mathbb P(Z(S_{k,t})=j\mid Z(0)=i)`.
The expected value of :math:`S_{k,t}` is :math:`t` and the variance of :math:`S_{k,t}` is :math:`t/k`.
This means that

	.. math::

		P(t) = \lim_{k\to\infty}R(k/t)^k. 

The Erlangization method for approximating transition probabilities is discussed in [AsmussenEtAl2002]_, [MandjesTaylor2016]_ and [StanfordEtAl2011]_ for models related to birth-and-death processes. 
Similar to the uniformization method, error arises in the Erlangization method since the state space may need to be truncated, and the infinite limit above needs to be approximated by a finite :math:`k`. 

See :ref:`birdepy.probability(method='Erlang')` for the implementation of this method in BirDePy. 

Inverse Laplace transform
^^^^^^^^^^^^^^^^^^^^^^^^^
The Laplace transform of the transition function :math:`p_{i,j}(t)` is 

	.. math:: 

		f_{i,j}(s) = \mathcal L[p_{i,j}](s) = \int_0^\infty p_{i,j}(t)e^{-st} d t. 

Let :math:`\frac{u_1}{v_1+}\frac{u_2}{v_2+}\frac{u_3}{v_3+}\cdots` be a short-hand notation for the *continued fraction* 
	
	.. math::

		\dfrac{u_1}{v_1+\dfrac{u_2}{v_2+\dfrac{u_3}{v_3 +\cdots }}},

where :math:`(u_i,~i=1,2,\dots)` and :math:`(v_i,~i=1,2,\dots)` are sequences of real numbers. 
As first reported in [Murphy1975]_ and detailed in [CrawfordSuchard2012]_, the Laplace transform takes the continued fraction form

	.. math::

		f_{i,j}(s) = \left\{ \begin{array}{ll}
		\left(\prod_{k=j+1}^i\mu_k\right)\frac{B_j(s)}{B_{i+1}(s)+}\frac{B_i(s)\,a_{i+2}}{b_{i+2}(s)+}\frac{a_{i+3}}{b_{i+3}(s)+}\cdots, & \text{for } j \le i,\\[0.5em]
		\left(\prod_{k=i}^{j-1}\lambda_k\right) \frac{B_i(s)}{B_{j+1}(s)+}\frac{B_j(s)\,a_{j+1}}{b_{j+2}(s)+}\frac{a_{j+3}}{b_{j+3}(s)+}\cdots, & \text{for } j \ge i,
		\end{array}\right.

where 

	.. math::

		B_0(s)&=1,\\
		B_1(s)&=b_1(s),\quad\text{and}\\
		B_k(s)&=b_k(s)B_{k-1}(s)+a_kB_{k-2}(s),\quad\text{for } k\ge 2,

with :math:`a_1=1` and :math:`a_j = -\lambda_{j-2}{\mu_{j-1}}` for :math:`j\ge2`, and :math:`b_1(s)=s+\lambda_0` and :math:`b_j(s)=s+\lambda_{j-1}+\mu_{j-1}` for :math:`j\ge2`. 
An advantage of the continued fraction form is that it can be evaluated using the Lentz algorithm to a user-specified error tolerance. 
The Laplace transform can then be numerically inverted to obtain :math:`p_{i,j}(t)`. 

See :ref:`birdepy.probability(method='ilt')` for the implementation of this method in BirDePy. 


Diffusion approximation
^^^^^^^^^^^^^^^^^^^^^^^
Define :math:`(Z_i^{(r)},~r \in\mathbb N)` to be a parametric family of CTMCs which evolve according to transition rates :math:`\lambda_z^{(r)} = r\lambda_{z/r}` and :math:`\mu_z^{(r)} = r\mu_{z/r}` with :math:`Z_i^{(r)}(0)=i`. 
For :math:`s \in[0, t]` define the diffusion-scaled process 

	.. math::

		\tilde Z_i(s) = \lim_{r\to\infty}\sqrt{r}\left(\hat Z_i^{(r)}(s) - \hat z_i(s)\right),

where :math:`\hat Z_i^{(r)}:=\frac{1}{r} Z_i^{(r)}(s)` and :math:`\hat z_i(s)` satisfies 

	.. math::

		\frac{ d }{d s}\hat z_i(s) = \lambda_{\hat z_i(s)} - \mu_{\hat z_i(s)},\qquad \hat z_i(0) = i.

Loosely speaking, Theorem~3.5 in [Kurtz1971]_ can be used to show that, under some regularity conditions, :math:`(\hat Z_i^{(r)}(s),~s\in[0,t])` converges weakly, as :math:`r\to\infty`, in the space of cadlag functions on :math:`[0, t]` to a zero-mean Gaussian diffusion :math:`(\tilde Z_i(s),~s\in[0,t])` with variance 

	.. math::
		\sigma_i^2(s) := \text{Var}(\tilde Z(s)) = M(s)^2\int_0^s \left(\lambda_{\hat z(\tau)} + \mu_{\hat z(\tau)}\right) M(\tau)^{-2} d \tau

for each :math:`s\in[0, t]`, 
where :math:`M(s):= \exp\left(\int_0^s B(\tau)d \tau\right)` with :math:`B(\tau) = H\left(\hat z(\tau)\right)` defined in terms of 

	.. math::

		H(z) = \frac{d}{d z}\Big(\lambda_z - \mu_z\Big).

In particular this implies that 

	.. math::

		p_{i,j}(t) \approx \mathbb P(\tilde Z_i(t)=j). 

Hence :math:`p_{i,j}(t)` can be approximated by the probability density of a normally distributed random variable with mean :math:`\hat z_i(t)` and variance :math:`\sigma_i^2(t)` as given above.
 
This approach is very closely related to the functional central limit theorem. 
The diffusion approach to approximate transition probabilities is used in [RossEtAl2009]_, and discussed at length in [Allen2008]_. 

See :ref:`birdepy.probability(method='da')` for the implementation of this method in BirDePy. 


Ornstein--Uhlenbeck approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The diffusion approximation discussed above can be substantially simplified if it is assumed that :math:`Z` is fluctuating about a steady state point :math:`z_{\text{eq}}` of :math:`\hat z`. 
Such a point occurs when the birth rate is equal to the death rate, i.e., :math:`z_{\text{eq}}` satisfies :math:`\lambda_{z_{\text{eq}}} = \mu_{z_{\text{eq}}}`.
In this case, as argued in [RossEtAl2006]_, the limiting Gaussian diffusion is an Ornstein--Uhlenbeck process. 
Therefore, :math:`p_{i,j}(t)` can be approximated by the density of a normally distributed random variable with mean :math:`\tilde z= z_{\text{eq}} + e^{H(z_{\text{eq}})t}(i - z_{\text{eq}})` and variance :math:`\tilde\sigma^2 = \frac{\lambda_{z_{\text{eq}}} + \mu_{z_{\text{eq}}}}{2H(z_{\text{eq}})}\left(e^{2H(z_{\text{eq}})t}-1\right)`. 
This method assumes that :math:`z_{\text{eq}}` is asymptotically stable (i.e., :math:`H(z_{\text{eq}})<0`), and as such favors values of :math:`z_{\text{eq}}` that minimize :math:`H(z_{\text{eq}})`. 

See :ref:`birdepy.probability(method='oua')` for the implementation of this method in BirDePy. 


Galton--Watson approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A PSDBDP can be approximated by a piecewise-linear birth-and-death process by decomposing time into sub-intervals of finite length, and letting the per-individual birth and death rates be constant over each time interval; more precisely, if :math:`z` is the population size at the beginning of a time interval, then we let :math:`\lambda := \frac{1}{z}\lambda_{z}` and :math:`\mu := \frac{1}{z}\mu_{z}` over that interval. 

Suppose we want to approximate :math:`p_{i,j}(t)` for some :math:`i,j\in\mathcal{S}` and :math:`t> 0`. Consider a linear birth-and-death process :math:`\mathring Z^{(b)}` with per-individual birth rate :math:`\lambda = \frac{1}{b}\lambda_{b}` and per-individual death rate :math:`\mu = \frac{1}{b}\mu_{b}`, where possible choices of :math:`b` include :math:`b=i`, :math:`b=j`, :math:`b=\max(i,j)`, :math:`b=\min(i,j)`, and :math:`b=(i+j)/2`. 
The transition probability :math:`p_{i,j}(t)` can then be approximated by 

	.. math::

		p_{i,j}(t) \approx \mathring p_{i,j}(t) :=\mathbb P(\mathring Z^{(b)}(t)=j\mid \mathring Z^{(b)}(0)=i).

What constitutes a good choice of :math:`b` will depend highly on :math:`i`, :math:`j` and the model under study. 
This approximation is particularly convenient since it is well known (see [Guttorp1991]_) that :math:`\mathring p_{i,0}(t) = \beta_1(t)^{i}`, and for :math:`j\ge1`,

	.. math:: 
		\mathring p_{i,j}(t) = \sum_{k=\max(0, i-j)}^{i-1} \binom{i}{k} \binom{j -1}{i-k-1} \beta_{1}(t)^k\Big[\big\{1-\beta_{1}(t)\big\}\big\{1-\beta_{2}(t)\big\}\Big]^{j-k} \beta_{2}(t)^{j-i+k},

where :math:`\beta_1(t)` and :math:`\beta_2(t)` are given by 

	.. math:: 

		\beta_1(t) = \left\{\begin{array}{cc} \mu\{\exp\big((\lambda-\mu)t\big)-1\}/\{\lambda \exp\big((\lambda-\mu)t\big)-\mu\} & \text{if } \lambda\ne \mu, \\ 
		\lambda t/(1+\lambda t) & \text{if } \lambda = \mu. 
		\end{array}\right.

and

	.. math:: 

		\beta_2(t) = \left\{ \begin{array}{cc} \lambda\beta_1(t)/\mu& \text{if } \lambda\ne \mu, \\ 
		\beta_1(\tau) & \text{if } \lambda = \mu. 
		\end{array}\right.

See :ref:`birdepy.probability(method='gwa')` for the implementation of this method in BirDePy. 

When the binomial coefficients above cause numerical problems or take a long time to compute, an alternative expression developed by [DavisonEtAl2021]_ using a saddlepoint approximation [Butler2007]_ may be used. See :ref:`birdepy.probability(method='gwasa')` for the implementation of this method in BirDePy. 


Simulation
^^^^^^^^^^
Using :func:`birdepy.simulate.discrete()` or :func:`birdepy.gpu_functions.discrete()` it is possible to obtain :math:`k` realizations of :math:`Z(t)` conditional on :math:`Z(0)=i`. 
The proportion of these realisations which equal :math:`j` can be used to approximate the transition probability :math:`p_{i,j}(t)`. That is, 

.. math::

	p_{i,j}(t) \approx \frac{1}{k}\sum_{n=1}^k 1_{\{\hat Z_n(t)=j\}}
	
where :math:`\hat Z_l(t)` are simulated realisations of :math:`Z(t)` and :math:`I_{\{A\}}` evaluates to 1 when :math:`A` is true and otherwise evaluates to 0. 

See :ref:`birdepy.probability(method='sim')` and :ref:`birdepy.gpu_functions.probability()` for the implementations of this method in BirDePy. 


Summary
^^^^^^^
The table below summarizes methods described on this page, and gives the label used to call them in :ref:`birdepy.estimate()`. 

.. list-table:: Methods for computing transition probabilities.
   :widths: 18 20 20
   :header-rows: 1

   * - Method  
     - Label
     - Brief description
   * - Matrix exponential
     - :code:`'expm'`
     - Uses :math:`P(t)=\exp(Qt)`.
   * - Uniformization
     - :code:`'uniform'`
     - Evaluates probability using an approximating discrete-time process. 
   * - Erlangization
     - :code:`'Erlang'`
     - Evaluates probability at an Erlang-distributed time. 
   * - Inverse Laplace transform 
     - :code:`'ilt'`
     - Numerically inverts Laplace transform.
   * - Diffusion approx.
     - :code:`'da'`
     -  Approx. true model by a general diffusion-scaled model.
   * - Ornstein--Uhlenbeck approx.
     - :code:`'oua'`
     - Approx.\ true model by a simple diffusion process.
   * - Galton--Watson approximation
     - :code:`'gwa'`
     - Approx.~true model with a linear model.
   * - Saddlepoint approximation 
     - :code:`'gwasa'`
     - As above combined with a saddlepoint approx
   * - Simulation
     - :code:`'sim'`
     - Average of Monte Carlo samples.


These methods are also described in detail in [HautphennePatch2021a]_, which can be downloaded `here <https://github.com/birdepy/paper/blob/main/birdepy.pdf>`_. If you use BirDePy for published research, then please cite this paper.


.. [Grassman1977] Grassman, W., 1977. Transient solutions in Markovian queues. *European Journal of Operations Research*, 1(6):396--402. 

.. [vanDijkEtAl2018] van Dijk, N. M., van Brummelen, S. P. J., & Boucherie, R. J. (2018). Uniformization: Basics, extensions and applications. *Performance evaluation*, 118, 8-32.

.. [Jensen1953] Jensen, A. (1953). Markoff chains as an aid in the study of Markoff processes. *Scandinavian Actuarial Journal*, 1953(sup1), 87-91.

.. [AsmussenEtAl2002] Asmussen, S., Avram, F. and Usabel, M., 2002. Erlangian approximations for finite-horizon ruin probabilities, *ASTIN Bulletin: The Journal of the IAA*, 32(2)267--281. 

.. [MandjesTaylor2016] Mandjes, M. and Taylor, P., 2016. The running maximum of a level-dependent quasi-birth-death process, *Probability in the Engineering and Informational Sciences*, 30(2):212--223. 

.. [StanfordEtAl2011] Stanford, D.A., Yu, K. and Ren, J., 2011. Erlangian approximation to finite time ruin probabilities in perturbed risk models, *Scandinavian Actuarial Journal*, 2011(1):38--58. 

.. [Murphy1975] Murhy, J. A., & O'donohoe, M. R. (1975). Some properties of continued fractions with applications in Markov processes. *IMA Journal of Applied Mathematics*, 16(1), 57-71.

.. [CrawfordSuchard2012] Crawford, F.W. and Suchard, M.A. Transition probabilities for general birth-and-death processes with applications in ecology, genetics, and evolution. *Journal of Mathematical Biology*, 65(3):553--580. 

.. [Kurtz1971] Kurtz, T.J., 1971. Limit theorems for sequences of jump Markov processes, *Journal of Applied Probability*, 8(2):344--356. 

.. [RossEtAl2009] Ross, J. V., Pagendam, D. E., & Pollett, P. K. (2009). On parameter estimation in population models II: multi-dimensional processes and transient dynamics. *Theoretical Population Biology*, 75(2-3), 123-132.

.. [Allen2008] Allen, L. J. (2008). An introduction to stochastic epidemic models. In *Mathematical Epidemiology* (pp. 81-130). Springer, Berlin, Heidelberg.

.. [RossEtAl2006] Ross, J. V., Taimre, T., & Pollett, P. K. (2006). On parameter estimation in population models. *Theoretical Population Biology*, 70(4), 498-510.

.. [Guttorp1991] Guttorp, P., 1991.  Statistical inference for branching processes, Vol. 122. Wiley-Interscience.

.. [DavisonEtAl2021] Davison, A.C., Hautphenne, S. and Kraus, A., 2021. Parameter estimation for discretely observed linear birth‐and‐death processes. *Biometrics*, 77(1), pp.186-196.

.. [Butler2007] Butler, R.W., 2007.  Saddlepoint approximations with applications (Vol. 22). Cambridge University Press.
