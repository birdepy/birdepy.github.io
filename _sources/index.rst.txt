============================================
BirDePy: Birth-and-death Processes in Python
============================================
`BirDePy` is a Python package for working with continuous time `birth-and-death processes <https://en.wikipedia.org/wiki/Birth-and-death_process>`_. It includes functions which can be used to approximate transition probabilities (:ref:`birdepy.probability()<birdepy.probability()>`), estimate parameter values from data (:ref:`birdepy.estimate()<birdepy.estimate()>`), simulate sample paths (:ref:`birdepy.simulate.discrete()<birdepy.simulate.discrete()>` and :ref:`birdepy.simulate.continuous()<birdepy.simulate.continuous()>`), and generate forecasts (:ref:`birdepy.forecast()<birdepy.forecast()>`). The main focus of the package is the estimation of parameter values from discretely-observed sample paths, however the much more straightforward case of estimation of parameter values from continuously observed sample paths is also included. 

You can install `BirDePy` using `pip <https://pypi.org/project/pip/>`_. To do this open Python and execute::

    pip install birdepy

You may also need to install some dependencies (listed below).  

`BirDePy` is developed on Github, the source code is found `here <https://github.com/birdepy/birdepy_project/tree/main/src/birdepy>`_. 

The package is also described in detail in [HautphennePatch2021a]_, which can be downloaded `here <https://github.com/birdepy/paper/blob/main/birdepy.pdf>`_. If you use BirDePy for published research, then please cite this paper.

If you use our new piecewise approximation simulation algorithm, then please also cite our paper [HautphennePatch2021b]_, which can be downloaded `here <https://github.com/birdepy/paper/blob/main/psdbdp_simulation.pdf>`_. 

Usage Example 
-------------
Suppose we are interested in the susceptible-infective-susceptible model. This is the Verhulst model with rate of spread given by :math:`\gamma`, recovery rate given by :math:`\nu`, population size given by :math:`1/\alpha` and :math:`\beta=0`: ::

   import birdepy as bd
   import numpy as np
   model = 'Verhulst'
   rate_of_spread = 0.75
   recovery_rate = 0.25
   population_size = 1000
   true_parameters = [rate_of_spread, recovery_rate, 1/population_size, 0]
   simulation_horizon = 100
   initial_number_infected = 10
   obs_times = np.arange(0, simulation_horizon+1, 1)

BirDePy's function :func:`birdepy.simulate.discrete()` can be used to simulate a possible sample path: ::

   number_infected = bd.simulate.discrete(true_parameters, model, initial_number_infected,
                                          obs_times, seed=2021)

The path can be plotted: ::

   import matplotlib.pyplot as plt

   plt.step(obs_times, number_infected, 'r', where='post', color='tab:purple')
   plt.ylabel('infected population')
   plt.xlabel('day')
   plt.show()

It looks like this:

.. image:: example_fig_1.png

If we assume :math:`\beta=0` is known, BirDePy's function :func:`birdepy.estimate()` can be used to estimate the other parameter values: ::

   initial_guess = [0.5]*3
   p_bounds = [[1e-6, 1]]*3
   est = bd.estimate(obs_times, number_infected, initial_guess, p_bounds,
                     framework='dnm', likelihood='da', model=model,
                     known_p=[0], idx_known_p=[3])
   print('Estimate:', est.p)
   print('Standard error:', est.se)

BirDePy produces an estimate of [0.7590, 0.2411, 0.0010] with standard errors [0.0566 0.0386 5.23e-05].
The argument of `framework` is set to 'dnm' so :ref:`direct numerical maximisation <Direct Numerical Maximization>` is used to find a maximum likeliood estimate. 
The argument of `likelihood` is set to 'da' so a :ref:`diffusion approximation <Diffusion approximation>` is used to approximate the likelihood function. 

If we suppose that :math:`\alpha=0.001` is known, a confidence region can be produced for the remaining two unknown parameters by setting `ci_plot` to 'True': ::

   est = bd.estimate(obs_times, number_infected, [0.5]*2, [[1e-6, 1]]*2,
                     framework='dnm', likelihood='da', model=model,
                     known_p=[1/population_size, 0], idx_known_p=[2, 3], ci_plot=True,
                     xlabel='$\gamma$', ylabel='$\\nu$')

It looks like this:

.. image:: example_fig_2.png

A forecast confidence interval can be made: ::

   future_t = np.arange(100, 151, 1)
   bd.forecast('Verhulst', number_infected[-1], future_t, est.p, cov=est.cov,
               p_bounds=[[1e-6, 1]]*2, idx_known_p=[2, 3], known_p=[0.001, 0], 
               interval='confidence',  xticks=np.arange(100, 151, 10))

It looks like this:

.. image:: example_fig_3.png

The confidence intervals contain future mean behavior with high probability. 

More Examples
-------------
For further examples see our examples GitHub repo `here <https://github.com/birdepy/birdepy_examples>`_.

It contains Python scripts and Python notebooks. 


Dependency List
---------------
Required:

* ``numpy``
* ``mpmath``
* ``scipy``
* ``matplotlib``
* ``gwr_inversion``

Optional:

* ``numba``
* ``cudatoolkit``

.. toctree::
   :caption: Theoretical Background
    
   processes
   simulation
   probabilities
   estimation
   forecasting

.. toctree::
   :caption: Core Functions API

   api_probability
   api_estimate
   api_simulate_discrete
   api_simulate_continuous
   api_forecast

.. toctree::
   :caption: CUDA Functions API

   api_cuda_simulate
   api_cuda_probabilities

.. toctree::
   :caption: Development

   contributing
   release_notes


.. [HautphennePatch2021a] Hautphenne, S., and Patch, B., (2021). Birth-and-death Processes in Python: The BirDePy Package. arXiv preprint arXiv:2110.05067.

.. [HautphennePatch2021b] Hautphenne, S. and Patch, B., 2021. Simulating population-size-dependent birth-and-death processes using CUDA and piecewise approximations. In *Proceedings of the International Congress on Modelling and Simulation 2021*. (to appear)
