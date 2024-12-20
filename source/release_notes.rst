=============
Release Notes
=============

1.0.0
^^^^^^
Date:  November 30 2024

Contributors: Brendan Patch 

Notes:

Released a new version to coincide with publication of the associated paper in Journal
of Statistical Software.

0.0.26
^^^^^^
Date:  July 30 2023

Contributors: Brendan Patch 

Notes:

Since 0.0.12 several minor maintenance fixes have been applied. In addition:

1. The repository where the package is stored has been set up with Semantic Release to enable more efficient CI/CD. 
2. Python requirement is now 3.10 so that json files are more easily made accessible to the package. 
3. The opt_method argument for :func:`bd.estimate()` has been added to the help string of the function. 
4. Python notebook versions of code examples have been created and linked to from the home page. 

0.0.12
^^^^^^
Date: April 23 2022

Contributors: Brendan Patch 

Notes: 

1. Update requirements.txt

0.0.11
^^^^^^
Date: April 23 2022

Contributors: Brendan Patch 

Notes: 

1. Update utility_em.py with change in bound from [] to None when function a_opt is defined to accomodate recent changes in SciPy.


0.0.10
^^^^^^
Date: October 20 2021

Contributors: Brendan Patch 

Notes: 

1. Updated docstrings with arXiv reference for our paper. 

0.0.9
^^^^^
Date: October 6 2021

Contributors: Sophie Hautphenne and Brendan Patch 

Notes: 

1. Substantial revision to the ABC algorithm in :func:`bd.estimate()`. 
2. Minor changes to examples throughout. 

0.0.8
^^^^^
Date: August 29 2021

Contributors: Brendan Patch 

Notes: 

1. Provided docstrings for all functions, added more comments throughout codebase, and changed the formatting of examples inside docstrings so that they have a copy button on the website. 
2. Fixed a bug in :func:`bdg.discrete()` and :func:`bdg.probability()` where the default seed was set to 1 instead of chosen randomly. 
3. Model 'Poisson' is now treated the same as other models in all functions. 
4. Fixed a bug where models with exactly 1 unknown parameter had standard error reported as an array inside of a list instead of as a list. 
5. Fixed a bug where attempting to compute the carrying capacity for models where this does not exist would throw an exception. Instead a nan is returned for this attribute in the output. 

0.0.7
^^^^^
Date: August 9 2021

Contributors: Brendan Patch 

Notes: 

1. Added parameter `tau` as a kwarg for :func:`birdepy.probability()` when `method` is set to 'sim'. 

0.0.6
^^^^^
Date: August 7 2021

Contributors: Brendan Patch 

Notes: 

1. Added scikit-learn>=0.24.2 to requirements.txt.


0.0.5
^^^^^
Date: August 7 2021

Contributors: Sophie Hautphenne and Brendan Patch 

Notes: 

1. Added gwr-inversion to requirements.txt.
2. Changed where gwr-inversion is imported.
3. Ensured that :func:`np.ix_` receives arrays containing np.int32 objects.
4. Fixed a bug in the GPU version of the Verhulst model.
5. Attribute se of the output of :func:`bd.estimate()` is now a list instead of an array.
6. Simulation based standard errors and confidence intervals now use the estimated parameter as the initial condition for the optimisation procedure.
7. Improved how having a terminal state equal to 0 is handled by :func:`bd.probability(method='gwasa')`.
8. Made it so an a value equal to 0 raises an exception in probability.gwasa.w_fun so that the high precision version of the function is invoked. 
9. Updated how constraints are handled when differential evolution is used as an optimizer. Constraints are now always specified as dictionaries. 
10. When polishing is used to enhance the output of differential evolution (which is the default action) it is done using method 'SLSQP' of :func:`scipy.minimize()` instead of method 'trust-constr' (which causes problems in other parts of the package).
11. Added the parameter 'stat' to the ABC algorithm so that 'mean' or 'median' can be used as a summary statistics for the posterior distribution. 
12. Transition probability approximation methods 'gwa' and 'gwasa' now use a midpoint to anchor the linear approximation underlying their output.
13. Fixed a bug in :func:`bd.forecast()` related to displaying xticks.
14. Changed the default labels of x and y axis in bd.forecast().
15. Added the 'export' parameter to :func:`bd.estimate()` and :func:`bd.forecast()` which allows plots to be exported as tex files. 
16. Removed a constant from the distance function in :func:`bd.estimate(framework='abc')`. 
17. Renamed :func:`bdg.discrete_gpu` to :func:`bdg.discrete`.
18. Renamed :func:`bdg.probability_gpu` to :func:`bdg.probability`.
19. Renamed attribute 'err' of :func:`bd.estimate()`'s output to 'val'. 
20. Renamed argument 'num_bs_samples' of :func:`bd.estimate()` to 'num_samples'. 


0.0.4
^^^^^
Date: July 21 2021

Contributors: Brendan Patch

Notes: 

1. Replaced iltcme.json with iltcme.py to fix a bug where iltcme.json would not load properly


0.0.3
^^^^^
Date: July 20 2021

Contributors: Sophie Hautphenne and Brendan Patch 

Notes: 

1. Many bug fixes. 
2. Added basic readme. 
3. Substantial revision of bd.forecast() function. 


0.0.2
^^^^^
Date: July 14 2021

Contributors: Sophie Hautphenne and Brendan Patch 

Notes: 

1. Many bug fixes. 
2. This version was not uploaded to PyPI. 


0.0.1
^^^^^

Date: July 6 2021

Contributors: Sophie Hautphenne and Brendan Patch 

Notes: 

1. Initial release for testing
