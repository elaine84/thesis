These are the source files for my PhD thesis, [Accelerating Markov chain Monte Carlo via parallel predictive prefetching][1]. 

Abstract
--------

We present a general framework for accelerating a large class of widely used Markov chain Monte Carlo (MCMC) algorithms. This dissertation demonstrates that MCMC inference can be accelerated in a model of parallel computation that uses speculation to predict and complete computational work ahead of when it is known to be useful. By exploiting fast, iterative approximations to the target density, we can speculatively evaluate many potential future steps of the chain in parallel. In Bayesian inference problems, this approach can accelerate sampling from the target distribution, without compromising exactness, by exploiting subsets of data. It takes advantage of whatever parallel resources are available, but produces results exactly equivalent to standard serial execution. In the initial burn-in phase of chain evaluation, it achieves speedup over serial evaluation that is close to linear in the number of available cores.

Dissertation committee
----------------------

Margo Seltzer, Ryan P. Adams, Eddie Kohler

References
----------

Elaine Angelino. Accelerating Markov chain Monte Carlo via parallel predictive prefetching.
PhD thesis, School of Engineering and Applied Sciences, Harvard University, 2014. 
[Harvard version][1] (git tag `submit`). [Living version][2].

[Fetching][3] is the associated implementation.	

Elaine Angelino, Eddie Kohler, Amos Waterland, Margo Seltzer, Ryan P. Adams.
[Accelerating MCMC via parallel predictive prefetching][4].
In *30th Conference on Uncertainty in Artificial Intelligence*, UAI ’14, 2014.


[1]: http://www.eecs.harvard.edu/~elaine/thesis-harvard.pdf
[2]: http://www.eecs.harvard.edu/~elaine/thesis-living.pdf
[3]: https://github.com/elaine84/fetching
[4]: http://auai.org/uai2014/proceedings/individuals/286.pdf

Cloning
-------

	$ git clone git@github.com:elaine84/thesis.git

