# ParallelFastNonLinearACD
A parallel version of Loshchilov, Schoenauer and Sebag's (2012) adaptive coordinate descent algorithm for non-linear optimization

## Key reference
Loshchilov, I., Schoenauer, M. , Sebag, M. (2011). **Adaptive Coordinate Descent.**

in N. Krasnogor et al. (eds.), Genetic and Evolutionary Computation Conference (GECCO) 2012, Proceedings, ACM.

http://hal.inria.fr/docs/00/58/75/34/PDF/AdaptiveCoordinateDescent.pdf

## Files
 * **ACD0** is a slightly cleaned up version of the original code from here: https://uk.mathworks.com/matlabcentral/fileexchange/37057-fast-adaptive-coordinate-descent-for-non-linear-optimization
    * Extra features include support for linear constraints (with clipping to the bound) and more CMA-ES like input arguments.
 * **ACD1** is a parallelized version of ACD0, with parallelization controlled by the `Order` input.
    * `Order` 0 uses just 2 points per iteration, and will often be slower than serial code.
    * `Order` 1 uses 6 points per iteration, and may be faster than serial code when the objective function is slow to evaluate, and you have more than 6 cores.
    * `Order` `n` uses `2^(2+n)-2` points per iteration.
    * In all cases, the parallelization is providing improved univariate search in each iteration, reducing the required number of iterations.
 * **ACD2** is like ACD1, but searches over pairs of coordinates, rather than coordinate by coordinate.
    * At `Order` 0, this uses 8 points per iteration, so may be perfect for 8 core desktops.
    * At `Order` `n`, this uses `(2^(2+n)-1)^2-1` points per iteration.
 * **ACD3** is like ACD2, but searches over triples of coordinates, rather than coordinate pair by coordinate pair.
    * At `Order` 0, this uses 26 points per iteration, so may be perfect for e.g. 32 core servers.
    * At `Order` `n`, this uses `(2^(2+n)-1)^3-1` points per iteration.    
