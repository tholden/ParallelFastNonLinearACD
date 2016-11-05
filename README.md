# ParallelFastNonLinearACD
A parallel version of Loshchilov, Schoenauer and Sebag's (2012) adaptive coordinate descent algorithm for non-linear optimization

## Key reference
Loshchilov, I., Schoenauer, M. , Sebag, M. (2011). **Adaptive Coordinate Descent.**

in N. Krasnogor et al. (eds.), Genetic and Evolutionary Computation Conference (GECCO) 2012, Proceedings, ACM.

http://hal.inria.fr/docs/00/58/75/34/PDF/AdaptiveCoordinateDescent.pdf

## Files
 * **ACD0** is a slightly cleaned up version of the original code from here: https://uk.mathworks.com/matlabcentral/fileexchange/37057-fast-adaptive-coordinate-descent-for-non-linear-optimization
   * Extra features include support for linear constraints (with clipping to the bound) and slightly different input arguments (an initial point, and sigma, the initial search standard deviation(s)).
 * **ACD** is a (potentially) parallelized version of ACD0, with parallelization controlled by the `Order`, `SearchDimension` and `Parallel` inputs.
   * With `SearchDimension` 1, the search is coordinate by coordinate.
     * At `Order` 0 this uses just 2 points per iteration, and will often be slower than serial code if `Parallel` is enabled.
     * At `Order` 1 this uses 6 points per iteration, and may be faster than serial code when the objective function is slow to evaluate, and you have more than 6 cores.
     * At `Order` `n` this uses `2^(2+n)-2` points per iteration.
     * In all cases, the parallelization is providing improved univariate search in each iteration, reducing the required number of iterations.
   * With `SearchDimension` 2, the search is over pairs of coordinates, rather than coordinate by coordinate.
     * At `Order` 0, this uses 8 points per iteration, so may be perfect for 8 core desktops.
     * At `Order` `n`, this uses `(2^(2+n)-1)^2-1` points per iteration.
   * With `SearchDimension` 3, the search is over triples of coordinates, rather than coordinate pair by coordinate pair.
     * At `Order` 0, this uses 26 points per iteration, so may be perfect for e.g. 32 core servers.
     * At `Order` `n`, this uses `(2^(2+n)-1)^3-1` points per iteration. 
   * With `SearchDimension` `d`, the search is over `d`-tuples of coordinates.
     * At `Order` 0, this uses `3^d-1` points per iteration.
     * At `Order` `n`, this uses `(2^(2+n)-1)^d-1` points per iteration.

## Function usage
`[ xMean, BestFitness, Iterations, NEvaluations ] = ACD( FitnessFunction, xMean, sigma, LB, UB, A, b, MaxEvaluations, StopFitness, HowOftenUpdateRotation, Order, SearchDimension, Parallel );`

Inputs:
 * `FitnessFunction`: The objective function, a function handle.
 * `xMean`: The initial point.
 * `Sigma`: The initial search radius. Either a scalar, or a vector of search radiuses by coordinate, with the same number of elements as xMean.
 * `MinSigma`: The minimum search radius. The search will stop when all coordinates of sigma are below this value. Either a scalar, or a vector of minimum search radiuses by coordinate, with the same number of elements as xMean.
 * `LB`: A lower bound on the search for `xMean`. Either empty, a scalar, or a vector of lower bounds by coordinate, with the same number of elements as xMean.
 * `UB`: A lower bound on the search for `xMean`. Either empty, a scalar, or a vector of upper bounds by coordinate, with the same number of elements as xMean.
 * `A`: The `A` matrix from the inequality `A*x <= b`. May be empty if `b` is also empty.
 * `b`: The `b` vector from the inequality `A*x <= b`. May be empty if `b` is also empty.
 * `MaxEvaluations`: The maximum number of total function evaluations. (Set to `Inf` if this is empty.)
 * `StopFitness`: The terminal fitness. (Set to `-Inf` if this is empty.)
 * `HowOftenUpdateRotation`: How often the rotation should be updated. On problems with slow objective functions, this should be equal to `1`. Larger values may speed up computation if the objective function is very fast.
 * `Order`: Determines the number of points to use to search along each direction. A (small) non-negative integer.
 * `SearchDimension`: Determines how many dimensions to search in simultaneously. A (small) positive integer.
 * `Parallel`: Determines whether to use a `parfor` loop to invoke the objective function. A logical.
 
 
 Ouputs:
  * `xMean`: The optimal point.
  * `BestFitness`: The value of the objective at that point.
  * `Iterations`: The number of iterations performed.
  * `NEvaluations`: The number of function evaluations performed. 
