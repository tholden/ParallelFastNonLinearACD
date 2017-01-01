
% ---------------------------------------------------------------
% Adaptive Coordinate Descent. To be used under the terms of the BSD license
% Author : Ilya Loshchilov, Marc Schoenauer, Michele Sebag, 2012.  
% e-mail: ilya.loshchilov@gmail.com marc.schoenauer@inria.fr michele.sebag@lri.fr 
% URL:http://www.lri.fr/~ilya
% REFERENCE: Loshchilov, I., Schoenauer, M. , Sebag, M. (2011). Adaptive Coordinate Descent. 
%    N. Krasnogor et al. (eds.)
%    Genetic and Evolutionary Computation Conference (GECCO) 2012,
%    Proceedings, ACM.  http://hal.inria.fr/docs/00/58/75/34/PDF/AdaptiveCoordinateDescent.pdf
% This source code includes the Adaptive Encoding procedure by Nikolaus Hansen, 2008
% ---------------------------------------------------------------


N = 10; % try 2,40,100,1000

addpath( './tests' );

ffunc = @frosenbrock;
%ffunc = @fsphere;
%ffunc = @felli;
%ffunc = @ftablet;
%ffunc = @fdiffpow;
%ffunc = @fschwefel;
%ffunc = @fcigar;
MAX_EVAL = 1000000*N;
x_a = -5.0; x_b = 5.0;
ftarget = 1e-10;
fcurrent = 1e+30;

Order = 2;
NonProductSearchDimension = 2;
ProductSearchDimension = 2;
Parallel = true;

howOftenUpdateRotation = 1; % at each iteration -> quadratic time complexity of the algorithm, but need less function evaluations to reach the optimum
% howOftenUpdateRotation = floor(dim/10); % every N/10 iterations -> linear time complexity of the algorithm, but need more function evaluations to reach the optimum
% howOftenUpdateRotation = dim; % every N iterations

nevaltotal = 0;
itertotal = 0;
tic
while (nevaltotal < MAX_EVAL) && (fcurrent > ftarget)
    maxeval_available = MAX_EVAL - nevaltotal;
    [xmean, fcurrent, iter, neval] = ACD(ffunc,zeros(N,1),2.5,0,[],[],[],[],maxeval_available,ftarget,howOftenUpdateRotation,Order,NonProductSearchDimension,ProductSearchDimension,Parallel);
    nevaltotal = nevaltotal + neval;
    itertotal = itertotal + iter;
end;

toc;
disp([num2str(itertotal) ' ' num2str(nevaltotal) ' ' num2str(fcurrent)]);

