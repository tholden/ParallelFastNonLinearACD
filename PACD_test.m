
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
ftarget = 1e-10;

tic;
[xmean, fcurrent, iter, neval] = PACD( @(XV,mu) ffunc( XV ), zeros(N,1),[],[],[],[],false);
% [xmean, fcurrent, iter, neval] = PACD( @(XV,mu) ParForParallelWrapper( ffunc, XV, mu, 0 ), zeros(N,1),[],[],[],[],false);
% [xmean, fcurrent, iter, neval] = PACD( @(XV,mu) TimedParallelWrapper( ffunc, XV, mu, 5 ), zeros(N,1),[],[],[],[],false);
% [xmean, fcurrent, iter, neval] = PACD( @(XV,mu) SerialWrapper( ffunc, XV, mu, 0 ), zeros(N,1),[],[],[],[],false);
toc;

disp([num2str(iter) ' ' num2str(neval) ' ' num2str(fcurrent)]);

