% ---------------------------------------------------------------
% Adaptive Coordinate Descent. To be used under the terms of the BSD license 
% Author : Ilya Loshchilov, Marc Schoenauer, Michele Sebag, 2012.  
% Further work: Tom Holden, 2016.
% e-mail: ilya.loshchilov@gmail.com marc.schoenauer@inria.fr michele.sebag@lri.fr 
% URL:http://www.lri.fr/~ilya
% REFERENCE: Loshchilov, I., Schoenauer, M. , Sebag, M. (2011). Adaptive Coordinate Descent. 
%    N. Krasnogor et al. (eds.)
%    Genetic and Evolutionary Computation Conference (GECCO) 2012,
%    Proceedings, ACM.  http://hal.inria.fr/docs/00/58/75/34/PDF/AdaptiveCoordinateDescent.pdf
% This source code includes the Adaptive Encoding procedure by Nikolaus Hansen, 2008
% ---------------------------------------------------------------

function [ xMean, BestFitness, Iterations, NEvaluations ] = ACD0( FitnessFunction, xMean, Sigma, MinSigma, LB, UB, A, b, MaxEvaluations, StopFitness, HowOftenUpdateRotation )

    xMean = xMean(:);
    
    N = length( xMean );
    
    if isempty( Sigma )
        Sigma = ones( N, 1 );
    end
    if isempty( MinSigma )
        MinSigma = ones( N, 1 ) * sqrt( eps );
    end
    if isempty( LB )
        LB = -Inf( N, 1 );
    end
    if isempty( UB )
        UB = Inf( N, 1 );
    end
    if length( MinSigma ) == 1
        MinSigma = repmat( MinSigma, N, 1 );
    end
    if length( LB ) == 1
        LB = repmat( LB, N, 1 );
    end
    if length( UB ) == 1
        UB = repmat( UB, N, 1 );
    end
    if isempty( A )
        A = zeros( 0, N );
    end
    if isempty( b )
        b = zeros( 0, 1 );
    end
    if isempty( StopFitness )
        StopFitness = -Inf;
    end
    %%% parameters
    k_succ = 1.1;       
    k_unsucc = 0.5;
    c1 = 0.5 / N;
    cmu = 0.5 / N;
    HowOftenUpdateRotation = floor(HowOftenUpdateRotation); %integer >=1

    Sigma = min( Sigma, ( UB - LB ) * 0.25 );

    BestFitness = 1e+30;
    NEvaluations = 0;

    B = eye(N,N);

    Iterations = 0;
    firstAE = true;
    ix = 0;
    somebetter = false;


    allx = NaN(N,2*N);
    allf = NaN(1,2*N);

    % -------------------- Generation Loop --------------------------------

    while (NEvaluations < MaxEvaluations) && (BestFitness > StopFitness)
        Iterations = Iterations + 1;
        ix = ix + 1;
        if ix > N
            ix = 1;
        end
        
        %%% Sample two candidate solutions
        dx = Sigma(ix,1) * B(:,ix); % shift along ix'th principal component, the computational complexity is linear
        x1 = clamp( xMean, -dx, LB, UB, A, b );       % first point to test along ix'th principal component
        x2 = clamp( xMean, dx, LB, UB, A, b );        % second point to test is symmetric to the first one on the same ix'principal component
        
        %%% Compute Fitness
        try
            Fit1 = FitnessFunction( x1 );
        catch
            Fit1 = NaN;
        end
        NEvaluations = NEvaluations + 1;
        try
            Fit2 = FitnessFunction( x2 );   
        catch
            Fit2 = NaN;
        end
        NEvaluations = NEvaluations + 1;

        %%% Who is the next mean point?  
        lsucc = false;
        if Fit1 < BestFitness
            BestFitness = Fit1;
            xMean = x1;
            lsucc = true;
        end
        if Fit2 < BestFitness
            BestFitness = Fit2;
            xMean = x2;
            lsucc = true;
        end

        %%% Adapt step-size sigma depending on the success/unsuccess of the previous search
        if lsucc % increase the step-size
            Sigma(ix,1) = Sigma(ix,1) * k_succ;     % default k_succ = 1.1
            somebetter = true;
        else            % decrease the step-size
            Sigma(ix,1) = Sigma(ix,1) * k_unsucc;   % default k_unsucc = 0.5
        end

        if all( Sigma < MinSigma )
            break
        end

        %%% Update archive 
        posx1 = ix*2 - 1;
        posx2 = ix*2;
        if isfinite( Fit1 )
            allx(:,posx1) = x1(:,1);
            allf(1,posx1) = Fit1;
        end
        if isfinite( Fit2 )
            allx(:,posx2) = x2(:,1);
            allf(1,posx2) = Fit2;
        end

        if (ix == N) && somebetter && all( isfinite( allf ) ) %% we update our rotation matrix B every N=dimension iterations
            somebetter = false;

            [~, arindex] = sort(allf,2,'ascend');
            allxbest = allx(:,arindex(1:N));
            if firstAE
                ae = ACD_AEupdateFAST([], allxbest, c1, cmu,HowOftenUpdateRotation);    % initialize encoding
                ae.B = B;
                ae.Bo = ae.B;     % assuming the initial B is orthogonal
                ae.invB = ae.B';  % assuming the initial B is orthogonal
                firstAE = false;
            else 
                ae = ACD_AEupdateFAST(ae, allxbest, c1, cmu,HowOftenUpdateRotation);    % adapt encoding 
            end
            B = ae.B;
        end    
        
        if rem(Iterations,1000) == 0
            disp([ num2str(Iterations) ' ' num2str(NEvaluations) ' ' num2str(BestFitness) ]);
        end
    end
end
