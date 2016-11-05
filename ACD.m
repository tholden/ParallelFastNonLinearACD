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

function [ xMean, BestFitness, Iterations, NEvaluations ] = ACD( FitnessFunction, xMean, sigma, LB, UB, A, b, MaxEvaluations, StopFitness, HowOftenUpdateRotation, Order, SearchDimension, Parallel )

    xMean = xMean(:);
    
    N = length( xMean );
    assert( N >= SearchDimension, 'The problem dimension should be weakly greater than the search dimension.' );    
    
    if isempty( sigma )
        sigma = ones( N, 1 );
    end
    if isempty( LB )
        LB = -Inf( N, 1 );
    end
    if isempty( UB )
        UB = Inf( N, 1 );
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
    %%% parameters
    k_succ = 1.1;       
    k_unsucc = 0.5;
    c1 = 0.5 / N;
    cmu = 0.5 / N;
    HowOftenUpdateRotation = max( 1, floor( HowOftenUpdateRotation ) ); %integer >=1

    sigma = min( sigma, ( UB - LB ) * 0.25 );

    BestFitness = 1e+30;
    NEvaluations = 0;

    B = eye(N,N);

    Iterations = 0;
    firstAE = true;
    ix = 0;
    somebetter = false;

    NoD = ceil( N / SearchDimension );
    
    UNPoints = 2 .^ ( 2 + Order ) - 1;
    NPoints = UNPoints .^ SearchDimension - 1;
    CellUalpha = repmat( { nsobol1( UNPoints )' }, 1, SearchDimension );
    [ CellUalpha{:} ] = ndgrid( CellUalpha{:} );
    alpha = cell2mat( cellfun( @(alphaCoord) reshape( alphaCoord( 2:end ), 1, NPoints ), CellUalpha, 'UniformOutput', false )' );
    
    allx = NaN( N, NPoints*NoD );
    allf = NaN( 1, NPoints*NoD );

    disp( [ 'Using ' num2str( NPoints ) ' points per iteration.' ] );
    if Parallel
        disp( 'This number should ideally be equal to or just below a multiple (possibly 1) of your number of parallel workers.' );
    end
    
    stream = RandStream( 'mt19937ar', 'Seed', 0 );
    ixPerm = randperm( stream, N );

    % -------------------- Generation Loop --------------------------------

    while (NEvaluations < MaxEvaluations) && (BestFitness > StopFitness)
        Iterations = Iterations + 1;
        ix = ix + 1;
        if ix > NoD
            ix = 1;
            ixPerm = randperm( stream, N );
        end

        qIndices = ( ix - 1 ) * SearchDimension + ( 1 : SearchDimension );
        qIndices( qIndices > N ) = qIndices( qIndices > N ) - N;
        
        qix = ixPerm( qIndices );
        %%% Sample NPoints candidate solutions
        dx = bsxfun( @times, sigma(qix,1)', B(:,qix) ); % shift along qix'th principal component, the computational complexity is linear
        
        x = zeros( N, NPoints );
        Fit = NaN( NPoints, 1 );

        if Parallel
            parfor iPoint = 1 : NPoints

                x( :, iPoint ) = clamp( xMean, dx * alpha( :, iPoint ), LB, UB, A, b );       % first point to test along qix'th principal component

                %%% Compute Fitness   
                try
                    Fit( iPoint ) = FitnessFunction( x( :, iPoint ) ); %#ok<PFBNS>
                catch
                end

            end
        else
            for iPoint = 1 : NPoints

                x( :, iPoint ) = clamp( xMean, dx * alpha( :, iPoint ), LB, UB, A, b );       % first point to test along qix'th principal component

                %%% Compute Fitness   
                try
                    Fit( iPoint ) = FitnessFunction( x( :, iPoint ) );
                catch
                end

            end
        end
        NEvaluations = NEvaluations + NPoints;

        %%% Who is the next mean point?  
        lsucc = false;
        [ minFit, minFitLoc ] = min( Fit );
        if minFit < BestFitness
            BestFitness = minFit;
            xMean = x( :, minFitLoc );
            lsucc = true;
        end

        %%% Adapt step-size sigma depending on the success/unsuccess of the previous search
        foundAlpha = max( abs( alpha( :, minFitLoc ) ) );
        
        if lsucc % increase the step-size
            sigma(qix,1) = sigma(qix,1) * foundAlpha * k_succ;     % default k_succ = 1.1
            somebetter = true;
        else            % decrease the step-size
            sigma(qix,1) = sigma(qix,1) * min( k_unsucc, foundAlpha * k_unsucc );   % default k_unsucc = 0.5
        end

        %%% Update archive 
        finiteFit = find( isfinite( Fit ) );
        finiteIndices = ( ix - 1 ) * NPoints + finiteFit;
        allx( :, finiteIndices ) = x( :, finiteFit );
        allf( 1, finiteIndices ) = Fit( finiteFit );
        
        if (ix == NoD) && somebetter && all( isfinite( allf ) ) %% we update our rotation matrix B every N=dimension iterations
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
