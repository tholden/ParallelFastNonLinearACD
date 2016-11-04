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

function [ xmean, bestFit, iter, neval ] = ACD3( fitnessfct, xmean, sigma, LB, UB, A, b, MAX_EVAL, stopfitness, howOftenUpdateRotation, Order )

    N = length( xmean );
    if isempty( sigma )
        sigma = ones( N, 1 );
    end
    if isempty( LB )
        LB = -Inf( N, 1 );
    end
    if isempty( UB )
        UB = Inf( N, 1 );
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
    howOftenUpdateRotation = floor(howOftenUpdateRotation); %integer >=1

    sigma = min( sigma, ( UB - LB ) * 0.25 );

    bestFit = 1e+30;
    neval = 0;

    B = eye(N,N);

    iter = 0;
    firstAE = true;
    ix = 0;
    somebetter = false;

    Dimension = 3;
    NoD = ceil( N / Dimension );
    
    UNPoints = 2 .^ ( 2 + Order ) - 2;
    NPoints = ( 1 + UNPoints ) .^ Dimension - 1;
    Ualpha = [ 0; nsobol1( UNPoints ) ];
    [ alphaX, alphaY, alphaZ ] = meshgrid( Ualpha, Ualpha, Ualpha );
    alpha = [ alphaX( 2:end ); alphaY( 2:end ); alphaZ( 2:end ) ];
    
    allx = NaN( N, NPoints*NoD );
    allf = NaN( 1, NPoints*NoD );

    disp( [ 'Using ' num2str( NPoints ) ' points per iteration.' ] );
    disp( 'This number should be higher than your number of parallel workers.' );
    
    stream = RandStream( 'mt19937ar', 'Seed', 0 );
    ixPerm = randperm( stream, N );

    % -------------------- Generation Loop --------------------------------

    while (neval < MAX_EVAL) && (bestFit > stopfitness)
        iter = iter + 1;
        ix = ix + 1;
        if ix > NoD
            ix = 1;
            ixPerm = randperm( stream, N );
        end

        qIndices = ( ix - 1 ) * Dimension + ( 1 : Dimension );
        qIndices( qIndices > N ) = qIndices( qIndices > N ) - N;
        
        qix = ixPerm( qIndices );
        %%% Sample NPoints candidate solutions
        dx = bsxfun( @times, sigma(qix,1)', B(:,qix) ); % shift along qix'th principal component, the computational complexity is linear
        
        x = zeros( N, NPoints );
        Fit = zeros( NPoints, 1 );

        parfor iPoint = 1 : NPoints

            x( :, iPoint ) = clamp( xmean, dx * alpha( :, iPoint ), LB, UB, A, b );       % first point to test along qix'th principal component
        
            %%% Compute Fitness   
            Fit( iPoint ) = fitnessfct( x( :, iPoint ) ); %#ok<PFBNS>
            
        end
        neval = neval + NPoints;

        %%% Who is the next mean point?  
        lsucc = false;
        [ minFit, minFitLoc ] = min( Fit );
        if minFit < bestFit
            bestFit = minFit;
            xmean = x( :, minFitLoc );
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
                ae = ACD_AEupdateFAST([], allxbest, c1, cmu,howOftenUpdateRotation);    % initialize encoding
                ae.B = B;
                ae.Bo = ae.B;     % assuming the initial B is orthogonal
                ae.invB = ae.B';  % assuming the initial B is orthogonal
                firstAE = false;
            else 
                ae = ACD_AEupdateFAST(ae, allxbest, c1, cmu,howOftenUpdateRotation);    % adapt encoding 
            end
            B = ae.B;
        end    
        
        if rem(iter,1000) == 0
            disp([ num2str(iter) ' ' num2str(neval) ' ' num2str(bestFit) ]);
        end
    end
end

function x = clamp( xmean, dx, LB, UB, A, b )
    dxpos = dx > 0;
    dxneg = dx < 0;
    Adx = A * dx;
    Adxpos = Adx > 0;
    alpha = min( [ 1; ( b( Adxpos ) - A( ( Adxpos ), : ) * xmean ) ./ Adx( Adxpos ); ( LB( dxneg ) - xmean( dxneg ) ) ./ dx( dxneg ); ( UB( dxpos ) - xmean( dxpos ) ) ./ dx( dxpos ) ] );
    x = xmean + alpha * dx;
end
