function RV = TimedParallelWrapper( objective_function, XV, DesiredNumberOfNonTimeouts, InitialTimeOut, varargin )
    persistent BestRunTime MaxRunTime

    N = size( XV, 2 );
    
    if isfinite( InitialTimeOut )
        if isempty( BestRunTime ) || isempty( MaxRunTime )
            Timeout = InitialTimeOut;
        else
            CurrentPool = gcp;
            TargetScale = DesiredNumberOfNonTimeouts ./ CurrentPool.NumWorkers;
            Timeout = max( BestRunTime * ( TargetScale + 2 ), MaxRunTime * ( TargetScale + 1 ) );
        end
    else
        Timeout = Inf;
    end
    
    fprintf( 'Using timeout: %g\n', Timeout );
    [ TPVOut, RunTimes ] = TimedParFor( @( i ) objective_function( XV( :, i ), varargin{:} ), 1:N, { Inf }, Timeout, false );
    RV = TPVOut{ 1 };

    oIndices = isfinite( RV ) & isfinite( RunTimes );
    
    sRV = RV( oIndices );
    RunTimes = RunTimes( oIndices );
    if ~isempty( sRV )
        [ ~, sIndices ] = sort( sRV );
        if length( sIndices ) >= DesiredNumberOfNonTimeouts
            sIndices = sIndices( 1:DesiredNumberOfNonTimeouts );
        else
            fprintf( 'Timeout appears to be too low. You may wish to modify the logic in TimedParallelWrapper.m.\nDesired %d, received %d.\n', DesiredNumberOfNonTimeouts, length( sIndices ) );
        end
        
        RunTimes = RunTimes( sIndices );
        NewMaxRunTime = max( RunTimes );
        NewBestRunTime = RunTimes( 1 );
        
        if isempty( MaxRunTime )
            MaxRunTime = NewMaxRunTime;
        else
            MaxRunTime = 0.95 * MaxRunTime + 0.05 * max( MaxRunTime * 0.8, NewMaxRunTime );
        end
        if isempty( BestRunTime )
            BestRunTime = NewBestRunTime;
        else
            BestRunTime = 0.95 * BestRunTime + 0.05 * max( BestRunTime * 0.8, NewBestRunTime );
        end
    end
end
