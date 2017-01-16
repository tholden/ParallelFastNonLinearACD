function RV = TimedParallelWrapper( objective_function, XV, DesiredNumberOfNonTimeouts, InitialTimeOut, varargin )
    persistent Timeout

    N = size( XV, 2 );
    
    if isfinite( InitialTimeOut )
        if isempty( Timeout )
            CTimeout = InitialTimeOut;
            fprintf( 'Initial timeout: %g\n', CTimeout );
        else
            CTimeout = Timeout;
        end
    else
        CTimeout = Inf;
    end
    
    [ TPVOut, RunTimes ] = TimedParFor( @( i ) objective_function( XV( :, i ), varargin{:} ), 1:N, { Inf }, CTimeout, false );
    RV = TPVOut{ 1 };

    oIndices = isfinite( RV ) & isfinite( RunTimes );
    
    sRV = RV( oIndices );
    RunTimes = RunTimes( oIndices );
    if ~isempty( sRV )
        [ ~, sIndices ] = sort( sRV );
        if length( sIndices ) >= DesiredNumberOfNonTimeouts
            sIndices = sIndices( 1:DesiredNumberOfNonTimeouts );
        else
            fprintf( 'Timeout appears to be too low. You may wish to modify the logic in ParallelWrapper.m.\nDesired %d, received %d.\n', DesiredNumberOfNonTimeouts, length( sIndices ) );
        end
        RunTimes = RunTimes( sIndices );
        MaxRunTime = max( RunTimes );
        BestRunTime = RunTimes( 1 );
        
        CurrentPool = gcp;
        TargetScale = DesiredNumberOfNonTimeouts ./ CurrentPool.NumWorkers;
        TimeoutTarget = max( BestRunTime * ( TargetScale + 2 ), MaxRunTime * ( TargetScale + 1 ) );
        if isempty( Timeout )
            Timeout = TimeoutTarget;
        else
            Timeout = 0.95 * Timeout + 0.05 * max( Timeout * 0.8, TimeoutTarget );
        end
    end
    fprintf( 'New timeout: %g\n', Timeout );
end
