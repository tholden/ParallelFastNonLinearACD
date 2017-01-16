function RV = SerialWrapper( objective_function, XV, DesiredNumberOfNonTimeouts, Timeout, varargin )
    N = size( XV, 2 );
    RV = Inf( 1, N );
    
    StartTime = tic;
    for i = 1 : N
        try
            RV( i ) = objective_function( XV( :, i ), varargin{:} );
        catch Error
            DisplayError( Error );
        end
        if i >= DesiredNumberOfNonTimeouts
            ElapsedTime = toc( StartTime );
            if ElapsedTime > Timeout
                break;
            end
        end
    end
end
