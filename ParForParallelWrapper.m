function RV = ParForParallelWrapper( objective_function, XV, DesiredNumberOfNonTimeouts, Timeout, varargin )
    N = size( XV, 2 );
    RV = Inf( 1, N );
    
    StartTime = tic;
    for j = 1 : DesiredNumberOfNonTimeouts : N
        parfor i = j : min( N, ( j + DesiredNumberOfNonTimeouts - 1 ) )
            try
                RV( i ) = objective_function( XV( :, i ), varargin{:} ); %#ok<PFBNS>
            catch Error
                DisplayError( Error );
            end
        end
        ElapsedTime = toc( StartTime );
        if ElapsedTime > Timeout
            break;
        end
    end
end
