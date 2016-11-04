function x = nsobol1( n )
    assert( n >= 2 );
    persistent p d
    if isempty( p )
        p = sobolset( 1, 'Skip', 2 );
    end
    if isempty( p ) || length( d ) < n
        d = norminv( net( p, n ) );
        d = d ./ d( 2 );
    end
    x = d( 1:n );
end
