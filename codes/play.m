function play( s, tempo, Fs )
% spline, 2008.5
if nargin < 3
	Fs = 8192;
end
fprintf( 'Loading...\n' );
freq = 27.5 * 2.0 .^ ( ( 0 : 87 ) / 12 );
freq( 89 ) = 0;
%ot = [ 1 0.8 0.7 0.4 0.5 0.2 0.45 0.2 0.25 0.25 0.4 0.15 0.2 0.15 0.1 0.05 ];
ot = [ 1 0.2 0.4 0.08 0.2 0.05 0.1 0.04 ];
ot( 1 ) = 1;
p = 2 * pi * ( 1 : length( ot ) );
af = 0.4;
dt = 60.0 / tempo * Fs;
w = zeros( 1, round( sum( s( :, 4 ) ) * dt ) + Fs );
bt = 0;
for i = 1 : size( s, 1 )
	t = 1 : round( s( i, 2 ) * dt );
	wt = t * 0;
	for j = 1 : length( ot )
		wt = wt + ot( j ) * sin( p( j ) * freq( s( i, 1 ) ) / Fs * t );
	end
	w( bt + t ) = w( bt + t ) + s( i, 3 ) * af .^ ( t / Fs ) .* wt;
	bt = bt + round( s( i, 4 ) * dt );
end
sound( 0.15 * w, Fs );