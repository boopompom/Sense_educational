 function [alpha, beta, ok] = nufft_alpha(N, J, K, alpha, beta)
%	Determine alpha and beta, as associated with the scaling factors
%	for min-max interpolation, from arguments.
%	in
%		alpha	[L,1]	Fourier series coefficients of scaling factors
%				or a string (see below)
%		beta		scale gamma=2pi/K by this in Fourier series
%				typically is K/N (me) or 0.5 (Liu)
%	out
%		alpha,beta
%
%	Copyright 2001-3-30	Jeff Fessler	The University of Michigan

%	if no arguments
if nargin < 3
	help(mfilename)
	error(mfilename)
return
end

if ~isvar('alpha') || isempty(alpha)
	alpha = [1];	% default Fourier series coefficients of scaling factors
end
if ~isvar('beta') || isempty(beta)
	beta = 0.5;	% default is Liu version for now
end


%
%	see if 'best' alpha is desired
%
if ischar(alpha)
	if streq(alpha, 'uniform')
		alpha = [1];
		beta = 0.5;
	else
		if streq(alpha, 'best')
			L = 0;
		elseif streq(alpha, 'best,L=1')
			L = 1;
		elseif streq(alpha, 'best,L=2')
			L = 2;
		else
			error 'unknown alpha argument'
		end
		[alpha, beta, ok] = nufft_best_alpha(J, L, K/N);
		if ~ok
			tmp = 'optimal alpha unknown for J=%d, K/N=%g, L=%d';
			warning(sprintf(tmp, J, K/N, L))
		end
	end
end
