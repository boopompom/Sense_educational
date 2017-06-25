function  f = FFT(k,Id,Kd,varargin)
% Modified for 12D reconframe data
% Based on Miki Lustig his nufft operator

% Verbose options
if nargin < 4
    f.verbose=0;
else
    f.verbose=1;
end 

% Image space dimensions
f.Id=Id;

% K-space dimensions
f.Kd=Kd;

f.adjoint = 0;
f=class(f,'FFT');

% end
end

