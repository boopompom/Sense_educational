function  fg = NUFFT(k,Id,Kd,varargin)
% Modified for 12D reconframe data
% Based on Miki Lustig his nufft operator

% Verbose options
if nargin < 4
    fg.verbose=0;
else
    fg.verbose=1;
end

% Number of data chunks
num_data=numel(k);   

% Image space dimensions
fg.Id=Id;

% K-space dimensions
fg.Kd=Kd;

% Mix the readouts and samples in advance
for n=1:num_data;fg.k=reshape(k,[Kd(1)*Kd(2) 1]);end;clear k % [coils=1]

% Input for nufft_init
Jd = [5,5];     % Kernel width of convolution
for n=1:num_data
    Nd=Id(1:2);
    Gd = [Nd*2];    % Overgridding ratio
    n_shift = Nd/2;
end

% Create a seperate struct for all the dimensions that need seperate trajectories

om=[real(fg.k)';imag(fg.k)']'*2*pi;
fg.st = nufft_init(om, Nd, Jd, Gd, n_shift,'minmax:kb');

fg.phase = 1;
fg.w=1;
fg.adjoint = 0;
fg.mode = 2;   % 2= complex image
fg=class(fg,'NUFFT');

% end
end

