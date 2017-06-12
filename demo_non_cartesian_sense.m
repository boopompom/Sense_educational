%% DEMO for Cartesian 2D SENSE for scientific MRI meeting 12 June
clear all;close all;clc;echo off
addpath(genpath(pwd))

% This is a demo to play around with certain factors involved
% in the generalized iterative SENSE reconstruction. The simulations
% are performed on the shepp-logan phantom or brain images. If you would
% rather have a brain image to look at set the parameter below to 1.
brain=0;

%% Set input parameters
NC=5;           % Number of receive coils 
M=256;          % Matrix size of image (has to be 2^n in this implementation)
csm_sigma=4;    % Width of gaussian CSM profiles, this affects the g-factor!!
R=2;            % Undersampling factor (has to be 2^n in this implementation)
verbose=1;      % 1=visuals, 0=none

% Some noise characteristics, note that changing these affects the results
% a lot
noise_sampling=.0005;         % Noise on samples % of kmax
noise_csm=.001;               % Noise on acquired csm % of Imax

% Check input parameters 
M=makedyadic(M);   % Make M 2^n
R=makedyadic(R);   % Make R 2^n

%% Create phantom and coil sensitivity maps
% Create phantom without any noise
x_true=phantom(M);

% If you wanted the brain phantom, do some interpolation 
if brain;load brain;[x1,y1]=meshgrid(1:512,1:512);[x2,y2]=meshgrid(1:M,1:M);
    x_true=interp2(x1,y1,abs(x_true),512/M*x2,512/M*y2,'linear');end

% Generate Gaussian shaped coil sensitivity maps
csm=generate_csm(NC,M,csm_sigma);

% Visualization
if verbose
    figure(1);imshow3(abs(cat(3,x_true,csm)),[],[1 (NC+1)]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Phantom and coil sensitivity maps (CSM)');
    pause();
end


%% Why does the classical formalism NOT work for non-Cartesian sampling?
% Lets examine the corresponding point-spread-function (PSF)
% This simply means that we perform an experiment where we put all ones on
% our samples and mimic undersampling by entering zero values at specific phase
% encoding lines

% Lets examine the PSF for a couple of undersampled Radial cases, dont
% bother about the syntax for now. Just look at the results
for Rn=1:4                   % Undersampling factor
    angles=0:2*pi/(M*pi/2/Rn-1):2*pi;
    ktraj=radial_trajectory(angles,M);
    weights=radial_weights(ktraj);
    W=DCF(weights); % density operator
    F=FG2D(ktraj,[M M 1],[M size(ktraj,2) 1]); % Fourier operator

    K_psf=ones(size(ktraj));
    I_psf(:,:,Rn)=F'*(W*(K_psf));  % 2D Fast inverse Fourier Transform
    I_psf(:,:,Rn)=I_psf(:,:,Rn)/max(max(abs(I_psf(:,:,Rn))));
end

if verbose;
    figure(2),imshow3(log(abs(I_psf)+1),[0 0.05],[1 4]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Point Spread Function for R=1 till R=4');
    pause();
end

% You can observe that the energy gets folded over almost all pixels to
% some kind of amount. This is what we often refer to as incoherent aliasing.
% If we follow the same approach for the Cartesian sampling case we will have 
% to invert the S^(H) * S matrix which will have size of M^4 and we have to do this for every pixel!
% This is computionally too expensive and thats why it doesnt work.

%% Define radial trajectory, weights and NUFFT operators
% First we have to define a couple of operators associated with radial
% sampling. This includes the non-uniform fourier transform, density
% compensation operator and the correct k-space trajectory and density
% weights. It is not important to understand the syntax for now.

% Define the radial angles, uniformly spaced 
angles=0:2*pi/(M*pi/2/R-1):2*pi; % [radians]

% Get the corresponding k-space trajectory
ktraj=radial_trajectory(angles,M);

% Get the corresponding density
weights=radial_weights(ktraj);

% Define the coil operator, fourier transform and density operator
W=DCF(weights); % density operator
S=CC(csm);      % Sense operator
F=FG2D(ktraj,[M M NC],[M size(ktraj,2) NC]); % Fourier operator

%% Lets mimic the MR experiment
% We take our phantom image (x_true), multiply it with the coil sensitivity maps (S operator)
% and subsequently perform a Fourier transform (F operator). This mimics the process at the
% scanner!
y=F*(S'*x_true);

% Lets mimic real-live condition which include noise on the measurements
% First define a function to add gaussian complex noise
add_noise=@(x,lvl)(x+lvl*randn(size(x))+1j*lvl*randn(size(x)));

% Add noise to our measurements "y"
y=add_noise(y/max(abs(y(:))),noise_sampling);

% Add noise to our measured coil sensitivity maps
csm=add_noise(csm,noise_csm);

% Perform an inverse fourier transform on the undersampled and noisy data
x_aliased=F'*(W*y);

% Now lets examine the individual aliased coil images.
% You can observe that every receive channel exhibitis different degrees of aliasing.
% This is the principle we use to unfold the images.

% Visualization
if verbose
    figure(3),imshow3(abs(cat(3,x_aliased/max(abs(x_aliased(:))),csm)),[],[2 NC]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Individual aliased coil images and noise corrupted csm');
    pause();
end

%% Lets perform the iterative sense recon
% So instead of the direct inversion we discussed in the previous demo, we
% apply the generalized SENSE formulation that aims to minimize the
% following cost functional: ||FSx-y||^2 
% Where x = our current coil combined image estimate   [M x M]
%       y = the measured k-space data [M x 1/R*M*pi/2]
%       F = Fourier transform operator
%       S = Coil sensitivity map operator
%
% What we do to solve this is we have to define a function handle in matlab
% that performs the forward operation and adjoint operation.
% Forward: take our current coil combined image estimate --> apply the coil sensitivities --> and do a
% fourier transform to mimic the k-space sampling.
% Adjoint: take the multi-coil k-space data, do an inverse Fourier
% transform and apply the coil sensitivities to spread the data across the
% receivers.
%
% Now we create a struct (params) that contains all the information to do the steps 
% described above. And then we send this to the Matlab lsqr function. This
% is all we need to solve a l2-norm cost function! We don't need to have
% any knowledge about the problem, as long as we can define the forward and
% adjoint operation and know that the problem is linear and convex. The
% code to perform the forward and adjoint operation is covered in
% iterative_sense.m

% Generate the structure
params={};
params.S=S;                    % Coil operator
params.F=F;                    % Fourier operator
params.W=DCF(sqrt(weights));   % Density operator
params.Id=size(csm);           % Image dimensions
params.Kd=size(y);             % K-space dimensions
params.y=params.W*y;           % Measured k-space data

% Send the struct to the matlab solver (lsqr)
func=@(x,tr) iterative_sense(x,params,tr);
[recon,~,~,~,~,~]=lsqr(func,params.y(:),1E-10,10);
x_recon=reshape(recon,params.Id(1:2));

if verbose
   tmp=S*x_aliased;
   figure(4);
   imshow3(cat(3,abs(tmp)/max(abs(tmp(:))),abs(x_recon)/max(abs(x_recon(:)))),[],[1 2]); 
   set(gcf, 'Position', get(0, 'Screensize'));
   title('Undersampled image and iterative sense image');
end


% END
