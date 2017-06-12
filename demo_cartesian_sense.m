%% DEMO for Cartesian 2D SENSE for scientific MRI meeting 12 June
clear all;close all;clc;echo off
addpath(genpath(pwd))

% This is a demo to play around with certain factors involved
% in the most simplistic Cartesian SENSE reconstruction. The simulations
% are performed on the shepp-logan phantom or brain images. If you would
% rather have a brain image to look at, than set the parameter below to 1.
brain=0;

%% Set input parameters
NC=5;           % Number of receive coils 
M=256;          % Matrix size of image (has to be 2^n in this implementation)
csm_sigma=2;    % Width of gaussian CSM profiles, this affects the g-factor!!
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

% If you wanted the brain phantom, do some interpolation for the dimensions
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

%% Why does the classical formalism work for cartesian sampling?
% Lets examine the corresponding point-spread-function (PSF)
% This simply means that we perform an experiment where we put all ones on
% our samples and mimic undersampling by entering zero values at specific phase
% encoding lines

% Lets examine the PSF for a couple of undersampled Cartesian cases
for Rn=1:4                   % Undersampling factor
    K_psf=zeros(M,M);        % K-space samples
    K_psf(1:Rn:end,:)=1;     % Uniform undersampling
    I_psf(:,:,Rn)=fftshift(ifft2(ifftshift(K_psf)))*Rn;  % 2D Fast inverse Fourier Transform
end

if verbose;
    figure(2),imshow3(abs(I_psf),[],[1 4]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Point Spread Function for R=1 till R=4');
    pause();
end

% You can observe that the energy gets folded over R pixels in the
% points spread function. That allows us to unfold this energy over these
% discrete numbers of pixels using the coil sensitivity maps.

%% Lets mimic the MR experiment
% We take our phantom image (x_true), multiply it with the coil sensitivity maps (csm)
% and subsequently perform a Fourier transform. This mimics the process at the
% scanner!
y=ifftshift(fft2(fftshift(repmat(x_true,[1 1 NC]).*csm)));

% Lets mimic real-live condition which include noise on the measurements
% First define a function to add gaussian complex noise
add_noise=@(x,lvl)(x+lvl*randn(size(x))+1j*lvl*randn(size(x)));

% Add noise to our measurements "y"
y=add_noise(y/max(abs(y(:))),noise_sampling);

% Add noise to our measured coil sensitivity maps
csm=add_noise(csm,noise_csm);

% Perform the undersampling with factor R by zerofilling lines (uniform
% fashion but keep the center line of k-space)
y_us=zeros(size(y));
idx=sort([floor(M/2)-R+1:-R:1 floor(M/2)+1:R:M]);
y_us(idx,:,:)=y(idx,:,:);

% Perform an inverse fourier transform on the undersampled and noisy data
x_folded=fftshift(ifft2(ifftshift(y_us)));

% Now lets examine the individual aliased coil images.
% You can observe that every receive channel exhibits different degrees of aliasing.
% This is the principle we use to unfold the images.

% Visualization
if verbose
    figure(3),imshow3(abs(cat(3,x_folded/max(abs(x_folded(:))),csm)),[],[2 NC]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Individual aliased coil images and noise corrupted csm');
    pause();
end

%% Lets perform the unfolding according to the formula below
% Now that we have our aliased multi-coil images and the corresponding
% spatial sensitivities we can try to unfold them according to the method
% described by Pruessman et al.
% This is the governing equation from the article (also shown in the
% video)
%
% Equation (1): p=[S^(H)*Psi^(-1)*S]^(-1) * [S^(H) *Psi^(-1)*m]
% p=unfolded pixel value    [R x 1]
% Psi=noise covariance matrix --> Identity matrix for now [NC x NC]
% S=senstivity matrix       [NC x R]
% H=conjugate transpose operation  
% m=folded pixel value in each receiver coil [NCx1]

% Assign psi
psi=eye(NC);

% Preallocate unfolded matrix p
p=zeros(M);

% Assign aliasing distance from PSF
dx=M/R;

% Loop over all pixels and perform the inversion described above
% remember that we performed the undersampling in the X-dimension
for px=1:M/R
    % Calculate coordinates where aliasing is on (known from the PSF)
    coords=mod(px+((1:R)-1).*dx-1,M)+1;
    for py=1:M;
        % Generate vector m from equation (1)
        m=permute(x_folded(px,py,:),[3 1 2]);
        % Generate matrix S 
        S=permute(csm(coords,py,:),[3 1 2]);
        % Do inversion to unfold
        p(coords,py)=(conj(transpose(S))*psi*S)\(conj(transpose(S))*m);
        % Same as this line but faster --> p(coords,y)=inv(conj(transpose(S))*psi*S)*conj(transpose(S))*m;
    end
end

% Visualization
if verbose;
    tmp=fftshift(ifft2(ifftshift(y)));
    true=abs(sum(conj(csm).*tmp,3))./sum(abs(csm).^2,3);
    figure(4),imshow3(cat(3,true/max(abs(true(:))),...
    abs(p/max(abs(p(:))))),[],[1 2]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Recon of fully sampled (noisy) data and unfolded undersampled image');
end

%% Why does the unfolding fail dramatically for higher accelerations (R~NC)?
% The matrix [S^(H)*psi*S] has to be inverted. When this is
% ill-conditioned, i.e. when the rows or colums of the matrix 
% are not independent anymore, i.e. the matrix is singular. In a physical
% context this would be when there is barely any difference in spatial
% sensitivities between the aliased pixels across multiple coils.

% One solution for this is to constrain the solution of the inversion of
% the matrix. This is often referred to as a regularized matrix inversion,
% this is what I believe the scanner reconstruction does as well. The most
% simplest form of regularization is zeroth order Tikhonov. I won't go into
% details about this, but an implementation is shown below. Note that I
% suspect the scanner to use a much more advanced regularized
% reconstruction, I didn't have time to try to approximate such a
% reconstruction. This implementation doesnt behave very good yet.

return
ops=optimset('Display','off');
lambda=0.00001;
for px=1:M/R
    % Calculate coordinates where aliasing is on (know from PSF)
    coords=mod(px+((1:R)-1).*dx-1,M)+1;
    for py=1:M;
        % Generate vector m from equation (1)
        m=permute(x_folded(px,py,:),[3 1 2]);
        % Generate matrix S 
        S=permute(csm(coords,py,:),[3 1 2]);
        % Create tikhonov regularization matrix
        T=lambda*eye(R);
        % Do regularized inversion to unfold
        [p_reg(coords,py),~]=lsqr([conj(transpose(S))*psi*S;T],[conj(transpose(S))*m;zeros(R,1)]);
    end
end

% Visualization
if verbose
    figure(5),imshow(abs(p_reg/max(abs(p_reg(:)))),[]);
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Tikhonov regularized unfolded image');
end

% END
