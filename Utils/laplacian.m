function res = laplacian(img)
% Perform laplacian on image N-dimensional

% Get dimensions and reshape
dims=size(img);
img=reshape(img,[dims(1) dims(2) prod(dims(3:end))]);

% Generate kernel
laplacianKernel = [-1,-1,-1;-1,8,-1;-1,-1,-1]/8;

% Apply kernel
for n=1:prod(dims(3:end))
    res(:,:,n) = imfilter(double(img(:,:,n)), laplacianKernel);
end

% Reshape back
res=reshape(res,dims);

% END
end