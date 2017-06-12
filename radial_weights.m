function weights = radial_weights(ktraj)

M=size(ktraj,1);
N=size(ktraj,2);
% Create a Ram-Lak filter
w=ones(M,1);
for i=1:M
    w(i)=abs(M/2 - (i - .5));
end
w=pi/N*w;
w=w/max(abs(w(:)));

% Partition into dynamics
weights=repmat(w,[1 N]);

% END
end