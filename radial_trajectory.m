function traj = radial_trajectory(angles,M)

% Calculate sampling point on horizontal spoke
x=linspace(0,M-1,M)'-(M-1)/2;

% Modulate the phase of all the successive spokes
k=zeros([M,numel(angles)]);
for l=1:numel(angles)
	k(:,l)=x*exp(1j*angles(:,l));
end

% Normalize
traj=k/M;

% END
end