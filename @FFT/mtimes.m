function res = mtimes(fg,data) 
% Tom Bruijnen - University Medical Center Utrecht - 201704 

% Check what dimensions require new trajectory coordinates
Id=fg.Id;
Kd=fg.Kd;

if fg.adjoint==1    % Cartesian k-space to Cartesian image domain

    % Preallocate response matrix
    res=zeros(Id(1:4));
    for c=1:size(data,4)
        res(:,:,:,c)=fftshift(ifft2(ifftshift(data(:,:,:,c)),Id(1),Id(2)));
    end

else         % Cartesian image domain to Cartesian k-space 

    % Preallocate response matrix
    res=zeros(Kd(1:4));
    for c=1:size(data,4)
        res(:,:,:,c)=ifftshift(fft2(fftshift(data(:,:,:,c)),Kd(1),Kd(2)));
    end
        
end

% END  
end 
