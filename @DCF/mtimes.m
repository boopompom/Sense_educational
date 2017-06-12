function output = mtimes(dcf,input)

dims=size(input);
% Examine dimensions of dcf
for j=1:numel(size(input))
    if size(dcf.w,j)==size(input,j)
        dims(j)=1;
    end
end

% Multiply k-space data with DCF
W=repmat(dcf.w,dims);
output=input.*W;


% END  
end  
