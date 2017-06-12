function output = mtimes(cc,data) 

cc.S=double(cc.S);

Id=size(data);Id(end+1:13)=1;
if cc.adjoint==1 % S
    output=sum(data.*conj(cc.S),3);

else % S^-1
        output(:,:,:)=repmat(data,[1 1 size(cc.S,3)]).*cc.S;
end


% END
end  
