function res = iterative_sense(x,params,transp_flag)

% Forward and backward operations
if strcmp(transp_flag,'transp') % Nonuniform k-space uncombined --> uniform image space combined  
    
    % Vector to matrix
    x=vec_to_matrix(x(1:prod(params.Kd)),params.Kd);

    % Fourier transform
    res=params.F'*(params.W*x);
    
    % Coil sensitivity maps operator
    res=params.S*res;
    
    % Matrix to vector
    res=matrix_to_vec(res);
    
elseif strcmp(transp_flag,'notransp') % uniform image combined --> nonuniform k-space uncombined
    

    % Vector to matrix
    x=vec_to_matrix(x,params.Id(1:2));
    
    % S operator
    x=params.S'*x;
    
    % Fourier transform
    res=params.W*(params.F*x);
    
    % Vectorize
    res=matrix_to_vec(res);    
   
end

% END
end