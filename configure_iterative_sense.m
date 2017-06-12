function [res,lsvec] = configure_iterative_sense(params)

% LSQR settings
ops=optimset('Display','off');

% Function handle for (regularized) lsqr
func=@(x,tr) iterative_sense(x,params,tr);
s=matrix_to_vec(params.y);

% LSQR
[tmp,~,~,~,~,lsvec]=lsqr(func,s,1E-10,params.N_iter);
res=reshape(tmp,params.Id(1:3));

% Append resvec if it converges reached prior to N_iter
if size(lsvec,1)<params.N_iter;lsvec(end+1:params.N_iter)=0;end

% END
end