function res = zpad2(x,varargin)
%  res = zpad(x,sx,sy)
%  Zero pads a 2D matrix around its center.
%
%
%  res = zpad(x,sx,sy,sz,st)
%  Zero pads a 4D matrix around its center
%
%  
%  res = zpad(x,[sx,sy,sz,st])
%  same as the previous example
%
%
% (c) Michael Lustig 2007

s=cell2mat(varargin);

    m = size(x);
    if length(m) < length(s)
	    m = [m, ones(1,length(s)-length(m))];
    end
	
    if sum(m==s)==length(m)
	res = x;
	return;
    end

    res = zeros(s);
    
    for n=1:length(s)
	    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s)=x;',cmd);
    eval(cmd);





