function res = dynamic_indexing(A,b,c,varargin)
% A is input matrix where to get all entries from besides dimension b.
% From dimension b only get element c
% So function performs res = A(:,:,2) for b=3 and c=2.
%
% If varargin is 1 then the backward problem is solved.
% So A(:,:,2) = res
% Tom Bruijnen

sz=size(A);

if b == 0  || b > numel(sz) || c > sz(b)
    if ~isempty(varargin)
        res=varargin{1};
    else
        res=A;
    end
    return
end

inds=repmat({1},1,ndims(A));
for n=1:numel(sz);inds{n}=1:sz(n);end
inds{b}=c;

if nargin <4
    res=A(inds{:});
else
    Adims=[];
    for n=1:numel(sz)
        if n~=b
            Adims=[Adims sz(n)];
        else
            Adims=[Adims 1];
        end
    end
    A(inds{:})=reshape(varargin{1},Adims);
    res=A;
end

% END
end