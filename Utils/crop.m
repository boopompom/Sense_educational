function y = crop(X,n)
% X is the 2D matrix, X should be even valued in size 
% n is the number of rows/columns to crop from the center

dims=size(X);
cp=dims/2;

y=X(cp(1)-floor(n/2)+1:cp(1)+ceil(n/2),cp(2)-floor(n/2)+1:cp(2)+ceil(n/2),:,:,:);
end