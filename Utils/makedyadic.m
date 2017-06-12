function res = makedyadic(int)
% Round number to next dyadic

p=1;

while (p < int)
    p=p*2;
end

res=p;
% END
end