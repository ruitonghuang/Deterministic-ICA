function res = G(T, v)
d = size(v,1);
res = zeros(d,1);
for ind = 1:d
    res(ind,1) = v'*(T(:,:,ind)*v);
end