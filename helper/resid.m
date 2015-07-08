function [theta,res] = resid(A,Per)
Ainv = inv(A);
d = size(A,1);

V  = zeros(d^2,d);

for ind = 1:d
    V(:,ind) = reshape(A(:,ind)*Ainv(ind,:), [d^2,1]);
end

PerV = reshape(Per,[d^2,1]);

theta = V\PerV;
residule = PerV - V*theta;
res = reshape(residule,[d,d]);
