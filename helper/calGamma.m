function res = calGamma(A,phi,psi)

vec1 = psi'*A;
vec1 = vec1.^2;
vec2 = phi'*A;
vec2 = vec2.^2;
d = size(A,1);
res = inf;
for ind1 = 1:d
    for ind2 = ind1+1:d
        temp = abs(vec1(ind1)/vec2(ind1) - vec1(ind2)/vec2(ind2));
        if res > temp
            res = temp;
        end
    end
end