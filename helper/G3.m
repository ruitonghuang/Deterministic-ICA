function A = G3(phi,y)

num = size(y,2);
c = y*(y'*phi);
A = (1/num^2)*c*c';