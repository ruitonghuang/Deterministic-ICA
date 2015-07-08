function A = G2(phi,y)

num = size(y,2);
c = (phi'*y)*(phi'*y)';
A = c*(y*y')/num^2;
