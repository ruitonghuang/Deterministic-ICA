function A = G1(phi,y)

num = size(y,2);
d = size(y,1);
v = phi'*y;
temp = zeros(d);

for ind = 1:num
    q = (1/sqrt(num))*v(1,ind)* y(:,ind);
    temp = temp +q*q';
end

A = temp;

%for ind = 1: num
%   c = (phi'*y(:,ind))^2;
%   A = A + c * (y(:,ind) * y(:,ind)'); 
%end
%A = (1/num)*A;