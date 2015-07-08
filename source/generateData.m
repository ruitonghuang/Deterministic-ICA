function [y,A] = generateData(d,num,lambda) %generate num dimension d samples

A = rand(d);
for ind  = 1:d
   A(:,ind) = A(:,ind)/norm(A(:,ind),2);
end
A = A + eye(d);
x = rand(d,num);
eps = randn(d,num);
y = A*x + lambda*eps;
