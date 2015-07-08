function res = calKappa(x)
d = size(x,1);
%num = size(x,2);

x4 = x.^4;
x2 = x.^2;
res = zeros(d,1);

for i = 1:d
    res(i,1) = mean(x4(i,:)) - 3*mean(x2(i,:))^2;
end