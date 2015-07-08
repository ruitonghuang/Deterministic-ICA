function res = HRadical(y)

d = size(y,1);
num = size(y,2);
m = floor(sqrt(num));
res = 0;
for ind = 1:d
    x = sort(y(ind,:));
    temp = 0;
    for i = 1:num-m
        temp = temp + log((num+1)/m*(x(i+m)-x(i)));
    end
    res = res + 1/(num-m)*temp;
end