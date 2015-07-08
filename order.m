function res = order(seq)

[seq, I1] =  sort(seq);
len = length(seq);

gap = seq(2:len) - seq(1:len-1);
truegap = zeros(1,len);
truegap(1) = gap(1);
truegap(len) = gap(len-1);
for ind = 2:len-1
    truegap(ind) = min(gap(ind-1),gap(ind));
end

[~,I2] = sort(truegap);

res = I1(I2);