function res = bpsk(t,init,i)

sig = randi(2,1,t) - 1;
res = zeros(1,t);
list = [2,5,7,11,13,19]'.^(1/2);
kappa = [2.96,0.68,13.65,21.72,69.34,125.76].^(1/4);
for ind = 1:t
    if sig(1,ind) ==0
        res(1,ind) = i*sin(init+list(i)*ind);
    else
        res(1,ind) = -i*sin(init+list(i)*ind);
    end
end
 res = res/list(i);

