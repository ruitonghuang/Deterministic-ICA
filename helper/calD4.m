function res = calD4(x)

d = size(x,1);
res = 0;
xmean = mean(x,2);

for ind1 = 1:d
    for ind2 = 1:d
        for ind3 = 1:d
            for ind4 = 1:d
                ind = [ind1,ind2,ind3,ind4];
                temp = abs(mean(x(ind1,:).*x(ind2,:).*x(ind3,:).*x(ind4,:)) - calNomial(ind,x));
                if res < temp
                    res = temp;
                end
            end
        end
    end
end

for ind1 = 1:d
    for ind2 = 1:d
        for ind3 = 1:d
            ind = [ind1,ind2,ind3];
            temp = abs(mean(x(ind1,:).*x(ind2,:).*x(ind3,:)) - calNomial(ind,x));
            if res < temp
                res = temp;
            end
        end
    end
end

for ind1 = 1:d
    for ind2 = 1:d
        ind = [ind1,ind2];
        temp = abs(mean(x(ind1,:).*x(ind2,:)) - calNomial(ind,x));
        if res < temp
            res = temp;
        end
    end
end

for ind1 = 1:d
   % ind = [ind1];
    temp = abs(mean(x(ind1,:)));
    if res < temp
        res = temp;
    end
end