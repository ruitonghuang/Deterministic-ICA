function res = calNomial(seq,x)

num_ind = size(seq,2);

if num_ind == 2
    if seq(1) == seq(2)
        res = mean(x(seq(1),:).^2);
    else 
        res = 0;
    end
elseif num_ind ==3
    if seq(1) == seq(2) && seq(1) == seq(3)
        res = mean(x(seq(1),:).^3);
    else
        res = 0; 
    end
elseif num_ind == 4
    if seq(1) == seq(2) && seq(1) == seq(3) && seq(1) == seq(4)
        res = mean(x(seq(1),:).^4);
    elseif seq(1) == seq(2) && seq(4) == seq(3)
        res = mean(x(seq(1),:).^2) * mean(x(seq(3),:).^2);
    elseif seq(1) == seq(3) && seq(4) == seq(2)
        res = mean(x(seq(1),:).^2) * mean(x(seq(2),:).^2);
    elseif seq(1) == seq(4) && seq(3) == seq(2)
        res = mean(x(seq(1),:).^2) * mean(x(seq(2),:).^2);
    else res = 0;
    end
else
    res = 0;
end
%d = max(seq);
%ind = zeros(1,d);
%for i = 1:size(seq,2)
%    ind(1,seq(1,i)) = ind(1,seq(1,i))+1;
%end

%res = 1;
%for i = 1:d
%    res = res/(ind(1,i)+1);
%end