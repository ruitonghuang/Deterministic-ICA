function error = calError_entropy(A,y)

%d = size(A,1);
%deter = abs(det(A));
%scal = nthroot(1/deter,d);
%A = scal*A;
d = size(y,1);
yhat = A\y;
co = HShannon_spacing_VKDE_initialization(1);
error = 0;
for i = 1:d
    error = error + HShannon_spacing_VKDE_estimation(yhat(i,:),co);
end

error  = error + log(abs(det(A)));