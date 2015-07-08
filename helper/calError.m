function [error,p] = calError(A1,A2)

d = size(A1,2);
c = zeros(1,d);
indvec = 1:d;
P = perms(indvec);

for ind = 1:d
    A1(:,ind) = A1(:,ind)/norm(A1(:,ind),2);
    A2(:,ind) = A2(:,ind)/norm(A2(:,ind),2);
end

error = 100;
for trial = 1: size(P,1)
   
    A = A1(:,P(trial,:));
    for ind = 1:d
        c(1,ind) = (A2(:,ind)'*A(:,ind))/(A(:,ind)'*A(:,ind));
    end
    errorTemp = norm(A2 - A*diag(c),'fro');
    if errorTemp < error
        error = errorTemp;
        if nargout >1
            p = P(trial,:);
        end
    end
end