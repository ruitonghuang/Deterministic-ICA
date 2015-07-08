function error = calError_approx(A1,A2)

d = size(A1,1);
for ind = 1:d
    A1(:,ind) = A1(:,ind)/norm(A1(:,ind),2);
    A2(:,ind) = A2(:,ind)/norm(A2(:,ind),2);
end
D = A2\A1;


v1= zeros(d,1);
v2= zeros(d,1);
v3= zeros(d,1);
v4= zeros(d,1);

for ind =  1:d
    v1(ind) = norm(D(ind,:),1)-1;
    v2(ind) = norm(D(:,ind),1)-1;
    v3(ind) = norm(D(ind,:),2)^2-1;
    v4(ind) = norm(D(:,ind),2)^2-1;
end

error = norm(v1,2)^2+ norm(v2,2)^2+ norm(v3,1)+ norm(v4,1); 