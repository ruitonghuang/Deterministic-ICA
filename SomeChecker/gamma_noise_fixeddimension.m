%%This is to check the value of Gamma vs coherence of matrix A for a fixed
%%dimension d

A_pertube = 0.3:-0.01:0.01;
d = 6;
num_trail = 20000;

eps = ones(d,1);
Gamma = zeros(length(A_pertube), num_trail);
for ind = 1:length(A_pertube)
    ratio = A_pertube(ind);
    Atemp = randn(6);
    Atemp(:,4) =Atemp(:,3) + ratio * eps;
    
    for i = 1:num_trail
        phi = randn(d,1);
        psi = randn(d,1);
        Gamma(ind,i) = 1/calGamma(Atemp,phi,psi);
       
    end
    
end
Gamma_mean = mean(Gamma,2);
figure;
Xaxis = A_pertube;
plot(Xaxis, Gamma_mean);
%axis([0 0.05 0 200]);
%boxplot(log(Gamma'));