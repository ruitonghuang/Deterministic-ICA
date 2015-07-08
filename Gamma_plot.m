num = 200;
num_A = 100;
d_max = 95;
Xaxis = 5:10:d_max;

gammaA1 = zeros(d_max,num_A, num);
gammaA2 = zeros(d_max,num_A, num);
gammaA3 = zeros(d_max,num_A, num);
gammaA4 = zeros(d_max,num_A, num);
gamma_R = zeros(d_max,num_A, num);

for d = Xaxis 
    for ind_A = 1: num_A
        A1 = randn(d);
        
        TT = randn(d); %randn(d,1)*ones(1,d)+cohe*
            [TT,~,~] = svd(TT);
            A = randn(d,1);
            A = 2*A/norm(A,2);            
            A = A*ones(1,d);
        A2 = A+0.3*TT;
        
        TT = randn(d); %randn(d,1)*ones(1,d)+cohe*
            [TT,~,~] = svd(TT);
            A = randn(d,1);
            A = 2*A/norm(A,2);            
            A = A*ones(1,d);
        A3 = A+0.05*TT;
        
        TT = randn(d); %randn(d,1)*ones(1,d)+cohe*
            [TT,~,~] = svd(TT);
            A = randn(d,1);
            A = 2*A/norm(A,2);            
            A = A*ones(1,d);
        A4 = A+0.005*TT; 
        
        B = randn(d);
        [R,~,~] = svd(B); 
        if 0
        for i = 1:d
                A1(:,i) = A1(:,i)/norm(A1(:,i));
        end
        for i = 1:d
                A2(:,i) = A2(:,i)/norm(A2(:,i));
        end
        for i = 1:d
                A3(:,i) = A3(:,i)/norm(A3(:,i));
        end
        for i = 1:d
                A4(:,i) = A4(:,i)/norm(A4(:,i));
        end
        end
        

        for ind = 1: num
            phi = randn(d,1);
            psi = randn(d,1);
           
            gammaA1(d,ind_A,ind) = 1/calGamma(A1,phi,psi);
            gammaA2(d,ind_A,ind) = 1/calGamma(A2,phi,psi);
            gammaA3(d,ind_A,ind) = 1/calGamma(A3,phi,psi);
            gammaA4(d,ind_A,ind) = 1/calGamma(A4,phi,psi);
            gamma_R(d,ind_A,ind) = 1/calGamma(R,phi,psi);
        end
    end
end

gammaA1_mean_mid = zeros(d_max,num_A);
gammaA1_min_mid = zeros(d_max,num_A);
gammaA1_max_mid = zeros(d_max,num_A);
gammaA2_mean_mid = zeros(d_max,num_A);
gammaA2_min_mid = zeros(d_max,num_A);
gammaA2_max_mid = zeros(d_max,num_A);
gammaA3_mean_mid = zeros(d_max,num_A);
gammaA3_min_mid = zeros(d_max,num_A);
gammaA3_max_mid = zeros(d_max,num_A);
gammaA4_mean_mid = zeros(d_max,num_A);
gammaA4_min_mid = zeros(d_max,num_A);
gammaA4_max_mid = zeros(d_max,num_A);

gammaR_max_mid = zeros(d_max,num_A);
gammaR_min_mid = zeros(d_max,num_A);
gammaR_mean_mid = zeros(d_max,num_A);

for ind = 2:d_max
    for ind2 = 1: num_A
        gammaA1_mean_mid(ind,ind2) = mean(gammaA1(ind,ind2,:));
        gammaA1_max_mid(ind,ind2) = max(gammaA1(ind,ind2,:));
        gammaA1_min_mid(ind,ind2) = min(gammaA1(ind,ind2,:));
        gammaA2_mean_mid(ind,ind2) = mean(gammaA2(ind,ind2,:));
        gammaA2_max_mid(ind,ind2) = max(gammaA2(ind,ind2,:));
        gammaA2_min_mid(ind,ind2) = min(gammaA2(ind,ind2,:));
        gammaA3_mean_mid(ind,ind2) = mean(gammaA3(ind,ind2,:));
        gammaA3_max_mid(ind,ind2) = max(gammaA3(ind,ind2,:));
        gammaA3_min_mid(ind,ind2) = min(gammaA3(ind,ind2,:));
        gammaA4_mean_mid(ind,ind2) = mean(gammaA4(ind,ind2,:));
        gammaA4_max_mid(ind,ind2) = max(gammaA4(ind,ind2,:));
        gammaA4_min_mid(ind,ind2) = min(gammaA4(ind,ind2,:));        
        
        gammaR_mean_mid(ind,ind2) = mean(gamma_R(ind,ind2,:));
        gammaR_min_mid(ind,ind2) = min(gamma_R(ind,ind2,:));
        gammaR_max_mid(ind,ind2) = max(gamma_R(ind,ind2,:));
    end
end

gammaA1_mean = zeros(d_max,1);
gammaA1_max = zeros(d_max,1);
gammaA1_min = zeros(d_max,1);
gammaA2_mean = zeros(d_max,1);
gammaA2_max = zeros(d_max,1);
gammaA2_min = zeros(d_max,1);
gammaA3_mean = zeros(d_max,1);
gammaA3_max = zeros(d_max,1);
gammaA3_min = zeros(d_max,1);
gammaA4_mean = zeros(d_max,1);
gammaA4_max = zeros(d_max,1);
gammaA4_min = zeros(d_max,1);

gammaR_mean = zeros(d_max,1);
gammaR_max = zeros(d_max,1);
gammaR_min = zeros(d_max,1);

for ind = 2:d_max
        gammaA1_mean(ind,1) = mean(gammaA1_min_mid(ind,:));
        gammaA1_min(ind,1) = min(gammaA1_min_mid(ind,:));
        gammaA1_max(ind,1) = max(gammaA1_min_mid(ind,:));
        
        gammaA2_mean(ind,1) = mean(gammaA2_min_mid(ind,:));
        gammaA2_min(ind,1) = min(gammaA2_min_mid(ind,:));
        gammaA2_max(ind,1) = max(gammaA2_min_mid(ind,:));    
        
        gammaA3_mean(ind,1) = mean(gammaA3_min_mid(ind,:));
        gammaA3_min(ind,1) = min(gammaA3_min_mid(ind,:));
        gammaA3_max(ind,1) = max(gammaA3_min_mid(ind,:));
        
        gammaA4_mean(ind,1) = mean(gammaA4_min_mid(ind,:));
        gammaA4_min(ind,1) = min(gammaA4_min_mid(ind,:));
        gammaA4_max(ind,1) = max(gammaA4_min_mid(ind,:));
        
        gammaR_mean(ind,1) = mean(gammaR_min_mid(ind,:));
        gammaR_min(ind,1) = min(gammaR_min_mid(ind,:));
        gammaR_max(ind,1) = max(gammaR_min_mid(ind,:));
end


figure;
Xaxis = 5:10:d_max;

plot(Xaxis,log(gammaR_mean(Xaxis,1)),'-r','LineWidth',2);
hold all;
plot(Xaxis,log(gammaA1_mean(Xaxis,1)),'*b','LineWidth',2);
hold all;
plot(Xaxis,log(gammaA2_mean(Xaxis,1)),'-.b','LineWidth',2);
hold all;
plot(Xaxis,log(gammaA3_mean(Xaxis,1)),'-b','LineWidth',2);
hold all;
plot(Xaxis,log(gammaA4_mean(Xaxis,1)),'-k','LineWidth',2);

h_legend = legend('1 / \gamma_{R}', '1 / \gamma_{A1}', '1/ \gamma_{A2}', '1/ \gamma_{A3}', '1/\gamma_{A4}', 'Fontsize', 14);
%set(h_legend, 'Fontsize', 14);
%title('Minimal Spacing (200 repetations )');
xlabel('Dimension', 'Fontsize', 14);
ylabel('log( 1/ Minimal Spacing)', 'Fontsize', 14)

%figure;
%boxplot([log(gamma_min_mid(7,:)); log(gamma_R_min_mid(7,:))]');