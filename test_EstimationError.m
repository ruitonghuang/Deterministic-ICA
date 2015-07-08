clear;
d = 3;
noise_ratio_list = 0.01:0.03:0.15;
sample_size_list = [200000]; 
num_sample = sample_size_list(1);
num_trail = 50;
%num_phi = 3; 
num_sample_test = 5000;
cohe_list = [2];
run_DICA = 1;
run_DHsu = 1;
run_FICA = 1;
run_FPCA = 1;
error_DICA = zeros(length(cohe_list),length(noise_ratio_list),num_trail,4);
error_HK = zeros(length(cohe_list),length(noise_ratio_list),num_trail,3);


% for different cohenrence 
for ind1 = 1:length(cohe_list)
    cohe = cohe_list(ind1);
    % for different noise ratio
    for ind2 = 1:length(noise_ratio_list)
        noise_ratio = noise_ratio_list(ind2);
        for ind3 = 1:num_trail
            % generate A
            TT = randn(d); %randn(d,1)*ones(1,d)+cohe*
            [TT,~,~] = svd(TT);
            A = randn(d);
            A(:,d) = A(:,d-1);
            for i = 1:d
                A(:,i) = 2*A(:,i)/norm(A(:,i));
            end
            A = A+cohe*TT; %
            
            %A = zeros(d);
            %genertate x and y
            x = zeros(d,num_sample);
            for i = 1:d
                x(i,:) = bpsk(num_sample,0,i); %rand(1,num_sample) -0.5;%
            end
            mu = zeros(d,1);
            sigmahalf = randn(d);
            sigma = sigmahalf*sigmahalf';            
            eps = generateGaussianNoise(num_sample, mu, sigma);
            y = A*x + noise_ratio*eps;
            % compute the Hession M, M1, M2, and T
            kappa = diag(calKappa(x));
            
            
            phi = randn(d,1);%[0.4455,1.16183,3.9472]';%
            D_phi = diag((phi'*A).^2);
            M_real = A*D_phi*kappa*A';
            [U_real,Sigma_real,~] =svd(M_real);
            B_real = U_real*Sigma_real^(1/2);
            R_real = (A*D_phi^(1/2)*abs(kappa)^(1/2))\B_real;

            psi1 = randn(d,1);% [2.2987, 2.4374, 0.8254]';%
            R_psi1 = diag((psi1'*R_real).^2);
            M1_real = -A*R_psi1*D_phi^(-1)*A';


            psi2 = randn(d,1);%[-0.9154,-0.5175,-0.217]'; %
            R_psi2 = diag((psi2'*R_real).^2);
            M2_real = -A*R_psi2*D_phi^(-1)*A';

            T_real = A*R_psi1*R_psi2^(-1)*A^(-1);

            % compute its estimation
            G1_phi = G1(phi,y);
            G2_phi = G2(phi,y);
            G3_phi = G3(phi,y);
            M = G1_phi-G2_phi-2*G3_phi;
            [U,D,~] = svd(M);
            X = U*D.^(1/2);
            R = (A*D_phi^(1/2)*abs(kappa)^(1/2))\X;

            R_psi1_int = diag((psi1'*R).^2);
            M1_int = A*R_psi1_int*D_phi^(-1)*A';

            G1_psi1 = G1(X'\psi1,y);
            G2_psi1 = G2(X'\psi1,y);
            G3_psi1 = G3(X'\psi1,y);
            M1 = G1_psi1-G2_psi1-2*G3_psi1;

            G1_psi2 = G1(X'\psi2,y);
            G2_psi2 = G2(X'\psi2,y);
            G3_psi2 = G3(X'\psi2,y);
            M2 = G1_psi2-G2_psi2-2*G3_psi2;

            T = M1/M2;
            % error 
            error_DICA(ind1,ind2,ind3,1) = norm(M-M_real,'fro')/norm(M_real,'fro');
            error_DICA(ind1,ind2,ind3,2) = norm(M1-M1_real,'fro')/norm(M1_real,'fro');
            error_DICA(ind1,ind2,ind3,3) = norm(M2-M2_real,'fro')/norm(M2_real,'fro');
            error_DICA(ind1,ind2,ind3,4) = norm(T-T_real,'fro')/norm(T_real,'fro');
            %-------------------------------------------------------------------------------------
            phi1 = randn(d,1);
            phi2 = randn(d,1);
            D_phi1 = diag((phi1'*A).^2);
            M_real1 = A*D_phi1*kappa*A';
            D_phi2 = diag((phi2'*A).^2);
            M_real2 = A*D_phi2*kappa*A';
            T_HKreal = M_real1/M_real2;
            
            G1_phi1 = G1(phi1,y);
            G2_phi1 = G2(phi1,y);
            G3_phi1 = G3(phi1,y);
            Mphi1 = G1_phi1-G2_phi1-2*G3_phi1;
            
            G1_phi2 = G1(phi2,y);
            G2_phi2 = G2(phi2,y);
            G3_phi2 = G3(phi2,y);
            Mphi2 = G1_phi2-G2_phi2-2*G3_phi2;
            T_HK = Mphi1/Mphi2;
            
            error_HK(ind1,ind2,ind3,1) = norm(Mphi1-M_real1,'fro')/norm(M_real,'fro');
            error_HK(ind1,ind2,ind3,2) = norm(Mphi2-M_real2,'fro')/norm(M_real2,'fro');
            error_HK(ind1,ind2,ind3,3) = norm(T_HK-T_HKreal,'fro')/norm(T_HKreal,'fro');
            
            disp({'coherence',ind1, 'Noise:',  ind2, 'trail:', ind3});   

        end
    end

end

error_DICA_mean = mean(error_DICA,3);
error_HK_mean = mean(error_HK,3);