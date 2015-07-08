clear;
d = 6;
noise_ratio_list = 0:0.3:0.6;
sample_size_list = 1000*2.^(0:7); 
num_trail = 150;
num_phi = 3; 
%num_sample_test = 5000;
cohe = 0.005;
run_DICA = 1;
run_DHsu = 1;

run_FICA = 1;
run_FPCA = 1;

%phi_all = zeros(d,length(sample_size_list), length(noise_ratio_list), num_trail,num_phi);
%psi1_all = zeros(d,length(sample_size_list), length(noise_ratio_list), num_trail,num_phi);
%psi2_all = zeros(d,length(sample_size_list), length(noise_ratio_list), num_trail,num_phi);

%x_all = zeros(d, 405*d^2, length(sample_size_list), length(noise_ratio_list), num_trail);
%y_all = zeros(d, 405*d^2, length(sample_size_list), length(noise_ratio_list), num_trail);
%x_all_test = zeros(d,num_sample_test, length(sample_size_list), length(noise_ratio_list), num_trail);
%y_all_test = zeros(d,num_sample_test, length(sample_size_list), length(noise_ratio_list), num_trail);

%Gamma = zeros(length(sample_size_list), length(noise_ratio_list), num_trail);
error_FICA = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_FPCA = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DICA = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DICASimp = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DICASymm = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DICA_recur = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DICASimp_recur = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DICASymm_recur = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DHsu = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DHsuSymm = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DHsu_recur = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);
error_DHsuSymm_recur = 100*ones(length(sample_size_list), length(noise_ratio_list), num_trail,num_phi,3);

time_FICA = 0;
time_FPCA = 0;
time_DICA = 0;
time_DICASimp = 0;
time_DICASymm = 0;
time_DHsu = 0;
time_DHsuSymm=0;
time_DICA_recur = 0;
time_DICASimp_recur = 0;
time_DICASymm_recur = 0;
time_DHsu_recur = 0;
time_DHsuSymm_recur=0;


for ind1 = 1: length(sample_size_list)
    num_sample = sample_size_list(ind1);  
    for ind2 = 1: length(noise_ratio_list)
        noise_ratio = noise_ratio_list(ind2);         
        flag = 0;
        ind3 = 1;
        while ind3 <= num_trail
            TT = randn(d); %randn(d,1)*ones(1,d)+cohe*
            
            [TT,~,~] = svd(TT);
            A  = TT;
            %A = randn(d,1);
            %A = 2*A/norm(A,2);            
            %A = A*ones(1,d);
            %A = A+cohe*TT; %
            %
            %A = randn(d);
            %while min(svd(A))>0.001%max(svd(A))/min(svd(A))>20
            %    A = randn(d,1)*ones(1,d)+cohe*randn(d);
            %end
            %A(:,2) = A(:,3)+ 0.01*randn(d,1)
            A = A/diag([norm(A(:,1)),norm(A(:,2)),norm(A(:,3)),norm(A(:,4)),norm(A(:,5)),norm(A(:,6))]);            
            x = zeros(d,num_sample);
            for i = 1:d
                x(i,:) = bpsk(num_sample,0,i); %rand(1,num_sample) -0.5;%
            end
            
            mu = zeros(d,1);
            sigmahalf = randn(d);
            sigma = sigmahalf*sigmahalf';
            eps = generateGaussianNoise(num_sample, mu, sigma);
            
            y = A*x + noise_ratio*eps;            
            %x_test = zeros(d,num_sample_test);
            %init = num_sample;
            %for i = 1:d
            %    x_test(i,:) =  bpsk(num_sample_test,init,i); %rand(1,num_sample_test) - 0.5;%
            %end
            %y_test = A*x_test;
     
            for j = 1:num_phi 
                TMP = eye(d);
                [V,y_white] = whitening(y);
                
                if(run_DICA)
                tic;
                Ahat_DICA = DICA(y);
                toc;
                time = toc;
                time_DICA = time_DICA+time;
                error1 = calError(Ahat_DICA,A);
                %error2 = calError2(Ahat_DICA,A,y);
                %error3 = calError_entropy(Ahat_DICA,y_test);
                error_DICA(ind1,ind2,ind3,j,1) = error1;
                %error_DICA(ind1,ind2,ind3,j,2) = error2;
                %error_DICA(ind1,ind2,ind3,j,3) = error3; 
                
                
                tic;
                Ahat_DICA_recur = V*DICA_recur(y_white,TMP);
                toc;
                time = toc;
                time_DICA_recur = time_DICA_recur+time;
                error1 = calError(Ahat_DICA_recur,A);
                %error2 = calError2(Ahat_DICA_recur,A,y);
                %error3 = calError_entropy(Ahat_DICA_recur,y_test);
                error_DICA_recur(ind1,ind2,ind3,j,1) = error1;
                %error_DICA_recur(ind1,ind2,ind3,j,2) = error2;
                %error_DICA_recur(ind1,ind2,ind3,j,3) = error3; 
                
                
                    
                %if 0
                %tic;
                %Ahat_DICASymm = DICASymm(y);
                %toc;
                %time = toc;
                %time_DICASymm = time_DICASymm+time;
                %error1 = calError(Ahat_DICASymm,A);
                %error2 = calError2(Ahat_DICASymm,A,y);
                %error3 = calError_entropy(Ahat_DICASymm,y_test);
                %error_DICASymm(ind1,ind2,ind3,j,1) = error1;
                %error_DICASymm(ind1,ind2,ind3,j,2) = error2;
                %error_DICASymm(ind1,ind2,ind3,j,3) = error3; 
                
                
                %tic;
                %Ahat_DICASymm_recur = V*DICASymm_recur(y_white,TMP);
                %toc;
                %time = toc;
                %time_DICASymm_recur = time_DICASymm_recur+time;
                %error1 = calError(Ahat_DICASymm_recur,A);
                %error2 = calError2(Ahat_DICASymm_recur,A,y);
                %error3 = calError_entropy(Ahat_DICASymm_recur,y_test);
                %error_DICASymm_recur(ind1,ind2,ind3,j,1) = error1;
                %error_DICASymm_recur(ind1,ind2,ind3,j,2) = error2;
                %error_DICASymm_recur(ind1,ind2,ind3,j,3) = error3; 
                %end
                
                
                tic;
                Ahat_DICASimp = DICASimp(y);
                toc;
                time = toc;
                time_DICASimp = time_DICASimp+time;
                error1 = calError(Ahat_DICASimp,A);
                %error2 = calError2(Ahat_DICA,A,y);
                %error3 = calError_entropy(Ahat_DICASimp,y_test);
                error_DICASimp(ind1,ind2,ind3,j,1) = error1;
                %error_DICASimp(ind1,ind2,ind3,j,2) = error2;
                %error_DICASimp(ind1,ind2,ind3,j,3) = error3;
                
                
                tic;
                Ahat_DICASimp_recur = V*DICASimp_recur(y_white,TMP);
                toc;
                time = toc;
                time_DICASimp_recur = time_DICASimp_recur+time;
                error1 = calError(Ahat_DICASimp_recur,A);
                %error2 = calError2(Ahat_DICA_recur,A,y);
                %error3 = calError_entropy(Ahat_DICASimp_recur,y_test);
                error_DICASimp_recur(ind1,ind2,ind3,j,1) = error1;
                %error_DICASimp_recur(ind1,ind2,ind3,j,2) = error2;
                %error_DICASimp_recur(ind1,ind2,ind3,j,3) = error3;
                
                end
                
                
                if(run_DHsu)
                tic;
                Ahat_DHsu = DHsu(y);
                toc;
                time = toc;
                time_DHsu = time_DHsu+time;
                error1 = calError(Ahat_DHsu,A);
                %error2 = calError2(Ahat_DHsu,A,y);
                %error3 = calError_entropy(Ahat_DHsu,y_test);
                error_DHsu(ind1,ind2,ind3,j,1) = error1;
                %error_DHsu(ind1,ind2,ind3,j,2) = error2;
                %error_DHsu(ind1,ind2,ind3,j,3) = error3;
                
                
                tic;
                Ahat_DHsu_recur = V*DHsu_recur(y_white,TMP);
                toc;
                time = toc;
                time_DHsu_recur = time_DHsu_recur+time;
                error1 = calError(Ahat_DHsu_recur,A);
                %error2 = calError2(Ahat_DHsu_recur,A,y);
                %error3 = calError_entropy(Ahat_DHsu_recur,y_test);
                error_DHsu_recur(ind1,ind2,ind3,j,1) = error1;
                %error_DHsu_recur(ind1,ind2,ind3,j,2) = error2;
                %error_DHsu_recur(ind1,ind2,ind3,j,3) = error3;
                
                
                
                %if 0
                %tic;
                %Ahat_DHsuSymm = DHsuSymm(y);
                %toc;
                %time = toc;
                %time_DHsuSymm = time_DHsuSymm+time;
                %error1 = calError(Ahat_DHsuSymm,A);
                %error2 = calError2(Ahat_DHsuSymm,A,y);
                %error3 = calError_entropy(Ahat_DHsuSymm,y_test);
                %error_DHsuSymm(ind1,ind2,ind3,j,1) = error1;
                %error_DHsuSymm(ind1,ind2,ind3,j,2) = error2;
                %error_DHsuSymm(ind1,ind2,ind3,j,3) = error3;
                
                
                %tic;
                %Ahat_DHsuSymm_recur = V*DHsuSymm_recur(y_white,TMP);
                %toc;
                %time = toc;
                %time_DHsuSymm_recur = time_DHsuSymm_recur+time;
                %error1 = calError(Ahat_DHsuSymm_recur,A);
                %error2 = calError2(Ahat_DHsuSymm_recur,A,y);
                %error3 = calError_entropy(Ahat_DHsuSymm_recur,y_test);
                %error_DHsuSymm_recur(ind1,ind2,ind3,j,1) = error1;
                %error_DHsuSymm_recur(ind1,ind2,ind3,j,2) = error2;
                %error_DHsuSymm_recur(ind1,ind2,ind3,j,3) = error3;
                %end
                end
                
                %%%%%%%%%%%%%%%%%%% run FICA    
                if (run_FICA)
                tic;
                Ahat_FICA = FICA(y);
                toc;
                time = toc;
                time_FICA = time_FICA+time;
                error1 = calError(Ahat_FICA,A);
                %error2 = calError2(Ahat_FICA,A,y);
                %error3 = calError_entropy(Ahat_FICA,y_test);
                error_FICA(ind1,ind2,ind3,j,1) = error1;
                %error_FICA(ind1,ind2,ind3,j,2) = error2;
                %error_FICA(ind1,ind2,ind3,j,3) = error3;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run Fourier PCA
                if (run_FPCA)
                tic;
                dap = 1.6;
                Ahat_FPCA = V*recursiveFPCA(eye(d),y_white',dap);
                toc;
                time = toc;
                time_FPCA = time_FPCA+time;
                error1 = calError(Ahat_FPCA,A);
                %error2 = calError2(Ahat_FPCA,A,y);
                %error3 = calError_entropy(Ahat_FPCA,y_test);
                error_FPCA(ind1,ind2,ind3,j,1) = error1;
                %error_FPCA(ind1,ind2,ind3,j,2) = error2;
                %error_FPCA(ind1,ind2,ind3,j,3) = error3;
                
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            end
            disp({'sample_size',ind1, 'Noise:',  ind2, 'trail:', ind3});            
            ind3 = ind3+1;
            %if flag == 1
             %   ind3 = ind3 -1;
            %    flag = 0;
           % end
        end
    end
end
save('bpsk_samplesize_R.mat');

%if error_FICA(ind1,ind2,ind3,1) > error1
                %    error_FICA(ind1,ind2,ind3,1) = error1;
                %end
                %if error_FICA(ind1,ind2,ind3,2) > error2
                %    error_FICA(ind1,ind2,ind3,2) = error2;
                %end
                %if error_FICA(ind1,ind2,ind3,3) > error3
                %    error_FICA(ind1,ind2,ind3,3) = error3;
                %end




