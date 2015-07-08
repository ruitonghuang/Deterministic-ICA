run_recur = 1;

error_FICA_mean_intermid = mean(error_FICA,4);
error_FICA_min_intermid = min(error_FICA,[],4);
error_FICA_max_intermid = max(error_FICA,[],4);
error_DICA_mean_intermid = mean(error_DICA,4);
error_DICA_min_intermid = min(error_DICA,[],4);
error_DICA_max_intermid = max(error_DICA,[],4);
error_DICASymm_mean_intermid = mean(error_DICASymm,4);
error_DICASymm_min_intermid = min(error_DICASymm,[],4);
error_DICASymm_max_intermid = max(error_DICASymm,[],4);
error_DICASimp_mean_intermid = mean(error_DICASimp,4);
error_DICASimp_min_intermid = min(error_DICASimp,[],4);
error_DICASimp_max_intermid = max(error_DICASimp,[],4);
error_DHsu_mean_intermid = mean(error_DHsu,4);
error_DHsu_min_intermid = min(error_DHsu,[],4);
error_DHsu_max_intermid = max(error_DHsu,[],4);
error_DHsuSymm_mean_intermid = mean(error_DHsuSymm,4);
error_DHsuSymm_min_intermid = min(error_DHsuSymm,[],4);
error_DHsuSymm_max_intermid = max(error_DHsuSymm,[],4);

error_FICA_mean_mean = mean(error_FICA_mean_intermid,3);
error_FICA_mean_max = mean(error_FICA_max_intermid,3);
error_FICA_mean_min = mean(error_FICA_min_intermid,3);
error_DICA_mean_mean = mean(error_DICA_mean_intermid,3);
error_DICA_mean_max = mean(error_DICA_max_intermid,3);
error_DICA_mean_min = mean(error_DICA_min_intermid,3);
error_DICASymm_mean_mean = mean(error_DICASymm_mean_intermid,3);
error_DICASymm_mean_max = mean(error_DICASymm_max_intermid,3);
error_DICASymm_mean_min = mean(error_DICASymm_min_intermid,3);
error_DICASimp_mean_mean = mean(error_DICASimp_mean_intermid,3);
error_DICASimp_mean_max = mean(error_DICASimp_max_intermid,3);
error_DICASimp_mean_min = mean(error_DICASimp_min_intermid,3);
error_DHsu_mean_mean = mean(error_DHsu_mean_intermid,3);
error_DHsu_mean_max = mean(error_DHsu_max_intermid,3);
error_DHsu_mean_min = mean(error_DHsu_min_intermid,3);
error_DHsuSymm_mean_mean = mean(error_DHsuSymm_mean_intermid,3);
error_DHsuSymm_mean_max = mean(error_DHsuSymm_max_intermid,3);
error_DHsuSymm_mean_min = mean(error_DHsuSymm_min_intermid,3);

error_FICA_max_mean = max(error_FICA_mean_intermid,[],3);
error_FICA_max_max = max(error_FICA_max_intermid,[],3);
error_FICA_max_min = max(error_FICA_min_intermid,[],3);
error_DICASimp_max_mean = max(error_DICASimp_mean_intermid,[],3);
error_DICASimp_max_max = max(error_DICASimp_max_intermid,[],3);
error_DICASimp_max_min = max(error_DICASimp_min_intermid,[],3);
error_DICASymm_max_mean = max(error_DICASymm_mean_intermid,[],3);
error_DICASymm_max_max = max(error_DICASymm_max_intermid,[],3);
error_DICASymm_max_min = max(error_DICASymm_min_intermid,[],3);
error_DICA_max_mean = max(error_DICA_mean_intermid,[],3);
error_DICA_max_max = max(error_DICA_max_intermid,[],3);
error_DICA_max_min = max(error_DICA_min_intermid,[],3);
error_DHsu_max_mean = max(error_DHsu_mean_intermid,[],3);
error_DHsu_max_max = max(error_DHsu_max_intermid,[],3);
error_DHsu_max_min = max(error_DHsu_min_intermid,[],3);
error_DHsuSymm_max_mean = max(error_DHsuSymm_mean_intermid,[],3);
error_DHsuSymm_max_max = max(error_DHsuSymm_max_intermid,[],3);
error_DHsuSymm_max_min = max(error_DHsuSymm_min_intermid,[],3);

if (run_recur==1)
error_FPCA_mean_intermid = mean(error_FPCA,4);
error_FPCA_min_intermid = min(error_FPCA,[],4);
error_FPCA_max_intermid = max(error_FPCA,[],4);
error_DICA_recur_mean_intermid = mean(error_DICA_recur,4);
error_DICA_recur_min_intermid = min(error_DICA_recur,[],4);
error_DICA_recur_max_intermid = max(error_DICA_recur,[],4);
error_DICASymm_recur_mean_intermid = mean(error_DICASymm_recur,4);
error_DICASymm_recur_min_intermid = min(error_DICASymm_recur,[],4);
error_DICASymm_recur_max_intermid = max(error_DICASymm_recur,[],4);
error_DICASimp_recur_mean_intermid = mean(error_DICASimp_recur,4);
error_DICASimp_recur_min_intermid = min(error_DICASimp_recur,[],4);
error_DICASimp_recur_max_intermid = max(error_DICASimp_recur,[],4);
error_DHsu_recur_mean_intermid = mean(error_DHsu_recur,4);
error_DHsu_recur_min_intermid = min(error_DHsu_recur,[],4);
error_DHsu_recur_max_intermid = max(error_DHsu_recur,[],4);
error_DHsuSymm_recur_mean_intermid = mean(error_DHsuSymm_recur,4);
error_DHsuSymm_recur_min_intermid = min(error_DHsuSymm_recur,[],4);
error_DHsuSymm_recur_max_intermid = max(error_DHsuSymm_recur,[],4);

error_FPCA_mean_mean = mean(error_FPCA_mean_intermid,3);
error_FPCA_mean_max = mean(error_FPCA_max_intermid,3);
error_FPCA_mean_min = mean(error_FPCA_min_intermid,3);
error_DICA_recur_mean_mean = mean(error_DICA_recur_mean_intermid,3);
error_DICA_recur_mean_max = mean(error_DICA_recur_max_intermid,3);
error_DICA_recur_mean_min = mean(error_DICA_recur_min_intermid,3);
error_DICASymm_recur_mean_mean = mean(error_DICASymm_recur_mean_intermid,3);
error_DICASymm_recur_mean_max = mean(error_DICASymm_recur_max_intermid,3);
error_DICASymm_recur_mean_min = mean(error_DICASymm_recur_min_intermid,3);
error_DICASimp_recur_mean_mean = mean(error_DICASimp_recur_mean_intermid,3);
error_DICASimp_recur_mean_max = mean(error_DICASimp_recur_max_intermid,3);
error_DICASimp_recur_mean_min = mean(error_DICASimp_recur_min_intermid,3);
error_DHsu_recur_mean_mean = mean(error_DHsu_recur_mean_intermid,3);
error_DHsu_recur_mean_max = mean(error_DHsu_recur_max_intermid,3);
error_DHsu_recur_mean_min = mean(error_DHsu_recur_min_intermid,3);
error_DHsuSymm_recur_mean_mean = mean(error_DHsuSymm_recur_mean_intermid,3);
error_DHsuSymm_recur_mean_max = mean(error_DHsuSymm_recur_max_intermid,3);
error_DHsuSymm_recur_mean_min = mean(error_DHsuSymm_recur_min_intermid,3);

error_FPCA_max_mean = max(error_FPCA_mean_intermid,[],3);
error_FPCA_max_max = max(error_FPCA_max_intermid,[],3);
error_FPCA_max_min = max(error_FPCA_min_intermid,[],3);
error_DICASymm_recur_max_mean = max(error_DICASymm_recur_mean_intermid,[],3);
error_DICASymm_recur_max_max = max(error_DICASymm_recur_max_intermid,[],3);
error_DICASymm_recur_max_min = max(error_DICASymm_recur_min_intermid,[],3);
error_DICA_recur_max_mean = max(error_DICA_recur_mean_intermid,[],3);
error_DICA_recur_max_max = max(error_DICA_recur_max_intermid,[],3);
error_DICA_recur_max_min = max(error_DICA_recur_min_intermid,[],3);
error_DICASimp_recur_max_mean = max(error_DICASimp_recur_mean_intermid,[],3);
error_DICASimp_recur_max_max = max(error_DICASimp_recur_max_intermid,[],3);
error_DICASimp_recur_max_min = max(error_DICASimp_recur_min_intermid,[],3);
ror_DHsu_recur_max_mean = max(error_DHsu_recur_mean_intermid,[],3);
error_DHsu_recur_max_max = max(error_DHsu_recur_max_intermid,[],3);
error_DHsu_recur_max_min = max(error_DHsu_recur_min_intermid,[],3);
error_DHsuSymm_recur_max_mean = max(error_DHsuSymm_recur_mean_intermid,[],3);
error_DHsuSymm_recur_max_max = max(error_DHsuSymm_recur_max_intermid,[],3);
error_DHsuSymm_recur_max_min = max(error_DHsuSymm_recur_min_intermid,[],3);
end

% mean_min: the mean of 200 trials; take the best of each trial




if 1==0
    error_tt = zeros(num_trail,num_phi);
    for i = 1:num_trail 
        for j = 1:num_phi
            A = randn(d,1)*ones(1,d)+0.1*randn(d);%+0.01*eye(d);
            tt = randn(d);
            error_tt(i,j) = calError(tt,A);
        end
        %error_tt(i,1) = min(error_tt(i,:));
    end
    error_tt_mean_mean = mean(mean(error_tt,2));
    error_tt_mean_min = mean(min(error_tt,[],2));
    error_tt_mean_max = mean(max(error_tt,[],2));
    error_tt_max_mean = max(mean(error_tt,2));
    error_tt_max_min = max(min(error_tt,[],2));
errortype = 1;
u = error_FICA_min_intermid(1,1,:,1,1);
u = [u, error_DHsu_min_intermid(1,1,:,1,1)];
u = [u, error_FPCA_min_intermid(1,1,:,1,1)];
u = [u, error_DICA_min_intermid(1,1,:,1,1)];
u = [u, error_DICASimp_min_intermid(1,1,:,1,1)];
u = reshape(u,[num_trail,5]);

figure;
boxplot(u, 'Labels', {'FastICA', 'HKICA', 'FPCA', 'DICA', 'MDICA'});

Xaixs = log10(sample_size_list);
%Xaixs = noise_ratio_list;

sam = 1:8;
noise = 3;

    figure;
    %errortype = 1;
    plot(Xaixs, error_FICA_mean_min(sam,noise,errortype),':r', 'LineWidth',2);
    hold all;
    %plot(Xaixs, error_DHsuSymm_mean_min(sam,:,errortype),'--b','LineWidth',2);
    %hold all;
    plot(Xaixs, error_DHsu_mean_min(sam,noise,errortype),':b','LineWidth',2);
    hold all;
    plot( Xaixs, error_DICASimp_mean_min(sam,noise,errortype),'--k', 'LineWidth',2);
    hold on;
    %plot( Xaixs, error_DICASymm_mean_min(sam,:,errortype),'-.k', 'LineWidth',2);
    %hold on;
    plot( Xaixs, error_DICA_mean_min(sam,noise,errortype),':k', 'LineWidth',2);
    hold on;
    if(run_recur ==1)
        plot(Xaixs, error_FPCA_mean_min(sam,noise,errortype),'-r', 'LineWidth',2);
        hold all;
        %plot(Xaixs, error_DHsuSymm_recur_mean_min(sam,:,errortype),'-.b','LineWidth',2);
        %hold all;
        plot(Xaixs, error_DHsu_recur_mean_min(sam,noise,errortype),'*-b','LineWidth',2);
        hold all;
        %plot( Xaixs, error_DICASymm_recur_mean_min(sam,:,errortype),'-.g', 'LineWidth',2);
        %hold on;
        plot( Xaixs, error_DICA_recur_mean_min(sam,noise,errortype),':g', 'LineWidth',2);
        hold on;
        plot( Xaixs, error_DICASimp_recur_mean_min(sam,noise,errortype),'--g', 'LineWidth',2);
    hold on;
    end
    plot( Xaixs, error_tt_mean_min*ones(1,length(Xaixs)),'*-k', 'LineWidth',2);
    hold on;
    %set(h_legend, 'Fontsize', 14);
    xlabel('log_{10}(Sample sizes)','Fontsize', 20);
    ylabel('Reconstruction Error','Fontsize', 20);
    title('The Reconstruction Error of A_1','Fontsize', 20);
    h_legend = legend('FICA', 'HKICA', 'MDICA',  'DICA', 'FPCA', 'HKICA.R', 'DICA.R', 'MDICA.R');
    %h_legend = legend('FICA', 'HKICA\_Symm', 'HKICA\_Re', 'DICA\_Mod', 'DICA\_Symm', 'DICA\_Re', 'Random');
    set(h_legend, 'Fontsize', 14);
    axis([3 5.1 0 2.5]);

    if 1==0
%figure;
noise = 2;
H  =[error_FICA_mean_intermid(sam,noise,:,1,errortype); error_DHsu_mean_intermid(sam,noise,:,1,errortype);...
    error_DICA_mean_intermid(sam,noise,:,1,errortype)];
K = reshape(H,[3,num_trail]);
K = [K',min(error_tt,[],2)];
%K = log(K);
%boxplot(K, 'Labels', { 'FICA','DHsu',  'DICA','Random'});
end
    error_FICA_mean_min4 = error_FICA_mean_min;
error_DHsuSymm_mean_min4 = error_DHsuSymm_mean_min;
error_DHsu_mean_min4 = error_DHsu_mean_min;
error_DICASimp_mean_min4 = error_DICASimp_mean_min;
error_DICASymm_mean_min4 = error_DICASymm_mean_min;
error_DICA_mean_min4 = error_DICA_mean_min;
error_tt_mean_min4 = error_tt_mean_min;
%save('error1_bpsk_d6_20000_0.05_3phi_200A.mat', 'error_FICA_mean_min4', 'error_DHsuSymm_mean_min4', 'error_DHsu_mean_min4',...
 %    'error_DICASimp_mean_min4', 'error_DICASymm_mean_min4',  'error_DICA_mean_min4', 'error_tt_mean_min4');
end