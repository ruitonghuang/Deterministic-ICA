load('error1.mat');
load('error2.mat');
load('error3.mat');
load('error4.mat');

errortype = 3;

figure;
subplot(1,4,1);
    plot(Xaixs, log(error_FICA_mean_min1(sam,:,errortype)),':r', 'LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsuSymm_mean_min1(sam,:,errortype)),'--b','LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsu_mean_min1(sam,:,errortype)),':b','LineWidth',2);
    hold all;
    plot( Xaixs, log(error_DICASimp_mean_min1(sam,:,errortype)),'--k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICASymm_mean_min1(sam,:,errortype)),'-.k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICA_mean_min1(sam,:,errortype)),':k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_tt_mean_min1*ones(1,length(Xaixs))),'-k', 'LineWidth',2);
    hold on;
    xlabel('Noise\_ratio c');
    title('Matrix Type 1.');
    ylabel('log(Reconstruction Error)');
    %h_legend = legend('FICA', 'HKICA\_Symm', 'HKICA\_Re', 'DICA\_Int', 'DICA\_Symm', 'DICA\_Re', 'Random');

    subplot(1,4,2);
    plot(Xaixs, log(error_FICA_mean_min2(sam,:,errortype)),':r', 'LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsuSymm_mean_min2(sam,:,errortype)),'--b','LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsu_mean_min2(sam,:,errortype)),':b','LineWidth',2);
    hold all;
    plot( Xaixs, log(error_DICASimp_mean_min2(sam,:,errortype)),'--k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICASymm_mean_min2(sam,:,errortype)),'-.k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICA_mean_min2(sam,:,errortype)),':k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_tt_mean_min2*ones(1,length(Xaixs))),'-k', 'LineWidth',2);
    hold on;
    xlabel('Noise\_ratio c');
    ylabel('log(Reconstruction Error)');
    title('Matrix Type 2.');
    
    subplot(1,4,3);
    plot(Xaixs, log(error_FICA_mean_min3(sam,:,errortype)),':r', 'LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsuSymm_mean_min3(sam,:,errortype)),'--b','LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsu_mean_min3(sam,:,errortype)),':b','LineWidth',2);
    hold all;
    plot( Xaixs, log(error_DICASimp_mean_min3(sam,:,errortype)),'--k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICASymm_mean_min3(sam,:,errortype)),'-.k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICA_mean_min3(sam,:,errortype)),':k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_tt_mean_min3*ones(1,length(Xaixs))),'-k', 'LineWidth',2);
    hold on;
    xlabel('Noise\_ratio c');
    title('Matrix Type 3.');
    ylabel('log(Reconstruction Error)');
    %h_legend = legend('FICA', 'HKICA\_Symm', 'HKICA\_Re', 'DICA\_Int', 'DICA\_Symm', 'DICA\_Re', 'Random');

    subplot(1,4,4);
    plot(Xaixs, log(error_FICA_mean_min4(sam,:,errortype)),':r', 'LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsuSymm_mean_min4(sam,:,errortype)),'--b','LineWidth',2);
    hold all;
    plot(Xaixs, log(error_DHsu_mean_min4(sam,:,errortype)),':b','LineWidth',2);
    hold all;
    plot( Xaixs, log(error_DICASimp_mean_min4(sam,:,errortype)),'--k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICASymm_mean_min1(sam,:,errortype)),'-.k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_DICA_mean_min4(sam,:,errortype)),':k', 'LineWidth',2);
    hold on;
    plot( Xaixs, log(error_tt_mean_min4*ones(1,length(Xaixs))),'-k', 'LineWidth',2);
    hold on;
    xlabel('Noise\_ratio c');
    ylabel('log(Reconstruction Error)');
    title('Matrix Type 4.');
    
    
    
    h_legend = legend('FICA', 'HKICA\_Symm', 'HKICA\_Re', 'DICA\_Int', 'DICA\_Symm', 'DICA\_Re', 'Random');
