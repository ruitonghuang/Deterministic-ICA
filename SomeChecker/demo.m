hdis = 0;
vdis = -1.4;
x = 0:0.005:15;
A = 0.5 - mod(floor(x),2);
B = cos(x);
figure; 
subplot(10,1,1); plot(x, A,'black','LineWidth',1.5);         % plot A
ylim([-1.5 1.5]);
y = ylabel('\fontsize{20} Hidden {\fontsize{70}\{ }','Rotation', 0);
set(y, 'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);
subplot(10,1,2); plot(x, B, 'black','LineWidth',1.5);    % plot B
ylim([-1.5 1.5]);

M1 = A - 2*B;                  % mixing 1
M2 = 2.6*A-5.1*B;            % mixing 2
%figure;
subplot(10,1,3); plot(x, M1, 'r','LineWidth',1.5);      % plot mixing 1
y = ylabel('\fontsize{20} Observed {\fontsize{70}\{ }','Rotation', 0);
set(y,'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);
subplot(10,1,4); plot(x, M2, 'r','LineWidth',1.5); % plot mixing 2

hdis = 0.005;
vdis = -1.4;
%seed = rng;
rng(seed)
ICA.opt_type = 'fastICA';
d = 2;
[~,W] = estimate_ICA([M1;M2],ICA,d); % compute and plot unminxing using fastICA	
sig = W*[M1;M2];
sig = rescale(sig);
subplot(10,1,5); plot(x, sig(1,:),'b','LineWidth',1.5);
ylim([-1.5 1.5]);
y = ylabel('\fontsize{17} FastICA Result{\fontsize{70}\{ }','Rotation', 0);
set(y, 'HorizontalAlignment','right', 'Units', 'Normalized', 'Position', [hdis, vdis, 0]);
subplot(10,1,6); plot(x, sig(2,:),'b','LineWidth',1.5);
ylim([-1.5 1.5]);
%y = ylabel('FastICA 2', 'FontSize',20,'FontWeight','bold','Rotation', 0);
%set(y, 'Units', 'Normalized', 'Position', [hdis, vdis, 0]);

A_DICA = DICA([M1;M2]);
sig_DICA = A_DICA\[M1;M2];
sig_DICA = rescale(sig_DICA);
subplot(10,1,7); plot(x,sig_DICA(2,:),'b','LineWidth',1.5);
ylim([-1.5 1.5]);
y = ylabel('\fontsize{17} DICA Result{\fontsize{70}\{ }','Rotation', 0);
set(y, 'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);
subplot(10,1,8); plot(x, sig_DICA(1,:),'b','LineWidth',1.5);
ylim([-1.5 1.5]);
%y = ylabel('DICA 2', 'FontSize',20,'FontWeight','bold', 'Rotation', 0);
%set(y, 'Units', 'Normalized', 'Position', [hdis, vdis, 0]);

A_HKICA = DHsu([M1;M2]);
sig_HKICA = A_HKICA\[M1;M2];
sig_HKICA = rescale(sig_HKICA);
subplot(10,1,9); plot(x,sig_HKICA(2,:),'b','LineWidth',1.5);
ylim([-1.5 1.5]);
y = ylabel('\fontsize{17} HKICA Result{\fontsize{70}\{ }','Rotation', 0);
set(y, 'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);
subplot(10,1,10); plot(x, sig_HKICA(1,:),'b','LineWidth',1.5);
ylim([-1.5 1.5]);
%y = ylabel('HKICA 2', 'FontSize',20,'FontWeight','bold', 'Rotation', 0);
%set(y, 'Units', 'Normalized', 'Position', [hdis, vdis, 0]);