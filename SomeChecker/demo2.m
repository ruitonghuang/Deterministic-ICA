load('demo2-noisy.mat');

%seed = rng;
rng(seed)
noise_ratio = 0.4;
d = 6;
num_sample = 20000;

TT = randn(d);
[TT,~,~] = svd(TT);
A  = TT;
x = zeros(d,num_sample);
for i = 1:d
    x(i,:) = bpsk(num_sample,0,i); 
end
mu = zeros(d,1);
sigmahalf = randn(d);
sigma = sigmahalf*sigmahalf';
eps = generateGaussianNoise(num_sample, mu, sigma);
y = A*x + noise_ratio*eps; 
Ahat_DICA = DICA(y);

hdis = -0.01;
vdis = 3.5;

sig = x;
for i = 1:6
    sig(i,:) = sig(i,:)/max(abs(sig(i,:))); 
end
figure; 
subplot(6,1,1); 
plot(1:50,sig(1,1:50),'black','LineWidth',1.5);
subplot(6,1,2);
plot(1:50,sig(2,1:50),'black','LineWidth',1.5);
subplot(6,1,3);
plot(1:50,sig(3,1:50),'black','LineWidth',1.5);
subplot(6,1,4);
plot(1:50,sig(4,1:50),'black','LineWidth',1.5);
subplot(6,1,5);
plot(1:50,sig(5,1:50),'black','LineWidth',1.5);
subplot(6,1,6);
plot(1:50,sig(6,1:50),'black','LineWidth',1.5);
qq = ylabel('\fontsize{20} Hidden','Rotation', 0);
set(qq, 'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);

sig = y;
for i = 1:6
    sig(i,:) = sig(i,:)/max(abs(sig(i,:))); 
end
figure; 
subplot(6,1,1); 
plot(1:50,sig(1,1:50),'r','LineWidth',1.5);
subplot(6,1,2);
plot(1:50,sig(2,1:50),'r','LineWidth',1.5);
subplot(6,1,3);
plot(1:50,sig(3,1:50),'r','LineWidth',1.5);
subplot(6,1,4);
plot(1:50,sig(4,1:50),'r','LineWidth',1.5);
subplot(6,1,5);
plot(1:50,sig(5,1:50),'r','LineWidth',1.5);
subplot(6,1,6);
plot(1:50,sig(6,1:50),'r','LineWidth',1.5);
qq = ylabel('\fontsize{20} Observed','Rotation', 0);
set(qq, 'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);

sig = Ahat_DICA\A*x;
for i = 1:6
    sig(i,:) = sig(i,:)/max(abs(sig(i,:))); 
end
figure; 
subplot(6,1,1); 
plot(1:50,sig(4,1:50),'b','LineWidth',1.5);
subplot(6,1,2);
plot(1:50,-sig(1,1:50),'b','LineWidth',1.5);
subplot(6,1,3);
plot(1:50,sig(2,1:50),'b','LineWidth',1.5);
subplot(6,1,4);
plot(1:50,sig(5,1:50),'b','LineWidth',1.5);
subplot(6,1,5);
plot(1:50,sig(6,1:50),'b','LineWidth',1.5);
subplot(6,1,6);
plot(1:50,sig(3,1:50),'b','LineWidth',1.5);
qq = ylabel('\fontsize{20} Reconstructed','Rotation', 0);
set(qq, 'HorizontalAlignment','right','Units', 'Normalized', 'Position', [hdis, vdis, 0]);

%save('demo2-noisy.mat');

