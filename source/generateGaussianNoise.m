function eps = generateGaussianNoise(num_sample, mu, sigma)
d = size(mu,1);
eps = zeros(d,num_sample);
for i = 1: num_sample
    eps(:,i) = mvnrnd(mu,sigma);
end
