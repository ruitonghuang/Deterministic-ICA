function Ahat = FICA(y,A_init)

ICA.opt_type = 'fastICA';
d = size(y,1);
if nargin == 1
    [~,What] = estimate_ICA(y,ICA,d);
elseif nargin ==2
    [~,What] = estimate_ICA(y,ICA,d,A_init);
end
Ahat = inv(What);

