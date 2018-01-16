function [ s_noise ] = noise_sig( s, RSB )

len = length(s);

b = randn(1, len);
alpha = sqrt((10^(-RSB/10))*(sum(s.^2)/sum(b.^2))); % prod scal;
b=alpha*b;

s_noise=s+b;

end
