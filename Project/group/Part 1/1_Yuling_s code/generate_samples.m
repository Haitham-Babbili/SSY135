function X = generate_samples(n,a)
% a = 1/sqrt(2);

X = (randn(n,1) + 1i*randn(n,1))*sqrt(a);

end