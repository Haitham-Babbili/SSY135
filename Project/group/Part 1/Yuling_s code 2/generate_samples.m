function X = generate_samples(n,Gp)
% a = 1/sqrt(2);
a = n^2/(norm(Gp)^2)/2;
X = (randn(n,1) + 1i*randn(n,1))*sqrt(a);

end