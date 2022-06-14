function c = generate_rayleigh_channel(X,G)

C = X.*G;
c = ifft(C);
%c_norm = c/sqrt(mean(abs(c).^2));
%c_norm = generate_normalized(c);

end