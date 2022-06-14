function c_spectrum = spectrumMethod(f_D,T_s,N_s,K_c)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
f_s = 1/T_s;
% create the frequency vector from origin to the edge of bandwidth
f = 0:f_s/N_s:f_D; % only consider the onesided, we flip it later to create the other side
% calculate S_c according to lab PM
S_c = (1/(pi*f_D)).*(1./sqrt(1-(f./f_D).^2));
% get G as sqrt of S_c
G = sqrt(S_c);
if G(end) == Inf
    G(end) = (1/(pi*f_D)).*(1./sqrt(1-(f_D+eps/f_D).^2));
end
% initialize the periodic version G_p with zeros
G_p = zeros(1,N_s);
% store the right side of S_c in the beggining of G_p and the left side at
% the end by flipping the right side version. They should be equal.
G_p(1:length(G)) = G;
G_p(end-length(G)+1:end) = flip(G);
% Generate samples from zero-mean, complex gaussian dist
a = sqrt(N_s^2/norm(G_p)^2); % find the constant a (using Parseval's theorem)
% generate samples from a complex gaussian distribution and multiply with G
X = (randn(N_s,1) + 1i*randn(N_s,1))*a/sqrt(2); 
C = X'.*G_p; % in time-domain this is the convolution with the gaussian
c_spectrum = ifft(C)+K_c; % if K_c is 0 then this is Rayleigh
end

