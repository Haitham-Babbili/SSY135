function c_filter = filterMethod(f_D,T_s,N_s,K_c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = N_s/10; % make sure N_s is much larger than N
% initiate g vector and store g(nTs) values in it
t = (-N*T_s:T_s:N*T_s);
g = besselj(1/4,2*pi*f_D*abs(t))/nthroot(abs(t),4);
g(t==0) = nthroot(pi*f_D,4)/gamma(5/4);

% make sure g is of unit energy
g=g/norm(g);

% generate samples from zero-mean unit-variance complex Gaussian dist, we
% must divide by sqrt(2) since adding a N(0,1) with a N(0,1) gives a N(0,2)
x = (randn(N_s,1) + 1i*randn(N_s,1))/sqrt(2);

% convolve x with g and remove the transients (can be done using 'same')
c_filter = conv(x,g,'same') + K_c; % if K_c is 0 then this is Rayleigh

end

