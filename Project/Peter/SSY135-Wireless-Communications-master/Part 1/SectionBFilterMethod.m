clear, clc, close all

Ts = 0.1e-3;                    %symbol time Ts in s
Ns = 1e6;                       %number of samples per channel
N = 500;                        %filterlength
fD = 55;                    	%doppler shift in Hz

                                    %create filter 
n=(-N:N)*Ts;
g = besselj(1/4,2*pi*fD*abs(n))./nthroot(abs(n),4);
g(N+1) = nthroot((pi*fD),4)/gamma(5/4);

g = g./sqrt(sum(g.^2));             %g should be unit energy

x = (randn(1,Ns)+1i*randn(1,Ns));
x = (x-mean(x))./std(x);                      %zero mean and unit variance
c = conv(x,g,'same');                   %convolve x with g, parameter 'same' removes the first and last N samples 
c = c./ sqrt(var(c));

figure(1);
plot(real(x), imag(x), '.'), title('title');

figure(2);
histogram(real(c),'Normalization','probability');
ylabel('probability density');
print -depsc pdf_filter_real_gauss;

figure(3);
histogram(imag(c),'Normalization','probability');
ylabel('probability density');
print -depsc pdf_filter_imag_gauss;

figure(4);
histogram(abs(x),'Normalization','probability');
ylabel('probability density');
print -depsc pdf_filter_rayleigh;

variance = var(c)              %should be 1
expectation = mean(c)          %should be 0