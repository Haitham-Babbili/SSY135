Ts = 0.1e-3;                    %symbol time Ts in s
Ns = 1e6;                       %number of samples per channel
N = 500;                        %filterlength
fD = 55;                    	%doppler shift in Hz

g=zeros(1,(2*N+1));             %create filter 
for n=-N:N
    if n == 0
        g(1,n+N+1) = nthroot(pi*fD,4)/gamma(5/4);
    else
        g(1,n+N+1) = besselj(1/4,2*pi*fD*abs(n*Ts)/nthroot(abs(n*Ts),4));
    end
end
g = g/sqrt(sum(g.*g));          % g should be unit energy

x = zeros(1,Ns);                %create samples
for m = 0:(Ns-1)
    x(1,m+1) = randn+1i*randn;
end
x = (x-mean(x))/std(x);

figure(1);
hold on
plot(real(x), imag(x), '.'), title('name');

c = conv(x,g,'same');
figure(2);
subplot(1,2,1)
histogram(real(c),100), title('real part');
subplot(1,2,2)
histogram(imag(c),100), title('imag part');
figure(3);
histogram(abs(x),150) ,title('Rayleigh pdf');

variance = var(c)              %should be 1
expectation = mean(c)          %should be 0