function c = channelByFilter(Ts,Ns,N,fD)
    %Ts: symbol time Ts in s
    %Ns: number of samples per channel
    %N: filterlength
    %fD: doppler shift in Hz

                                    %create filter 
n=(-N:N)*Ts;
g = besselj(1/4,2*pi*fD*abs(n))./nthroot(abs(n),4);
g(N+1) = nthroot((pi*fD),4)/gamma(5/4);

g = g./sqrt(sum(g.^2));             %g should be unit energy

x = (randn(1,Ns)+1i*randn(1,Ns));
x = (x-mean(x))./std(x);                      %zero mean and unit variance
c = conv(x,g,'same');                   %convolve x with g, parameter 'same' removes the first and last N samples 
c = c./ sqrt(var(c));

end

