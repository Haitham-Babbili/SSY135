function c = channelByFilter(Ts,Ns,N,fD)
    %Ts: symbol time Ts in s
    %Ns: number of samples per channel
    %N: filterlength
    %fD: doppler shift in Hz

g=zeros(1,(2*N+1));             %create filter 
for n=-N:N
    if n == 0
        g(1,n+N+1) = nthroot(pi*fD,4)/gamma(5/4);
    else
        g(1,n+N+1) = besselj(1/4,2*pi*fD*abs(n*Ts)/nthroot(abs(n*Ts),4));
    end
end
g = g/sqrt(sum(g.*g));          %g should be unit energy

x = zeros(1,Ns);                %create samples
for m = 0:(Ns-1)
    x(1,m+1) = randn+1i*randn;
end
x = (x-mean(x))/std(x);         %zero mean and unit variance

c = conv(x,g,'same');           %convolve x with g 

end

