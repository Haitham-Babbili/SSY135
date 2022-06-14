function c = channelBySpectrum(Ts, Ns, fD)
G=zeros(1,Ns);             %create filter 
for k=0:Ns-1
    for m=-100:100
    G(1,k+1) = G(1,k+1) + sqrt(Sc((k/Ns-m)/Ts,fD));
    end
end

X = zeros(Ns,Ns);                %create samples
for m = 0:(Ns-1)
    X(:,m+1) = a*(randn(Ns,1)+1i*randn(Ns,1));
end
%X = X-mean(X);                  % zero mean
%X = (X-mean(X))/std(X);        %zero mean and unit variance
C = G*X;
size(C)
c = ifft(C);
end

