function c = channelBySpectrum(Ts, Ns, fD)

%Ts = 1e-2;                                  %sampling time Ts in [s]
fs = 1/Ts;                                  %samplingrate in [Hz]
%Ns = 5e3;                                   %number of samples
%fD = 40;                                    %doppler shift in [Hz]

G=zeros(1,Ns);                              %initialising periodic extension with zeros
for k=0:Ns-1
      G(1,k+1) = sqrt(Gp(k*fs/Ns,fD,fs));   %calculate values of G
end

sumg = 0;                                   %summation of G for calculating a
for k = 0:(Ns-1)
    sumg = sumg + Gp(k*fs/Ns,fD,fs);        %Gp(k*fs/Ns,fD,fs); %%somethin is here wrong
end
a = sqrt(Ns/sumg);                          %calculate this value by parsevall equatiion 

X = zeros(Ns,Ns);                           %initialising samples with zeros
for k=0:Ns-1                                %draws complex random numbers
X(:,k+1) = (a/sqrt(2))*(randn(Ns,1)+1i*randn(Ns,1)); 
end

C = G*X;                                    %multiplicate X with G
c = ifft(C);                                %apply an inverse DFT to C
end