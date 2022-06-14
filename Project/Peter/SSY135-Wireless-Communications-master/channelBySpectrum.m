function c = channelBySpectrum(Ts, Ns, fD)

%Ts = 1e-2;                                     %sampling time Ts in [s]
fs = 1/Ts;                                      %samplingrate in [Hz]
%Ns = 5e3;                                      %number of samples
%fD = 40;                                       %doppler shift in [Hz]

G=zeros(1,Ns);                                  %initialising periodic extension with zeros
for k=0:Ns-1
      G(1,k+1) = Gp(k*fs/Ns,fD,fs);             %calculate values of G
end

sumg = 0;                                       %summation variable for calculating a
for k = 0:(Ns-1)
    sumg = sumg + Gp(k*fs/Ns,fD,fs)^2;          
end
a = sqrt(1/(sumg*2))*Ns;                          %calculate this value by parsevall equatiion 

                                                %draws complex random numbers
X= a*(randn(1,Ns)+1i*randn(1,Ns)); 

C = X.*G;                                        %multiplicate X with G
cn = ifft(C);                                    %apply an inverse DFT to C
c = cn/sqrt(var(cn));
end