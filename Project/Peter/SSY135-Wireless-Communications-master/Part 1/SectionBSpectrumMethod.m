clear, clc, close all

Ts = 1e-2;                                  %sampling time Ts in [s]
fs = 1/Ts;                                  %samplingrate in [Hz]
Ns = 1e6;                                   %number of samples
fD = 40;                                    %doppler shift in [Hz]

G=zeros(1,Ns);                              %initialising periodic extension with zeros
for k=0:Ns-1
      G(1,k+1) = Gp(k*fs/Ns,fD,fs);         %calculate values of G
end

sumg = 0;                                   %summation of G for calculating a
for k = 0:(Ns-1)
    sumg = sumg + Gp(k*fs/Ns,fD,fs)^2;      
end
a = Ns*sqrt(1/(2*sumg));                    %calculate this value by parsevall equatiion 

                                            %draws complex random numbers
X = a*(randn(1,Ns)+1i*randn(1,Ns)); 

C = X.*G;                                    %multiplicate X with G elementwise
c = ifft(C);                              %apply an inverse DFT to C

var_c = var(c)                              %Check that E[|c(nTs)|^2] = 1

%plots
figure(1)
histogram(abs(c),'Normalization','probability');
ylabel('probability density');
print -depsc c_spectrum_rayleigh;

figure(2);
y = 0:(fs/Ns):fs-fs/Ns;                     %periodic extension of G for intervall [0:fs]
subplot(2,2,1);
plot(y,G), title('periodic extension of G');
subplot(2,2,2);                             %should look like Rayleigh distribution
histogram(abs(c),40), title('Rayleigh pdf of c');
subplot(2,2,3);                             %should look like Gaussian distribution
histogram(real(c),50), title('pdf of real part of c');
subplot(2,2,4);                             %should look like Gaussian distribution
histogram(imag(c),50), title('pdf of imag part of c');
print -deps epsFig;
