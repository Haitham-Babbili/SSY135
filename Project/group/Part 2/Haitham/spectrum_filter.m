function Gp = spectrum_filter(Ns,Ts,fd)

% Gp=zeros(1,Ns);
%            
% k=0:Ns-1;
f= 0:fd;

if abs(f) <= fd
    sc = 1./(pi*fd.*sqrt(1-(f./fd).^2));
else
    sc = 0;
end

G = sqrt(sc);

Gp = [G(length(G)/2+1:end) zeros(1,Ns-length(G)) G(1:length(G)/2)]; 


sc = abs(G).^2; 

figure(2)
plot(f,real(G))
xlabel('Frequency [Hz]')
ylabel('G(f)')
title('PSD of G(f)')
end

