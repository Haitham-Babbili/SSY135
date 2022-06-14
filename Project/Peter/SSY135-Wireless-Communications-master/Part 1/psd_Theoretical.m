% PSD of Gf
V=30;
v=(1000/3600) .*V;
fc=2e9;
c=3e8;
Ts = 0.1e-3;
fs=1/Ts;
Ns=10000;
lamda=fc/c; % wavelength in metre
fD=(v*fc)./c ; %fD is maximum when cos0 =1
% fD= 55.5556 Hertz

% To avoid aliasing (fs=1/Ts) > 2*fD
f = -fD:fs/Ns:fD;
Sf = 1./(pi*fD.*sqrt(1-(f./fD).^2));
        Sf(1) = Sf(end); % This is to Get rid of infinite spread
Gf = sqrt(Sf);
        T_psd = abs(Gf).^2; 
        figure(1),clf, plot(f,real(Gf))
        xlabel('Frequency [Hz]')
        ylabel('G(f)')
        title('Theoretical PSD of G(f)')
        
        figure(2),clf, plot(f,Sf)
        xlabel('Frequency [Hz]')
        ylabel('S(f)')
        title('Theoretical PSD of S(f)')
      