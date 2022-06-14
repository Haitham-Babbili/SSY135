clear all;
clc;
close all;

 
fc= 2e9;        % carrier frequency
v=8.33;         % speed of transition
lamda= 3e8/fc;  % wave length
fd= v/lamda;    % dopplar effect
Ts=0.1e-3;      % sample time
fs= 1/Ts;       % sample frequency
Ns=100000;       % number of sample for the channel in frequency
N= 1000;         % number of sample for the channel in time 
kc=1;           % power of line of site components
L=1;            % number of tabs
M=64;

%% b part
n=-N:N;
l=1:L;

g = g_filter(N,Ts,fd); % generat g filter in time domain
x = (randn(1, Ns) + 1i * randn(1, Ns)) * sqrt(1/2);  % generat x(t)
% x = (randn(1, Ns) + 1i * randn(1, Ns)) * sqrt(1/a);  % generat x(t)
c_t = conv(x,g, 'same');  % the channel by convolve g and x(t)
E_ct = mean(abs(c_t).^2,2); % energy of c E(c^2) = 1
meu = mean(real(c_t),2); % expected values of c E(c) = 0

figure(1);
hold on
histogram(abs(x)) 
title('x distribution (Rayleigh of g)');

 
G = spectrum_filter(Ns,Ts,fd); % generat G in frequency domain
a= Ns^2/(norm(G)^2)/2;
X=(randn(1, Ns) + 1i * randn(1, Ns)) * sqrt(1/a);  % generat x(t)
C=G.*X;  % the channel by convolve G and x in frequency it's multiplication
c= ifft(C,Ns,2);  % IDFT for the channel

E_c = mean(abs(c).^2,2); % energy of c E(c^2) = 1
mean_c_spectrum = mean(c,'all')
var_c_spectrum = var(c,0,'all')

%% c part


for l=1:L
    filter=['g', 'G'];
    
    switch filter
        case 1
             c_t(l,:)= conv(x,g,'same');
     
     channel=['Re','Ric'];
     switch channel
         case 1
             % Rayleigh in time 
             c_energy= abs(c_t(l,:)).^2; % unit energy
             c(l,:) = c_t(l,:)./sqrt(mean(c_energy)); % normalize the channel
         case 2
             % Rician
             c_Ri_t(l,:)= c_t(l,:) + kc;
             c_Ri_t_energy= abs(c_Ri(l,:)).^2; % unit energy
             c(l,:) = c_Ri(l,:)./sqrt(mean(c_Ri_energy)); % normalize the channel
     end
     

    case 2
       C(l,:)= G.*x;
       c=ifft(C,Ns,2);
       channel=['Re','Ric'];
       switch channel
         case 1
             C_energy= abs(c(l,:)).^2; % unit energy
             c(l,:) = c(l,:)./sqrt(mean(C_energy)); % normalize the channel
         case 2
             C_Ri(l,:)= G.*x + kc;
             c_Ri_energy= abs(C_Ri(l,:)).^2; % unit energy
             c(l,:) = C_Ri(l,:)./sqrt(mean(C_Ri_energy));% normalize the channel
       end
    end
    
end




% % Theoretical Doppler PSD
% fd_vec = linspace(-fd,fd,100); % Doppler frequency as a vector
% Sc = 1./(pi*fd.*sqrt(1-(fd_vec./fd).^2)); % general low for psd 
% 
% % Estimated PSD
% % PSD_est_c = pwelch(c_t); % psd for g in time domain
% PSD_est_C = pwelch(c_t(1,:), [], [], [], fs); % psd for time domain channel in frequency



% Theoretical Doppler PSD
fd_vec = linspace(-fd,fd,100); % Doppler frequency as a vector
Sc = 1./(pi*fd.*sqrt(1-(fd_vec./fd).^2)); % general low for psd 
Sc_dB= 10*log10(Sc);

% Estimated PSD
% PSD_est_c = pwelch(c_t); % psd for g in time domain
[PSD_est_C,w] = pwelch(c(1,:), [], fs); % psd for time domain channel in frequency
PSD_est_dB= 10.*log10(PSD_est_C);
% 
% [PSD_spectrum,w] = pwelch(c,1000,[],Ns);
% PSD_spectrum = PSD_spectrum /(sum(PSD_spectrum)*(fs/Ns));
figure(8); 
hold on;
plot(PSD_est_C, PSD_est_dB)
plot(linspace(-60, fd, length(PSD_est_dB)),PSD_est_dB)
% plot(w/pi,PSD_est_dB)
legend('Spectrum method', 'Filter method', 'Theory')
title('Power Spectral Density Estimate')
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
figure(9);
plot(linspace(-60, fd, length(Sc_dB)),Sc_dB, 'k', 'LineWidth', 3)
legend('Spectrum method', 'Filter method', 'Theory')
title('Theoretical Power Spectral Density')
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');


% figure(5)
% mesh(cc,abs(Cf))
% 
% %% d part
% 
% % Rayleigh in time 
% c_energy= abs(c(l,:)).^2; % unit energy
% c(l,:) = c(l,:)./sqrt(mean(c_energy)); % normalize the channel
% 
% % Rayleigh in frequency
% C_energy= abs(C(l,:)).^2; % unit energy
% C(l,:) = C(l,:)./sqrt(mean(C_energy)); % normalize the channel
% 
% 
% % % Rician in time
% % c_Ri(l,:) = conv(x(l,:),g,'same');
% % c_Ri(l,:)= c_Ri(l,:) + kc;
% % c_Ri_energy= abs(c_Ri(l,:)).^2; % unit energy
% % c_Ri(l,:) = c_Ri(l,:)./sqrt(mean(c_Ri_energy)); % normalize the channel
% 
% % Rician in frequency
% C_Ri(l,:)= G.*x + kc;
% C_Ri_energy= abs(C_Ri(l,:)).^2; % unit energy
% C_Ri(l,:) = C_Ri(l,:)./sqrt(mean(C_Ri_energy));% normalize the channel


ACF = besselj(0,2*pi*fd.*n.*Ts);
% ACF_estmat = autocorr(c(1,:).',10);
ACF_estmat_filter = real(xcorr(c_t(1,:),N,'unbiased'));
ACF_estmat_spectrom = real(xcorr(c(1,:).',N,'unbiased'));

figure(5) 
plot(n.*Ts,ACF)
title('Theoretical ACF')
legend('filter method','spectrum method','theory')
xlabel('t/s')
ylabel('autocorrelation')

figure; hold on;
plot(n.*Ts,ACF_estmat_spectrom)
title('est ACF')


c = [c; zeros(M-L,Ns)]; % Zero padding
Cf = fft(c,M,1);



