clc
clear all 
close all

%% 1. SIMULATION TASK

fc = 2*10^9;
v = 8.3; % 30km/h -> m/s
Ts = 0.1*10^(-3);
c = 3*10^8; % light speed [m/s]

%% A. Filter method
fs = 1/Ts; % sample rate
fd = v*fc/c; % doppler freq 
Ns = 10^4;
%% filter: g vector
N = 1000;
%t = (-N*Ts:Ts:N*Ts);
t = (0: Ts : (2*N+1)*Ts);
g_t = besselj(1/4,2*pi*fd*abs(t))./sqrt(sqrt(abs(t)));
g_t(1) = sqrt(sqrt(pi*fd))/gamma(5/4);
g_vec = g_t./sqrt(mean(abs(g_t).^2));
g_vec_eng = mean(abs(g_vec).^2)   %average energy
%% samples: x
a = 1/sqrt(2); % normalized
%a = mean(g_vec);
x = (randn(Ns,1) +1i*randn(Ns,1))*a; 
%% channel: c
c_1 = conv(x,g_vec);
c_dis = c_1(N+1:Ns+N);  % dicard transients
c_ray = generate_normalized(c_dis);  % normalized channel : SHOULD BE NORMALIZED???

mean_c_ray = mean(c_ray)   % mean value
eng_c_ray = mean(abs(c_ray).^2) % average energy

%% B. Spectrum method
Gp = generate_filter(fd,fs);        
%Gp_norm = Gp/norm(Gp);
%a_try = (Ns/norm(Gp))^2;
%a_try = sqrt(1/(2*var(Gp_norm)));
a_try = Ns^2/(norm(Gp)^2)/2;
X = generate_samples(Ns,a_try);

%c_spectrum = generate_rayleigh_channel(X,Gp);
C_spectrum = X'.*Gp;
c_spectrum = ifft(C_spectrum);
%eng_c_spectrum = mean(abs(c_spectrum).^2)
mean_c_spectrum = mean(c_spectrum,'all')
var_c_spectrum = var(c_spectrum,0,'all')
%% C. Generate
Ts_c = 0.1*10^(-3);
fs_c = 1/Ts_c;
L = 3;  %taps
M = 64;
tau_c = L*Ts_c;  %delay spread

for i =1:L
    G_l = generate_filter(fd,fs);
    a_G_l = sqrt(1/(2*var(G_l)));
    X_l = generate_samples(N,a_G_l);
    c_l{i} = generate_rayleigh_channel(X_l,G_l);
end
C_padding = channel_response(c_l,L,M,N);

figure;
mesh(0:M-1, 0:N-1, abs(C_padding'));
figure;
surf(0:N-1,1:M, abs(C_padding), 'MeshStyle','row')
view(2)


%% D. Simulation 
%% Rician Flat Fading
kc_ric = 1;  % kc >= 0
c_rff= c_ray + kc_ric; 
c_ric = generate_normalized(c_rff);

mean_c_ric = mean(c_ric)   % mean value
eng_c_ric = mean(abs(c_ric).^2) % average energy

f_ric = (-fd:fd);
Sc_theory = (1/(pi*fd)).*(1./sqrt(1-(f_ric./fd).^2)); 

[S_estimate,w] = pwelch(c_ric);
figure;
plot(w/pi,10*log10(S_estimate))
title('estimated PSD')
xlabel('w/pi')
ylabel('sepctrum [dB]')
figure;
plot(f_ric,10*log10(Sc_theory))
title('theoretical PSD')
xlabel('frequency [Hz]')
ylabel('spectrum [dB]')

%% verify kc in 2 channel generating method 
%% and estimate PDF & CDF

kc = [0,1,10];
sigma = 1;  % since the variance of channels are 1
mu = 0;

for i = 1:3
   K_factor(i) = kc(i)^2/(2*sigma^2);
   c_ric_filter{i} = generate_normalized(c_ray+kc(i));
   c_ray_spectrum = generate_rayleigh_channel(X,Gp);
   c_ric_spectrum{i} = generate_normalized(c_ray_spectrum+kc(i));
   figure;
   histogram(abs(c_ric_spectrum{i}));
   title(['Estimated PDF with kc= ' num2str(kc(i))])
   
   cdfx = cumsum(c_ric_spectrum{i})/sum(c_ric_spectrum{i});
   figure;
   plot(cdfx)
   title(['Estimated CDF with kc= ' num2str(kc(i))])
   
   %%%% MISS: THEORETICAL PDF & CDF and compare to estimated one %%%%%%
end


%% Estimate  ACF
delta_t = (0: Ts : (2*N+1)*Ts);
Ac_theory = besselj(0,2*pi*fd*delta_t);

%% 2. TIME AND FREQUENCY-VARYING RICIAN FADING CHANNELS
Ts_2 = 0.1*10^(-3);
fs_2 = 1/Ts_2;
N_2 = 300;
M_2 = 64;

L_2 = [1 2 3];
fdTs = [0.1 0.005];
kc_2 = [0 1 10];

% num = 1;

C_ric_padding_2 = cell(length(fdTs),length(kc_2),length(L_2));

for i = 1:length(fdTs)
   for j = 1:length(kc_2)
      for k = 1:length(L_2)
          for h = 1:L_2(k)
              G_2 = generate_filter(fdTs(i),fs_2);
              a_G_2 = sqrt(1/(2*var(G_2)));
              X_2 = generate_samples(N_2,a_G_2);
             
              c_ray_2 = generate_rayleigh_channel(X_2,G_2);
              c_ric_mtx_2{h} = c_ray_2 + kc_2(j);

          end
          C_ric_padding_2{i,j,k} = channel_response(c_ric_mtx_2,L_2(k),M_2,N_2);   
          figure;
          subplot(1,2,1);
          mesh(0:M_2-1, 0:N_2-1, abs(C_ric_padding_2{i,j,k}'));
          title(['fdTs = ' num2str(fdTs(i)) ', kc = ' num2str(kc(j)) ', L = ' num2str(L_2(k))])
          subplot(1,2,2);
          surf(0:N_2-1, 1:M_2, abs(C_ric_padding_2{i,j,k}), 'MeshStyle','row')
          view(2)
          title(['fdTs = ' num2str(fdTs(i)) ', kc = ' num2str(kc(j)) ', L = ' num2str(L_2(k))])

      end
      
   end
end

    






