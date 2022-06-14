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
N = 185;
%t = (-N*Ts:Ts:N*Ts);
t = (0: Ts : (2*N+1)*Ts);
g_t = besselj(1/4,2*pi*fd*abs(t))./sqrt(sqrt(abs(t)));
g_t(1) = sqrt(sqrt(pi*fd))/gamma(5/4);

g_vec = g_t/norm(g_t);
%var_g_vec = var(g_vec)
g_vec_eng = mean(abs(g_vec).^2)   %average energy
%% samples: x
a = 1/sqrt(2); % normalized
%a = mean(g_vec);
x = (randn(Ns,1) +1i*randn(Ns,1))*a; 
%% channel: c
c_1 = conv(x,g_vec);
c_filter = c_1(N+1:Ns+N);  % dicard transients
% c_ray = generate_normalized(c_dis);  % normalized channel : SHOULD BE NORMALIZED???

mean_c_ray = mean(c_filter)   % mean value
eng_c_ray = mean(abs(c_filter).^2) % average energy

%% B. Spectrum method
Gp = generate_filter(fd,fs);        
X = generate_samples(Ns,Gp);

C_spectrum = X'.*Gp;
c_spectrum = ifft(C_spectrum);

mean_c_spectrum = mean(c_spectrum,'all')
var_c_spectrum = var(c_spectrum,0,'all')
%% C. Generate
L = 3;  %taps
M = 64;
tau_c = L*Ts;  %delay spread

for i =1:L
    G_l = generate_filter(fd,fs);
    X_l = generate_samples(N,G_l);
    c_l{i} = generate_rayleigh_channel(X_l,G_l);
end
C_padding = channel_response(c_l,L,M,N);

figure; %%%%%% INCORRECT !!!!!!!
subplot(1,2,1);
mesh(0:M-1, 0:N-1, abs(C_padding'));
subplot(1,2,2);
surf(0:N-1,1:M, abs(C_padding), 'MeshStyle','row')
view(2)


%% D. Simulation 
%% Rician Flat Fading

[PSD_filter,psdf] = pwelch(c_filter,[],[],[],fs,'centered');
PSD_filter = PSD_filter /(sum(PSD_filter)*(fs/Ns));
[PSD_spectrum,psdf2] = pwelch(c_spectrum,[],[],[],fs,'centered');
PSD_spectrum = PSD_spectrum /(sum(PSD_spectrum)*(fs/Ns));
% theoretical PSD

f = (-fd:fs/Ns:fd);

S = 1./(pi*fd*(sqrt(1-(f/fd).^2)));
S(abs(f)>=fd) = 0;

%S=S+fliplr(S);

%PSD_Xaxis= linspace(-f_D,f_D,2000);

figure;
plot(f,pow2db(S))
xlim([-155,155])
hold on
plot(psdf,pow2db(PSD_filter))

plot(psdf2,pow2db(PSD_spectrum))
legend('theory','filter','spectrum');



%% verify kc in 2 channel generating method 
%% and estimate PDF & CDF

kc = [0,1,10];
sigma = 1;  % since the variance of channels are 1
mu = 0;
range = 0:0.01:30;

for i = 1:3
   K_factor(i) = kc(i)^2/(2*sigma^2);
   c_ric_filter = c_filter + kc(i);
   c_ric_spectrum = c_spectrum+kc(i);
   
   ric_dist = fitdist(abs(c_ric_spectrum'),'Rician');
   pdf_rician = pdf(ric_dist,range);
   cdf_rician = cdf(ric_dist,range);
   
   ric_filter_dist = fitdist(abs(c_ric_filter),'Rician');
   pdf_rician_filter = pdf(ric_filter_dist,range);
   cdf_rician_filter = cdf(ric_filter_dist,range);
   
   theory_dist = makedist('Rician','s',kc(i),'sigma',sigma);
   pdf_theory = pdf(theory_dist,range);
   cdf_theory = cdf(theory_dist,range);
   
   % choose either filter or spectrum method to plot pdf and cdf, then compare with
   % theoretical
   figure;subplot(3,1,1);
   plot(range,pdf_rician,range,pdf_theory)
   title(['[Spectrum method] Estimate PDF with kc= ' num2str(kc(i))])
   legend('Ricain','Theoretical')
   
   subplot(3,1,2); % it looks totally the same as spectrum method. 
   plot(range,pdf_rician_filter,range,pdf_theory)
   title(['[Filter method] Estimate PDF with kc= ' num2str(kc(i))])
   legend('Ricain','Theoretical')
   %histgrom(c_ric_filter)   another way to plot PDF
   
   subplot(3,1,3);
   cdfplot(abs(c_ric_spectrum))
   hold on
   plot(range,cdf_theory)
   title(['CDF with kc = ' num2str(kc(i))] )
   legend('Rician','Theoretical')
   
end


%% Estimate  ACF
%delta_t = (0: Ts : (2*N+1)*Ts);

delta_t = linspace(-Ns*Ts,Ns*Ts,Ns);
ACF_theory = besselj(0,2*pi*fd.*delta_t);
ACF_c = xcorr(abs(c_ric_spectrum));
ACF_c = ACF_c./ACF_c(Ns);
ACF_ric = ACF_c(Ns:2*Ns-1);  %% not pretty sure

figure;
plot(delta_t,ACF_theory)
hold on
plot(delta_t,ACF_ric)
legend('theory','Rician')

%% 2. TIME AND FREQUENCY-VARYING RICIAN FADING CHANNELS
Ts_2 = 0.1*10^(-3);
fs_2 = 1/Ts_2;
N_2 = 300;
M_2 = 64;

L_2 = [1 2 3];
fdTs = [0.1 0.005];
kc_2 = [0 1 10];

num = 1;

% C_ric_padding_2 = cell(length(fdTs),length(kc_2),length(L_2));

for i = 1:length(fdTs)
   for j = 1:length(kc_2)
      for k = 1:length(L_2)
          for h = 1:L_2(k)
              G_2 = generate_filter(fdTs(i),fs_2);
              a_G_2 = sqrt(1/(2*var(G_2)));
              X_2 = generate_samples(N_2,a_G_2);
             
              c_ray_2 = generate_rayleigh_channel(X_2,G_2);
              
              % c_ray_2 = c_ray_2/norm(c_ray_2);
              c_ric_mtx_2{h} = c_ray_2 + kc_2(j);
              c_ric_norm_2{h} = c_ric_mtx_2{h}/norm(c_ric_mtx_2{h});

          end
          C_ric_padding_2 = channel_response(c_ric_norm_2,L_2(k),M_2,N_2);   
          figure;
          subplot(1,2,1);
          mesh(0:M_2-1, 0:N_2-1, abs(C_ric_padding_2'));
          title(['fdTs = ' num2str(fdTs(i)) ', kc = ' num2str(kc(j)) ', L = ' num2str(L_2(k))])
          subplot(1,2,2);
          surf(0:N_2-1, 1:M_2, abs(C_ric_padding_2), 'MeshStyle','row')
          view(2)
          title(['fdTs = ' num2str(fdTs(i)) ', kc = ' num2str(kc(j)) ', L = ' num2str(L_2(k))])

      end
      
   end

end

    






