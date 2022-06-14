%% Part 1
clc, clear all, close all

fc=2e9; %in Hz
Ts = 0.1e-3;
Fs = 1/Ts;
v = 30; %in km/h
Ns = 1e4;
N = 300;
fD = (v/(3.6*physconst('LightSpeed')))*fc;
max_Ts = 1/(2*fD);
if max_Ts < Ts
    error('Ts too large, aliasing occurs')
end

%% Rayleigh channel
c_filter = channelByFilter(Ts,Ns,N,fD);
c_spectrum = channelBySpectrum(Ts,Ns,fD);
%% psd

[p_filter,f] = pwelch(c_filter,[],[],[],Fs,'centered');
[p_spectrum,f_] = pwelch(c_spectrum,[],[],[],Fs,'centered');


%% Theoretical
sc = zeros(1,Fs);
for i = -Fs/2:0
    sc(Fs/2+i+1) = Gp(i+Fs,fD,Fs)^2;
end
%% Compare PSD
figure
plot(f,pow2db(p_filter)), hold on
plot(f_,pow2db(p_spectrum)),hold on
plot(linspace(-Fs/2,Fs/2-1,length(sc)),pow2db(sc))
hold off
xlim([-200 200])
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
legend('filter','spectrum','theoretical')

%% pdfs and cdfs

% -- For filter method ---
% c + kc for rician
c_rician_0_filter = c_filter; % kc = 0;
c_rician_1_filter = c_filter +1*ones(1,Ns); % + kc = 1
c_rician_10_filter = c_filter +10*ones(1,Ns); % + kc = 10

% Get pdf
[f_r_0,xi_r_0] = ksdensity(abs(c_rician_0_filter),'Support','positive','Function','pdf');
[f_r_1,xi_r_1] = ksdensity(abs(c_rician_1_filter),'Support','positive','Function','pdf');
[f_r_10,xi_r_10] = ksdensity(abs(c_rician_10_filter),'Support','positive','Function','pdf');

% -- For spectrum method ---
% c + kc for rician
c_rician_0_spectrum = c_spectrum; % kc = 0;
c_rician_1_spectrum = c_spectrum +1*ones(1,Ns); % + kc = 1
c_rician_10_spectrum = c_spectrum +10*ones(1,Ns); % + kc = 10

% Get pdf
[f_r_0_s,xi_r_0_s] = ksdensity(abs(c_rician_0_spectrum),'Support','positive','Function','pdf');
[f_r_1_s,xi_r_1_s] = ksdensity(abs(c_rician_1_spectrum),'Support','positive','Function','pdf');
[f_r_10_s,xi_r_10_s] = ksdensity(abs(c_rician_10_spectrum),'Support','positive','Function','pdf');

%Create rician theoretical distributions
theory_rician_0 = makedist('Rician','s',0);
theory_rician_1 = makedist('Rician','s',1);
theory_rician_10 = makedist('Rician','s',10);


figure
sgtitle('pdf')
subplot(3,1,1)
    plot(xi_r_0,f_r_0,'r','LineStyle','--'), hold on, title('Kc = 0') % filter
    plot(xi_r_0_s,f_r_0_s,'b','LineStyle',':') %spectrum
    plot(linspace(0,max(xi_r_0)),theory_rician_0.pdf(linspace(0,max(xi_r_0))))% theory
    legend('filter','spectrum','theoretical')
    
subplot(3,1,2)
    plot(xi_r_1,f_r_1,'r','LineStyle','--'), hold on, title('Kc = 1') % filter
    plot(xi_r_1_s,f_r_1_s,'b','LineStyle',':') %spectrum
    plot(linspace(0,max(xi_r_1)),theory_rician_1.pdf(linspace(0,max(xi_r_1)))) % theory

subplot(3,1,3)
    plot(xi_r_10,f_r_10,'r','LineStyle','--'), hold on, title('Kc = 10') % filter
    plot(xi_r_10_s,f_r_10_s,'b','LineStyle',':') %spectrum
    plot(linspace(0,max(xi_r_10)),theory_rician_10.pdf(linspace(0,max(xi_r_10))))% theory

% Get cdf
[f_r_0,xi_r_0] = ksdensity(abs(c_rician_0_filter),'Support','positive','Function','cdf');
[f_r_1,xi_r_1] = ksdensity(abs(c_rician_1_filter),'Support','positive','Function','cdf');
[f_r_10,xi_r_10] = ksdensity(abs(c_rician_10_filter),'Support','positive','Function','cdf');

% Get cdf
[f_r_0_s,xi_r_0_s] = ksdensity(abs(c_rician_0_spectrum),'Support','positive','Function','cdf');
[f_r_1_s,xi_r_1_s] = ksdensity(abs(c_rician_1_spectrum),'Support','positive','Function','cdf');
[f_r_10_s,xi_r_10_s] = ksdensity(abs(c_rician_10_spectrum),'Support','positive','Function','cdf');

figure
sgtitle('cdf')
subplot(3,1,1)
    plot(xi_r_0,f_r_0,'r','LineStyle','--'), hold on, title('Kc = 0') % filter
    plot(xi_r_0_s,f_r_0_s,'b','LineStyle',':') %spectrum
    plot(linspace(0,max(xi_r_0)),theory_rician_0.cdf(linspace(0,max(xi_r_0))))% theory
    legend('filter','spectrum','theoretical')
    
subplot(3,1,2)
    plot(xi_r_1,f_r_1,'r','LineStyle','--'), hold on, title('Kc = 1') % filter
    plot(xi_r_1_s,f_r_1_s,'b','LineStyle',':') %spectrum
    plot(linspace(0,max(xi_r_1)),theory_rician_1.cdf(linspace(0,max(xi_r_1))))% theory

subplot(3,1,3)
    plot(xi_r_10,f_r_10,'r','LineStyle','--'), hold on, title('Kc = 10') % filter
    plot(xi_r_10_s,f_r_10_s,'b','LineStyle',':') %spectrum
    plot(linspace(0,max(xi_r_10)),theory_rician_10.cdf(linspace(0,max(xi_r_10))))% theory
%% Autocorrelation

% Do autocorrelation for each kc and for both methods
corr_filter_0 = xcorr(abs(c_rician_0_filter),N,'unbiased');
corr_spectrum_0 = xcorr(abs(c_rician_0_spectrum),N,'unbiased');
corr_filter_1= xcorr(abs(c_rician_1_filter),N,'unbiased');
corr_spectrum_1 = xcorr(abs(c_rician_1_spectrum),N,'unbiased');
corr_filter_10 = xcorr(abs(c_rician_10_filter),N,'unbiased');
corr_spectrum_10 = xcorr(abs(c_rician_10_spectrum),N,'unbiased');

Ac = @(N,kc) besselj(0,2*pi*fD*N) + kc^2; % Theoretical Ac for Rician channel (C_Rician = C_Rayleigh + kc, C_Rician complex N~(0,1)) => ACF_rician(dt) = ACF_[Rayleigh](dt)+ kc^2)
x = (-N:N)*Ts;
figure
subplot(3,1,1)
    plot(x,corr_filter_0), hold on , xlabel('lags'), ylabel('autocorr'), title('kc = 0')
    plot(x,corr_spectrum_0)
    plot(x,Ac(x,0));
subplot(3,1,2)
    plot(x,corr_filter_1), hold on, xlabel('lags'), ylabel('autocorr'), title('kc = 1')
    plot(x,corr_spectrum_1)
    plot(x,Ac(x,1));
subplot(3,1,3)
    plot(x,corr_filter_10), hold on, xlabel('lags'), ylabel('autocorr'), title('kc = 10')
    plot(x,corr_spectrum_10)
    plot(x,Ac(x,10));
