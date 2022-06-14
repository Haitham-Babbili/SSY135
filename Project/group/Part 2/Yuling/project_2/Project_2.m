clc
clear all
close all
% Project II

% given parameters
fc = 2e9;
B = 1e6; 
Ts = 1/B;
v = 15;    % m/s
c = 3e8;
tau = 5.4*10^(-6);
Ptx = 0.1;  % w
fd = v*fc/c;  % doppler frequency
M = 4; 

% chosen parameters
N = 2^6;
Ncp = 2^4; % N > Ncp, Ncp >= L-1
N_sym = round(N/Ncp); % number of symbols
N0 = 2*2.07*10^(-20);  % 2*2.07*10^(-20) W/Hz % becuase the BW is duple side so we don't multiply by 2
L = 6; 

a = 1/sqrt(2);

pl = 101; %dB
Ptx_dBm = 10*log10(Ptx/10^(-3)); % W -> dBm
Prx_dBm = Ptx_dBm - pl;
Prx = 10^(-3) * 10^(Prx_dBm/10); % dBm -> W  ??????
E = Prx * (N + Ncp) * Ts;
EbN0 = E*log2(M)/N0;
sigma = sqrt(1/EbN0);

%% B. Simulation Task
%% QPSK test

const = [1+1i, 1-1i, -1-1i, -1+1i] * a;
tau_test = [0 5 10 15 20 25];
P_test = [Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L];

figure;
plot(const,'ro')
hold on
for j = 1:L
    s_test_tx = randi([0 1], 1, N*log2(M)); % bits
    
    z_test_tx = transmitter(s_test_tx, N, Ncp);
    
    [r_test_cp, ~] = Fading_Channel(z_test_tx, tau_test, fd*Ts, P_test);
    
    len_test = length(r_test_cp);
    if i == 1
        a1 = ones(len_test,2);
        h_test = a1 + a1*1i;
        
    else 
        a1 = zeros(len_test,2);
        h_test = a1 + a1*1i;
    end
    
    noise_test = sigma * (randn(length(r_test_cp),1) + 1i * randn(length(r_test_cp),1));
    r_test_cp = r_test_cp + noise_test;
    [s_test_est_sym, s_test_est_bits] = receiver(r_test_cp, h_test, const, N, Ncp, j);
    
    errors_test(j) = sum(s_test_est_bits ~= s_test_tx);
    
    plot(s_test_est_sym, 'b.')
    title('Simulation test')
    hold on
end
BER_test = sum(errors_test)/(N_sym * N * 2); %total bits = N_sym * N * 2
disp(['BER_test = ' num2str(BER_test)])

%% with time-varying frequency selective channel
M = 4;   % QPSK
const = [1+1i, 1-1i, -1-1i, -1+1i] * a;
tau_samples = [0 5 10 15 20 25];
PDP = [Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L];
errors = zeros(N_sym, 1);

figure;
plot(const, 'ro')
hold on
for i = 1: N_sym
    s_tx = randi([0 1], 1, N*log2(M)); % bits
    
    z_tx = transmitter(s_tx, N, Ncp); % symbols with cp
    % why the length of z_tx is not the same as r_cp???
    [r_cp, h] = Fading_Channel(z_tx, tau_samples, fd*Ts, PDP); % r_cp: with cp
    
    % generate noise
    noise = randn(length(r_cp),1) + 1i * randn(length(r_cp),1);
    % Wn = fft(noise);
    r_cp1 = r_cp + noise; 
    [s_est_sym, s_est_bits] = receiver(r_cp1, h, const, N, Ncp, i); % symbols
    
    errors(i) = sum(s_est_bits ~= s_tx); 
    
    plot(s_est_sym, 'b.')
    title('Time-varying frequency selective channel')
    hold on
end

BER = sum(errors)/(N_sym * N * 2); %total bits = N_sym * N * 2
disp(['BER = ' num2str(BER)])



