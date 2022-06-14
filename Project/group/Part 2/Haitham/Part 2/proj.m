%% Proj part 2
clc
clear all

%=== parameters ====%

f_c = 2e9; % 2GHz frequency carrier
v = 15; % 30km/h velocity
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
B= 1e5; % bandwidth
f_s = 1e6; % 0.1 ms sample interval
T_s = 1/f_s; % sampling frequency
fdTs=f_D*T_s;
td= 4.5e-6; %delay 
N0 = 2.07e-20*B; % multiply by BW du to unit energy in [W]
N = 2^6; % sequence length
L = ceil(td/T_s); % number of channel taps
Ncp = L-1; % cyclic prefix
tau = (0:1:L-1); % Path delays in samples
M = 4;
% P = (1/length(tau))*ones(1,length(tau)); % Flat power delay profile
Pt= 0.1; % power Tx in W
P = [Pt/L Pt/L Pt/L Pt/L Pt/L];
Pt_dB=10*log10(Pt/10^(-3));
PL_dB=101; %path loss indB
PL= 10.^(PL_dB/10); % path loss in W
Prx_dBm = Pt_dB - PL_dB;
Prx = 10^(-3) * 10^(Prx_dBm/10); % dBm -> W 
E = Prx * (N + Ncp) * T_s;
EbN0 = E*log2(M)/N0;
sigma = sqrt(1/EbN0);


samp_num = 40; % Number of time samples
if (N + Ncp)*fdTs*10 > 1 || N*T_s > 1/(2*f_D) % Check constraint and Check if N*Ts is smaller than coherence time
    error('N or Ncp to large')
end

To = (N+Ncp)*T_s;
n = 1:N+Ncp;


%% ==== TX ==== %

symbs = zeros(N,samp_num);
bits_s = zeros(2*N,samp_num);
y = zeros(N,samp_num);
for i = 1:samp_num
    
[z]= qpsk_seq_gen(N,E,T_s,Ncp);
    
    % == Channel  == %
    [z, h] = Fading_Channel(z, tau, fdTs, P);
%     z = z(1:end-L+1); % Remove delayed symbols
    
    %c_output = awgn(z,i,'Measured');
%     z = z*sqrt(10^(-10.1)); %pathloss add 
    z = z*sqrt(PL); % Pathloss
    c_output = z + (randn(length(z),1) +1i*rand(length(z),1))*sqrt(N0);   % add noise
   
    % == Reciever  == %
    Rx = c_output(Ncp+1:end); %remove cp
    
    ch = fft(h(1,:),N); % FFT of channel
    
    y(:,i) = fft(Rx,N)/sqrt(N); % FFT at received symbols
    y(:,i) =y(:,i)./ch.';  % zero forcing 
    
end

% scatterplot 


scatterplot(y(1,:)); 
scatterplot(y(25,:)); 
scatterplot(y(50,:));
% scatterplot(y(75,:)); 
% scatterplot(y(100,:)); 

%%




