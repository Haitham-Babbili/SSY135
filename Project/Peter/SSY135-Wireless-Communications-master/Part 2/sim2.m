%% Project 2
clc, clear all, close all
% == Parameters == %
delay_spread = 5.4e-6;
fc = 2e9;
v = 15; % speed of receiver
fD = (v/(physconst('LightSpeed')))*fc;
B = 1e6;
Noise = 2.07e-20*B; % in [W]
fs = 1e6;
Ts = 1/fs; % Sample time
fDTs = fD*Ts; % Normalized Doppler frequency
E = 1;
N = 100;
L = ceil(delay_spread/Ts);
Ncp = L-1;
tau = (0:1:L-1); % Path delays in samples
P = (1/length(tau))*ones(1,length(tau)); % Flat power delay profile
M = 30; % Number of time samples
if (N + Ncp)*fDTs*10 > 1 % Check constraint
    error('N or Ncp to large')
end
if N*Ts > 1/(2*fD) % Check if N*Ts is smaller than coherence time
    error('N too large')
end
To = (N+Ncp)*Ts;
n = 1:N+Ncp;
%z = (m-1)*To;
%% == Simulation Parameters == %%
%snr_db = 1:1:30; % snr in db
%BER = zeros(1, length(snr_db)); % Create bit error vector
%maxNumErrs = 50; % Simulate until we get atleast this amount of errors
%maxNum = 1e4; % Stop if maxNum of bits have been simulated
% == Simulation  == %
%for i=1:length(snr_db)
%totError = 0; %Number of errors
%num = 0; % Number of processed bits
%while((totError < maxNumErrs) && (num < maxNum))

% == Transmitter  == %
symbs = zeros(N,M);
bits_s = zeros(2*N,M);
y = zeros(N,M);

for m = 1:M
    [symbs(:,m), bits_s(:,m)] = qpsk_sequence(N,E);
    symbs_ifft = sqrt(N)*ifft(symbs(:,m),N);
    z = [symbs_ifft(end-Ncp+1:end);symbs_ifft];
    
    % == Channel  == %
    [z, h] = Fading_Channel(z, tau, fDTs, P);
    z = z(1:end-L+1); % Remove delayed symbols
    
    %c_output = awgn(z,i,'Measured');
    z = z*sqrt(10^(-10.1)); % Pathloss
    c_output = z + (randn(length(z),1) +1i*rand(length(z),1))*sqrt(Noise);   % add noise
    
    % == Reciever  == %
    r2 = c_output(Ncp+1:end); %remove cp
    
    ch = fft(h(1,:),N); % FFT of channel
    
    y(:,m) = fft(r2,N)/sqrt(N); % FFT of received symbols
    y(:,m) =y(:,m)./ch.';  % zero forcing equalizer with the first tap of the channel frequency response.
    
end

% symb -> bits
%bits_r = nrSymbolDemodulate(y,'QPSK','DecisionType','Hard');
% calculate bit error rate
%bit_errors = biterr(bits_s,bits_r);
%BER(i) = bit_errors/(2*N);
%BER = bit_errors/(2*N);


scatterplot(y(1,:)); % scatterplot for the first subcarrier
scatterplot(y(25,:)); 
scatterplot(y(50,:));
scatterplot(y(100,:)); 


% figure
% for i= 1:4
%     subplot(2,2,i)
%     scatter(real(y(i,:)),imag(y(i,:)));
% end

%figure
%plot(snr_db,BER,'-bo')
%grid on
%xlabel('SNR [dB]');
%ylabel('BER');
