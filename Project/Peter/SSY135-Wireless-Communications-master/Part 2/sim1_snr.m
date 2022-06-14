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
E = linspace(0.0001,0.03,100);
N = 20;
L = ceil(delay_spread/Ts);
Ncp = L-1;
tau = (0:1:L-1); % Path delays in samples
P = (1/length(tau))*ones(1,length(tau)); % Flat power delay profile
M = 30; % Number of time samples
PL = 10^(-10.1);
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
pow = nan(length(E),1);
SER = nan(length(E),1);
% == Simulation  == %
for i=1:length(E)
    % == Transmitter  == %
    symbs = zeros(N,M);
    bits_s = zeros(2*N,M);
    y = zeros(N,M);
    
    for m = 1:M
        [symbs(:,m), bits_s(:,m)] = qpsk_sequence(N,E(i));
        symbs_ifft = sqrt(N)*ifft(symbs(:,m),N);
        z = [symbs_ifft(end-Ncp+1:end);symbs_ifft];
        
        % == Channel  == %
        [z, h] = Fading_Channel(z, tau, fDTs, P);
        z = z(1:end-L+1); % Remove delayed symbols
        
        z = z*sqrt(10^(-10.1)); % Pathloss
        c_output = z + (randn(length(z),1) +1i*rand(length(z),1))*sqrt(Noise);   % add noise
        
        % == Reciever  == %
        r2 = c_output(Ncp+1:end); %remove cp
        
        ch = fft(h(1,:),N); % FFT of channel
        
        y(:,m) = fft(r2,N)/sqrt(N); % FFT of received symbols
        y(:,m) =y(:,m)./ch.';  % zero forcing equalizer with the first tap of the channel frequency response.
    end
    const = [1+1j,1-1j,-1+1j,-1-1j]*sqrt(0.5*E(i));
    received_ofdm_syms = reshape(y,[N*M 1]);
    metric = abs(repmat(received_ofdm_syms,1,4) - repmat(const, length(received_ofdm_syms), 1)).^2; % compute the distance to each possible symbol
    [~, indx_r] = min(metric, [], 2); % find the closest for each received symbol
    msg_r = nan(length(indx_r),1); % All ofdm symbs in one vector
    for j=1:length(indx_r)
        msg_r(j) = const(indx_r(j));
    end
    msg_r = reshape(msg_r,[N M]);
    %  == Calculate symbol error rate == %
    [number,ratio] = symerr(round(msg_r,4),round(symbs,4));
    SER(i) = ratio;
    % == Power of signal == %
    pow(i) = mean(mean(abs(y).^2))*N/Noise; %averaging over all subcarriers, for large number of OFDM symbols is it the same than for single sucarriers.
    %pow(i) = qfuncinv(ratio/2)^2;
end

Snr = @(x)((x)*(10^(-10.1)))*N/Noise;   %x or x^2 depends on implementation of QPSK_sequence: 
                                        %x for sqrt(E)*symbol; x^2 for E*symbol
figure
semilogy(10*log10(E/Noise),SER,'b-'), hold on
xlabel('SNR [dB]'), ylabel('Symbol error rate')
figure
plot(E,10*log10(pow),'b-'), hold on
plot(E,10*log10(Snr(E)),'r--');
legend('Simulated SNR','Theory SNR','location','northwest')