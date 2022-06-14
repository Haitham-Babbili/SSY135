% Part 2
% settings
clc, clear all
addpath('./functions')

% initiate parameters
f_c = 2e9; % 2GHz frequency carrier
BW = 1e6; % 1MHz bandwidth
Ts = 1/BW; % symbol time
% f_s = 1/T_s; % sampling frequency
N0 = 2.07e-20*BW; % noise power (times the doublesided-bandwidth due to the unit)
v = 15; % velocity in m/s
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
fdTs = f_D*Ts; % Normalized Doppler frequency
t_ds = 5.4e-6; % delay spread
pathLoss = 10^(-101/10); % -101dB = 10*log10(Pr/Pt)
tau = [0 1 2 3 4 5];

L = ceil(t_ds/Ts); % L*T must be greater or equal to the delay spread, also an integer.
Ncp = L-1; % minimum number of cyclic prefix L-1, you can have more but is unnecessary.
N = 2^7; % number of samples, have it so that it is in a power of 2 (easy fft)
% we want the case when (N + Ncp)fDTs << 1 (basically much lower than coherence time)
if (N + Ncp)*fdTs > 0.1
    error('N is too large')
end

%% Simualation
mode = 1; % 0 is for no fading channel, 1 is for using fading channel
M = 100; % number of OFDM symbols we want to transmit
SNR = zeros(1,1000);
avgSER = zeros(1,1000);
E = linspace(1e-5,1e-1,1000); % E = Ptx*Ts
for j = 1:length(E) % (loop with varying E values to plot SER over SNR)
    Ptx = E(j)/Ts; % average transmit power
    P = [Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L]; % PDP, divide equally on all taps
    QPSKconst = [1+1i 1-1i -1+1i -1-1i]*(E(j)/sqrt(2)); % create the qpsk constellation
    SNR(j) = 10*log10(E(j)/N0);
    OFDMsymbols = zeros(N,M);
    SER = zeros(1,M);
    for i = 1:M
        
        % transmitter part
        msg = randi([1 4],N,1);
        s = QPSKconst(msg); % the channel input of L channels
        z = sqrt(N)*ifft(s); % apply the ifft to get signal in time domain
        signal = [z(end-Ncp+1:end) z]; % add cyclic prefix
        
        % generate Rayleigh fading channel (needed for second parts)
        [r, h] = Fading_Channel(signal, tau, fdTs, P); % simulate channel with premade function
        r = r(1:end-L+1); % remove the delayed symbols
        % Reciever part
        switch mode
            case 0 % here we only have awgn
                RXsignal = signal.';
                C = ones(1,length(RXsignal)-Ncp);
            case 1
                RXsignal = r;
                hm=h(1,:); % time domain response since (N + Ncp)fDTs << 1.
                C = fft(hm,N);
            otherwise
                error('set mode to correct mode')
        end
        
        awgn = N0/sqrt(2)*(randn(size(RXsignal))+1i*randn(size(RXsignal)));
        RXsignal = RXsignal*pathLoss + awgn; % add complex AWGN samples on recieved signal
        RXsignal = RXsignal(Ncp+1:end); % remove cyclic prefix
        
        if length(RXsignal) ~= N
            error('RXsignal is not of length N')
        end
        y = fft(RXsignal);
        s = y./C.'; % compute the equalization (remove effects of the channel such as rotation)
        OFDMsymbols(:,i) = s;
        % check minimum distance
        distance = abs(repmat(s,1,4) - repmat(QPSKconst, length(s), 1)).^2; %compute the distance to each possible symbol
        [~, idx] = min(distance, [], 2); % find the constellation index for symbol alternative at minimum distance for every recieved symbol
        % check SER
        SER(i) = length(find(abs(msg-idx)))/length(idx);
    end
    avgSER(j) = mean(SER);
end
semilogy(SNR,avgSER)
%% Plot some scatterplots
% scatterplot(OFDMsymbols(1:50,1))
% scatterplot(OFDMsymbols(48:90,3))
% scatterplot(OFDMsymbols(75:end,4))
% scatterplot(OFDMsymbols(66:120,6))
% scatterplot(OFDMsymbols(6:72,7))
% scatterplot(OFDMsymbols(32:89,8))
% scatterplot(OFDMsymbols(20:77,10))