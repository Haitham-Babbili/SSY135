%% Part 2
% settings
clc, clear all
addpath('./functions')

% initiate parameters
f_c = 2e9; % 2GHz frequency carrier
BW = 1e6; % 1MHz bandwidth
Ts = 1/BW; % symbol time
% f_s = 1/Ts; % sampling frequency
N0 = 2.07e-20*BW; % noise power (times the doublesided-bandwidth due to the unit)
v = 15; % velocity in m/s
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
fdTs = f_D*Ts; % Normalized Doppler frequency
t_ds = 5.4e-6; % delay spread
pathLoss = 10^(-101/10); % -101dB = 10*log10(Pr/Pt)

L = ceil(t_ds/Ts); % L*T must be greater or equal to the delay spread, also an integer.
tau = [0 1 2 3 4 5]; % The path delay in samples
Ncp = L-1; % minimum number of cyclic prefix L-1, you can have more but is unnecessary.
N = 2^7; % number of samples, have it so that it is in a power of 2 (easy fft)
% we want the case when (N + Ncp)fDTs << 1 (basically much lower than coherence time)
if (N + Ncp)*fdTs > 0.1
    error('N is too large')
end

%% Simualation
mode = 1; % 0 is for no fading channel, 1 is for using fading channel
M = 200; % number of OFDM symbols we want to transmit
E = linspace(1e-6,1e-2,200); % E = Ptx*Ts
SNR = nan(1,length(E));
SNR_subcarrier = nan(N,length(E));
SER = nan(1,length(E));
SERtheo = nan(1,length(E));
for i = 1:length(E) % (loop with varying E values to plot SER over SNR)
    Ptx = E(i)/Ts; % average transmit power (If we use this in P, the SER gets super smooth instead of when its 1)
    P = [1/L 1/L 1/L 1/L 1/L 1/L]; % PDP, divide equally on all taps
    QPSKconst = [1+1i 1-1i -1+1i -1-1i]*(sqrt(E(i))/sqrt(2)); % create the qpsk constellation
    SNR(i) = 10*log10(E(i)/N0);
    SERtheo(i) = 2.*qfunc(sqrt(10^(2*E(i)/10))); % this is not working
    OFDMsymbols = zeros(N,M);
    QPSKsymbols = zeros(N,M);
    for m = 1:M
        % transmitter part
        msg = randi([1 4],N,1);
        QPSKsymbols(:,m) = QPSKconst(msg); % the channel input of L channels
        z = sqrt(N)*ifft(QPSKsymbols(:,m),N); % apply the ifft to get signal in time domain
        signal = [z(end-Ncp+1:end);z]; % add cyclic prefix
        
        % generate Rayleigh fading channel (needed for second parts)
        [r, h] = Fading_Channel(signal, tau, fdTs, P); % simulate channel with premade function
        r = r(1:end-L+1); % remove the delayed symbols
        % Reciever part
        switch mode
            case 0 % here we only have awgn
                RXsignal = signal.';
                C = zeros(1,length(RXsignal)-Ncp);
                C(1) = 1;
                awgn = N0/sqrt(2)*(randn(size(RXsignal))+1i*randn(size(RXsignal)));
                RXsignal = RXsignal*sqrt(pathLoss) + awgn; % add complex AWGN samples on recieved signal
                RXsignal = RXsignal(Ncp+1:end); % remove cyclic prefix
                y = fft(RXsignal)/sqrt(N);
                s = y; % for AWGN channel we dont need to equalize
            case 1
                RXsignal = r;
                hm=h(1,:); % time domain response since (N + Ncp)fDTs << 1.
                C = fft(hm,N);
                awgn = sqrt(N0/N)*(randn(size(RXsignal))+1i*randn(size(RXsignal)));
                RXsignal = RXsignal*sqrt(pathLoss) + awgn; % add complex AWGN samples on recieved signal
                RXsignal = RXsignal(Ncp+1:end); % remove cyclic prefix
                y = fft(RXsignal)/sqrt(N);
                s = y./C.'; % compute the equalization (remove effects of the channel such as rotation)
            otherwise
                error('set mode to correct mode')
        end

        if length(s) ~= N
            error('s is not of length N')
        end
        
        OFDMsymbols(:,m) = s; % store the symbols in this matrix.
    end
        % check minimum distance
        recievedOFDM = reshape(OFDMsymbols,[N*M 1]);
        distance = abs(repmat(recievedOFDM,1,4) - repmat(QPSKconst, length(recievedOFDM), 1)).^2; %compute the distance to each possible symbol
        [~, idx] = min(distance, [], 2); % find the constellation index for symbol alternative at minimum distance for every recieved symbol
        
        % check SER
        recievedSymbols = nan(length(idx),1); % Store the recieved OFDM symbols in one vector
        for k=1:length(idx)
        recievedSymbols(k) = QPSKconst(idx(k));
        end
        recievedSymbols = reshape(recievedSymbols,[N M]);
        [~,ratio] = symerr(round(recievedSymbols,4),round(QPSKsymbols,4));
        SER(i) = ratio;
        SNR_subcarrier(:,m) = mean(abs(s).^2,2)*N/N0; % average subcarrier SNR
end

%% Plot SER vs SNR
figure()
semilogy(SNR,SER,'r'), grid on, hold on
semilogy(SNR,SERtheo,'b')
xlabel('SNR [dB]'), ylabel('SER')
%% Plot some scatterplots
scatterplot(OFDMsymbols(2,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 2')
scatterplot(OFDMsymbols(6,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 6')
scatterplot(OFDMsymbols(12,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 12')
scatterplot(OFDMsymbols(29,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 29')
scatterplot(OFDMsymbols(40,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 40')
scatterplot(OFDMsymbols(53,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 53')
scatterplot(OFDMsymbols(78,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 78')