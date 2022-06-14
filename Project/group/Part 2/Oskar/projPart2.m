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
pathLoss = 10^(-101/10); % -101dB = 10*log10(Prx/Ptx)
L = ceil(t_ds/Ts); % L*T must be greater or equal to the delay spread, also an integer.
tau = [0 1 2 3 4 5]; % The path delay in samples
Ncp = L-1; % minimum number of cyclic prefix L-1, you can have more but is unnecessary.
N = 2^7; % number of samples, have it so that it is in a power of 2 (easy fft)
% we want the case when (N + Ncp)fDTs << 1 (basically much lower than coherence time)
if (N + Ncp)*fdTs > 0.1
    error('N is too large')
end

%% Simualation
mode = 0; % 0 is for no fading channel, 1 is for using fading channel
M = 30; % number of OFDM symbols we want to transmit
Ptx = 20:5:90; % create a vector of different transmission power and space 
% it evenly in dB scale to get a smooth curve (more points make it edgy due to small variances)
E = 10.^(Ptx/10)*Ts; % E = Ptx*Ts [Joule = Watt*Seconds]
SNR = nan(1,length(E));
SNR_subcarrier = nan(N,length(E));
SER = nan(1,length(E));
SERtheo = nan(1,length(E));
for i = 1:length(E) % (loop with varying E values to plot SER over SNR)
    P = 1/L*ones(1,L); % PDP, divide equally on all taps
    QPSKconst = [1+1i 1-1i -1+1i -1-1i]*(sqrt(E(i))/sqrt(2)); % create the qpsk constellation (doesn't matter if gray)
    OFDMsymbols = zeros(N,M);
    QPSKsymbols = zeros(N,M);
    msg = zeros(N,M);
    for m = 1:M
        % transmitter part
        msg(:,m) = randi([1 4],N,1);
        QPSKsymbols(:,m) = QPSKconst(msg(:,m)); % the channel input of L channels
        z = sqrt(N)*ifft(QPSKsymbols(:,m),N); % apply the ifft to get signal in time domain
        signal = [z(end-Ncp+1:end);z]; % add cyclic prefix
        
        % Channel part
        switch mode
            case 0 % here we only have awgn
                RXsignal = signal.';
                C = [1 zeros(1,length(RXsignal)-Ncp-1)]; % this is how it looks but we don't need it
                awgn = sqrt(N0/2)*(randn(size(RXsignal))+1i*randn(size(RXsignal)));
                RXsignal = RXsignal*sqrt(pathLoss) + awgn; % add pathloss and complex AWGN samples on recieved signal
                RXsignal = RXsignal(Ncp+1:end); % remove cyclic prefix
                y = fft(RXsignal)/sqrt(N);
                s = y; % for AWGN channel we dont need to equalize
            case 1
                % generate Rayleigh fading channel (needed for second parts)
                [r, h] = Fading_Channel(signal, tau, fdTs, P); % simulate channel with premade function
                r = r(1:end-L+1); % remove the delayed symbols
                RXsignal = r;
                hm=h(1,:); % time domain response since (N + Ncp)fDTs << 1.
                C = fft(hm,N);
                awgn = sqrt(N0/2)*(randn(size(RXsignal))+1i*randn(size(RXsignal)));
                RXsignal = RXsignal*sqrt(pathLoss) + awgn; % add pathloss and complex AWGN samples on recieved signal 
                RXsignal = RXsignal(Ncp+1:end); % remove cyclic prefix
                y = fft(RXsignal)/sqrt(N);
                s = y./C.'; % compute the equalization (remove effects of the channel such as rotation)
            otherwise
                error('set mode to correct mode')
        end
        OFDMsymbols(:,m) = s; % store the symbols in this matrix.
    end
    % Reciever part
    recievedOFDM = reshape(OFDMsymbols,[N*M 1]);
    % check minimum distance
    distance = abs(repmat(recievedOFDM,1,4) - repmat(QPSKconst, length(recievedOFDM), 1)).^2; %compute the distance to each possible symbol
    [~, idx] = min(distance, [], 2); % find the constellation index for symbol alternative at minimum distance for every recieved symbol
    % calculate SER (alternative way is via symerr function)
    originalMSG = reshape(msg,[N*M 1]);
    SER(i) = length(find(originalMSG ~= idx))/length(idx);
    SNR(i) = 10*log10(E(i)*pathLoss/(N0)); % theoretical SNR in dB
    SERtheo(i) = 2*qfunc(sqrt(E(i)*pathLoss/(N0))); % theoretical SER
    SNR_subcarrier(:,i) = 10*log10(mean(abs(OFDMsymbols).^2,2)/(N0)); % average subcarrier SNR in dB
end
avgSNR = mean(SNR_subcarrier,1); % average SNR in dB
%% Plot SER vs SNR
switch mode
    case 0
        figure()
        semilogy(SNR,SER,'r'), grid on, hold on
        semilogy(SNR,SERtheo,'b')
        ylim([0.001 1]), xlim([-5 20])
        xlabel('SNR [dB]'), ylabel('SER')
        legend('simulated SER','theoretical AWGN SER')
        title('SER vs SNR for AWGN channel')
        saveas(gcf,'SERvsSNR_awgn','epsc')
    case 1
        figure()
        semilogy(SNR,SER,'r'), grid on, hold on
        semilogy(SNR,SERtheo,'b')
        ylim([0.001 1]), xlim([-5 30])
        xlabel('SNR [dB]'), ylabel('SER')
        title('SER vs SNR for fading channel')
        legend('simulated SER','theoretical AWGN SER')
        saveas(gcf,'SERvsSNR_fading','epsc')
end

%% Plot some scatterplots
scatterplot(OFDMsymbols(2,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 2')
saveas(gcf,'ScatterN2','epsc')
scatterplot(OFDMsymbols(6,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 6')
saveas(gcf,'ScatterN6','epsc')
scatterplot(OFDMsymbols(29,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 29')
saveas(gcf,'ScatterN29','epsc')
scatterplot(OFDMsymbols(53,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 53')
saveas(gcf,'ScatterN53','epsc')
scatterplot(OFDMsymbols(78,:)), grid on
title('Scatterplot of recieved symbols on subcarrier N = 78')
saveas(gcf,'ScatterN78','epsc')
%% Simulate with reducing Ncp lengths
E = 1e-1;
snr = 10*log10(E*pathLoss/N0);
P = 1/L*ones(1,L); % PDP, divide equally on all taps
QPSKconst = [1+1i 1-1i -1+1i -1-1i]*(sqrt(E)/sqrt(2)); % create the qpsk constellation
OFDMsymbols = zeros(N,M);
QPSKsymbols = zeros(N,M);
msg = zeros(N,M);
Ncp = L-1:-1:0;
SER = nan(1,length(Ncp));
for i = 1:length(Ncp)
    for m = 1:M
        % transmitter part
        msg(:,m) = randi([1 4],N,1);
        QPSKsymbols(:,m) = QPSKconst(msg(:,m)); % the channel input of L channels
        z = sqrt(N)*ifft(QPSKsymbols(:,m),N); % apply the ifft to get signal in time domain
        signal = [z(end-Ncp(i)+1:end);z]; % add cyclic prefix
        % channel part
        [r, h] = Fading_Channel(signal, tau, fdTs, P); % simulate channel with premade function
        r = r(1:end-L+1); % remove the delayed symbols
        RXsignal = r;
        hm=h(1,:); % time domain response since (N + Ncp)fDTs << 1.
        C = fft(hm,N);
        awgn = sqrt(N0/2)*(randn(size(RXsignal))+1i*randn(size(RXsignal)));
        RXsignal = RXsignal*sqrt(pathLoss) + awgn; % add complex AWGN samples on recieved signal
        RXsignal = RXsignal(Ncp(i)+1:end); % remove cyclic prefix
        y = fft(RXsignal)/sqrt(N);
        s = y./C.'; % compute the equalization (remove effects of the channel such as rotation)
        OFDMsymbols(:,m) = s; % store the symbols in this matrix.
    end
    % Reciever part
    recievedOFDM = reshape(OFDMsymbols,[N*M 1]);
    % check minimum distance
    distance = abs(repmat(recievedOFDM,1,4) - repmat(QPSKconst, length(recievedOFDM), 1)).^2; %compute the distance to each possible symbol
    [~, idx] = min(distance, [], 2); % find the constellation index for symbol alternative at minimum distance for every recieved symbol
    % calculate SER
    originalMSG = reshape(msg,[N*M 1]);
    SER(i) = length(find(abs(originalMSG-idx)))/length(idx);
    scatterplot(recievedOFDM)
    axis([-1.5e-5 1.5e-5 -1.5e-5 1.5e-5])
    title(['Scatterplots when N_{cp} = ' num2str(Ncp(i))])
    saveas(gcf,['scatter_Ncp' num2str(Ncp(i))],'epsc')
end
figure()
plot(Ncp,SER,'rx'), grid on 
title('SER vs length of N_{cp} with fixed SNR')
xlabel('N_{cp} lengths'), ylabel(['SER for SNR= ' num2str(snr) '[dB]'])
saveas(gcf,'SER_Ncp','epsc')
