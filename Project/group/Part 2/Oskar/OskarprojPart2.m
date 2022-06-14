%% Part 2
% settings
clc, clear all
addpath('./functions')

% initiate parameters
f_c = 2e9; % 2GHz frequency carrier
BW = 1e6; % 1MHz bandwidth
Ts = 1/BW; % symbol time
% f_s = 1/T_s; % sampling frequency
N0 = 2.07e-20*BW; % noise power (times the bandwidth due to the units given)
v = 15; % velocity in m/s
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
fdTs = f_D*Ts; % Normalized Doppler frequency
t_ds = 5.4e-6; % delay spread
pathLoss = 10^(-101/10); % -101dB = 10*log10(Pr/Pt)
QPSKconst = [1+1i 1-1i -1+1i -1-1i]/sqrt(2); % create the qpsk constellation
L = ceil(t_ds/Ts); % L*T must be greater or equal to the delay spread, also an integer.
Ncp = L-1; % minimum number of cyclic prefix L-1, you can have more but is unnecessary.
N = 256; % number of samples, have it so that it is in a power of 2 (easy fft)
tau = [0 1 2 3 4 5];
Ptx = 0.1; % average transmit power
P = [Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L Ptx/L]; % PDP, divide equally on all taps
mode = 0; % 0 is for no fading channel, 1 is for using fading channel

% we want the case when (N + Ncp)fDTs << 1 (basically much lower than coherence time)
if (N + Ncp)*fdTs > 0.1 || (N+Ncp)*Ts > 1/(2*f_D*10)
    error('N is too large')
end

%% transmitter part
s = QPSKconst(randi([1 4],N,1)); % the channel input of L channels
z = sqrt(N)*ifft(s); % apply the ifft to get signal in time domain
signal = [z(end-Ncp+1:end) z]; % add cyclic prefix

%% generate Rayleigh fading channel (needed for second parts)
[r, h] = Fading_Channel(signal, tau, fdTs, P); % simulate channel with premade function
r = r(1:end-L+1); % remove the delayed symbols
%% Reciever part
switch mode
    case 0 % here we only have awgn
        RXsignal = signal';
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
scatterplot(s)

% check minimum distance
distance = abs(repmat(s,1,4) - repmat(QPSKconst, length(s), 1)).^2; %compute the distance to each possible symbol
[~, idx] = min(distance, [], 2); % find the constellation index for symbol alternative at minimum distance for every recieved symbol
