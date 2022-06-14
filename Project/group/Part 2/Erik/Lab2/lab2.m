clc
clear all
close all

%% Parameters
% Constelation

M=4; %decides the constelation

%N % lecture slides


fc= 2*10^9;
fs= 1*10^6;

string_length=5; % length of data string
bits= randsrc(1,string_length,[0 1]);
modulation= get_modulation(M);

CP=2;
N=1; %the amount of channels
fT=1/N; % spacing
first_zero =1/(N+CP); % distance to first zero



%parameters for fading_channel.p

P = [0.5 0.5]; % Power delay profile 
tau = [0 4]; % Path delays in samples 
%M = 300; % Number of time samples
fd = 100; % Doppler frequency 
Ts = 1e-6; % Sample time 
fdTs = 1e-4; % Normalized Doppler frequency



%% Transmitter 

msg=buffer(bits,log2(M)).'; %bits converted to  messages

%verifying
% msgr=msg.'
% received_bits= msgr(:)'

m_idx= bi2de(msg,'left-msb')+1; %messages converted to symbols
sym=modulation(m_idx).';

%verifying
% msg_received= sym2msg(sym,modulation,M)
% msgr=msg_received.'
% received_bits= msgr(:)'

%  coded = repmat(sym, 1, 3)
%  coded = reshape(coded, numel(coded), 1)

% if mod(numel(serial) / channels, 1) ~= 0
%     serial = [serial ; zeros(channels - mod(numel(serial), channels), 1)];
% end
% parallel = reshape(serial, channels, numel(serial) / channels);
if mod(length(sym) / N,1) ~= 0
    sym = [sym ; zeros(N - mod(length(serial),N),1)];
end



%% Channel
% Pass the input through the channel
%[channel_output, tap] = Fading_Channel(channel_input, tau, fdTs,P);


%% Reciever


msg_received= sym2msg(sym,modulation,M);
msgr=msg_received.';
received_bits= msgr(:)';
