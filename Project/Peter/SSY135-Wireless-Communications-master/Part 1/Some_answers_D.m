% Given parameters
% V is in km/hr
% v is is m/s
% c= speed of light
% fc= Carrier Frequency
% fD = Doppler Frequency, in Hertz
% Ts=Sampling time
V=30;
v=(1000/3600) .*V;
fc=2e9;
c=3e8;
lamda=fc/c; % wavelength in metre
fD=(v*fc)./c ; %fD is maximum when cos0 =1
% fD= 55.5556 Hertz
% To avoid aliasing (fs=1/Ts) > 2*fD

% Highest fs is at lowest Ts and highest Ts is at lowest fs 

fs_l=2*fD; % fs_h is lowest frequency to avoid Aliasing, in hertz
%fs_l = 111.11Hertz, approximately

Ts_h= (1/fs_l); % Ts_h is highest Ts, in seconds
%Ts_h = 9ms

%Given Ts=0.1ms, means fs is 10kHz. Then, this still meets the requirement.
Ts=0.0001;
fs =1/Ts; % in hertz
% Coherence time (ct) = 1/Doppler Spread
ct= 1/fD;
% Coherence time is 18ms
