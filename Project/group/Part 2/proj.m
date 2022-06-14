%% Proj part 2
clc, clear all

% set given parameters for both filter- & spectrum method.
f_c = 2e9; % 2GHz frequency carrier
v = 15; % 30km/h velocity
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
f_s = 1e6; % 0.1 ms sample interval
T_s = 1/f_s; % sampling frequency
fdTs=f_D*T_s
N0=2*2.07e-14;
tau= 4.5e-6;

Ncp=tau/T_s;

