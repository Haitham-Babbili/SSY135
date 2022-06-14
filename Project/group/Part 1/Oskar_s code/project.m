%% Proj part 1
clc, clear all
% 1) simulate a frequency-flat Rayleigh fading gain process.
addpath('./functions')
% set given parameters for both filter- & spectrum method.
f_c = 2e9; % 2GHz frequency carrier
v = 30/3.6; % 30km/h velocity
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
T_s = 0.1e-3; % 0.1 ms sample interval
maximumToleratedTs = 1/(2*f_D); % calculate the max Ts for no aliasing
if T_s > maximumToleratedTs % this must be satisfied
    error('The chosen sample interval time is too small, aliasing will occur')
end


%% Create the PSD plots
% theoretical PSD
f_s = 1/T_s;
f = (-80:0.08:80);
S = 1./(pi*f_D*(sqrt(1-(f/f_D).^2)));
S(abs(f)>=f_D) = 0; % remove parts outside the bandwidth
% simulation PSD
f_simulation = linspace(-f_D,f_D,N_s);
PSD_spectrum = pwelch(c_spectrum,1000,[],N_s);
PSD_spectrum = PSD_spectrum /(sum(PSD_spectrum)*(f_s/N_s));
PSD_filter = pwelch(c_filter,1000,[],N_s);
PSD_filter = PSD_filter /(sum(PSD_filter)*(f_s/N_s));
figure()
plot(f,10*log10(S),'k'), hold on
plot(f_simulation,10*log10(PSD_spectrum),'r')
plot(f_simulation,10*log10(PSD_filter),'g')
%% 
K_c = [0 1 10];
xAxis = 0:20;
for i = 1:3
    theoreticalRician = makedist('Rician','s',K_c(i),'sigma',1);
    spectrumDistribution = fitdist(abs(c_spectrum'+K_c(i)),'Rician');
    filterDistribution = fitdist(abs(c_filter+K_c(i)),'Rician');
end
%% Part 2
N = 300;
M = 64;
L = [1 2 3];
K_c = [0 1 10];
fDTs = [0.1 0.005];

for i = 1:3
    for j = 1:3
        for k = 1:2
            
        end
    end
end



