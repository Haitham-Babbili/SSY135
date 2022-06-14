%% Proj part 1
clc, clear all, close all
% 1) simulate a frequency-flat Rayleigh fading gain process.

% set given parameters for both filter- & spectrum method.
f_c = 2e9; % 2GHz frequency carrier
v = 30/3.6; % 30km/h velocity
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
T_s = 0.1e-3; % 0.1 ms sample interval
maximumToleratedTs = 1/(2*f_D); % calculate the max Ts for no aliasing
if T_s > maximumToleratedTs % this must be satisfied
    error('The chosen sample interval time is too small, aliasing will occur')
end

%% filter method:
N_s = 10000; % choose amount of samples (idk)
N = 1000; 

% initiate g vector and store g(nTs) values in it
g=abs(-N*T_s:T_s:T_s*N);

g(:)=besselj(1/4,2*pi*f_D*g(:))./g(:).^(1/4);
g(N+1) = pi*f_D^(1/4)/gamma(5/4);


% make sure g is of unit energy
g=g/norm(g);

% generate samples from zero-mean unit-variance complex Gaussian dist, we
% must divide by sqrt(2) since adding a N(0,1) with a N(0,1) gives a N(0,2)
x = (randn(N_s,1) + 1i*randn(N_s,1))/sqrt(2);

% convolve x with g and remove the transients (can be done using 'same')
c_filter = conv(x,g,'same');

% check the to see if it's correctly made distribution
% plot(c_filter), grid on, axis equal % looks like a complex gaussian with zero mean 
% ExpectedValues = [mean(c_filter) var(c_filter)] % it's close to 0 and 1

%% Spectrum method

f_s = 1/T_s; % sampling rate
N_s = 2000; % number of samples
f = -f_D:f_s/N_s:f_D; % create a discrete frequency vector within bandwidth
S_c = zeros(1,length(f)); % initate S_c as vector
% calculate S_c as given by (3)
for i = 1:length(f)
    S_c(i) = (1/(pi*f_D))*(1/sqrt(1-(f(i)/f_D)^2));
end
% S_c(1) equals to infinity due to the term 1/sqrt(1-(-fD/fD)^2) goes to 
% infinity, therefore this must be fixed.
S_c(1) = S_c(end);

% choose G as the sqrt of S
G = sqrt(S_c);

% create G such that the last half of G is the start of G then a zero 
% padding and then the first half of G to create a "periodic" G
padding = zeros(1,N_s-length(S_c)); % make sure this is right size
G = [G(length(G)/2:end) padding G(1:length(G)/2)]; % now G_p should be of length N_s

% Generate samples from zero-mean, complex gaussian dist
var_x = N_s^2/(norm(G)^2); % calculate the variance in order to remove it
X = sqrt(var_x)*(randn(N_s,1) + 1i*randn(N_s,1)); 
C = X'.*G;
c_spectrum = ifft(C);
% variance = var(c_spectrum,0,'all'); % it's close to 1

%%[PSD_filter,psdf] = pwelch(c_filter,[],[],[],f_s,'centered');
PSD_filter = PSD_filter /(sum(PSD_filter)*(f_s/N_s));
[PSD_spectrum,psdf2] = pwelch(c_spectrum,[],[],[],f_s,'centered');
PSD_spectrum = PSD_spectrum /(sum(PSD_spectrum)*(f_s/N_s));
% theoretical PSD
f_s = 1/T_s;
f = (-f_s:f_s/N_s:f_s);
%f=[-f,fliplr(f)];
S = 1./(pi*f_D*(sqrt(1-(f/f_D).^2)));
S(abs(f)>=f_D) = 0;

%S=S+fliplr(S);

%PSD_Xaxis= linspace(-f_D,f_D,2000);

figure()
plot(f,pow2db(S))
xlim([-155,155])
hold on
plot(psdf,pow2db(PSD_filter))

plot(psdf2,pow2db(PSD_spectrum))
legend('theory','filter','spectrum');


figure()
plot(f(2:end),c_spectrum)


figure()
plot(f(2:end),c_filter')


N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

plot(freq,10*log10(psdx))
%% create rician
k_c = 1;

if k_c == 0
    c = c_filter;
else
    c = c_filter + k_c;
end