%% Proj part 1
clc, clear all
% 1) simulate a frequency-flat Rayleigh fading gain process.
addpath('./functions')
% set given parameters for both filter- & spectrum method.
f_c = 2e9; % 2GHz frequency carrier
v = 30/3.6; % 30km/h velocity
f_D = v/physconst('LightSpeed')*f_c; % calculate the doppler frequency
T_s = 0.1e-3; % 0.1 ms sample interval
f_s = 1/T_s; % sampling frequency
maximumToleratedTs = 1/(2*f_D); % calculate the max Ts for no aliasing
if T_s > maximumToleratedTs % this must be satisfied
    error('The chosen sample interval time is too small, aliasing will occur')
end

% choose a reasonable N_s, not too small
N_s = 10000;

% create the Rayleigh flat fading channels.
c_spectrum = spectrumMethod(f_D,T_s,N_s,0);
expectedStuff = [mean(c_spectrum) var(c_spectrum)];
c_filter = filterMethod(f_D,T_s,N_s,0);
expectedStuff2 = [mean(c_filter) var(c_filter)];

%% Create the PSD plots
% theoretical PSD
f = (-60:f_s/N_s:60);
S = 1./(pi*f_D*(sqrt(1-(f/f_D).^2)));
S(abs(f)>=f_D) = 0; % remove parts outside the bandwidth
% simulation PSD
[PSD_spectrum,psdf1] = pwelch(c_spectrum,[],[],[],f_s,'centered');
PSD_spectrum = PSD_spectrum /(sum(PSD_spectrum)*(f_s/N_s));
[PSD_filter,psdf2] = pwelch(c_filter,[],[],[],f_s,'centered');
PSD_filter = PSD_filter /(sum(PSD_filter)*(f_s/N_s));
figure(1)
plot(f,pow2db(S)), hold on, grid on, xlim([-155,155])
plot(psdf1,pow2db(PSD_spectrum))
plot(psdf2,pow2db(PSD_filter))
title('PSD plot of Rayleigh flat fading for both methods and theoretical')
xlabel('Frequency [Hz]'), ylabel('PSD [dB/Hz]')
legend('theoretical','filter method','spectrum method');
%% PDF and CDF 
K_c = [0 1 10];
xAxis = 0:0.01:20;
% initialize to store every iteration
PDF_theory = cell(3,1);
PDF_spectrum = cell(3,1);
PDF_filter = cell(3,1);
CDF_theory = cell(3,1);
CDF_spectrum = cell(3,1);
CDF_filter = cell(3,1);
for i = 1:3
    % create the Rician distributions
    theoreticalRician = makedist('Rician','s',K_c(i),'sigma',1);
    spectrumDistribution = fitdist(abs(spectrumMethod(f_D,T_s,N_s,K_c(i))'),'Rician');
    filterDistribution = fitdist(abs(filterMethod(f_D,T_s,N_s,K_c(i))),'Rician');
    
    % Get the pdf for these distributions
    PDF_theory{i} = pdf(theoreticalRician,xAxis);
    PDF_spectrum{i} = pdf(spectrumDistribution,xAxis);
    PDF_filter{i} = pdf(filterDistribution,xAxis);

    % Get the cdf for these distributions
    CDF_theory{i} = cdf(theoreticalRician,xAxis);
    CDF_spectrum{i} = cdf(spectrumDistribution,xAxis);
    CDF_filter{i} = cdf(filterDistribution,xAxis);
end    
% Plot the PDF 
figure
subplot(3,1,1)
plot(xAxis,PDF_filter{1},xAxis,PDF_spectrum{1},xAxis,PDF_theory{1})
title(['PDF with K_c = ' num2str(K_c(1))])
legend('Filter method','Spectrum method','Theoretical')

subplot(3,1,2)
plot(xAxis,PDF_filter{2},xAxis,PDF_spectrum{2},xAxis,PDF_theory{2})
title(['PDF with K_c = ' num2str(K_c(2))])
legend('Filter method','Spectrum method','Theoretical')

subplot(3,1,3)
plot(xAxis,PDF_filter{3},xAxis,PDF_spectrum{3},xAxis,PDF_theory{3})
title(['PDF with K_c = ' num2str(K_c(3))])
legend('Filter method','Spectrum method','Theoretical')
% Plot CDF
figure()
subplot(3,1,1);
plot(xAxis,CDF_spectrum{1},xAxis,CDF_filter{1},xAxis,CDF_theory{1})
title(['CDF with K_c = ' num2str(K_c(1))] )
legend('Spectrum simulation','Filter simulation','Theoretical')

subplot(3,1,2);
plot(xAxis,CDF_spectrum{2},xAxis,CDF_filter{2},xAxis,CDF_theory{2})
title(['CDF with K_c = ' num2str(K_c(2))] )
legend('Spectrum simulation','Filter simulation','Theoretical')

subplot(3,1,3);
plot(xAxis,CDF_spectrum{3},xAxis,CDF_filter{3},xAxis,CDF_theory{3})
title(['CDF with K_c = ' num2str(K_c(3))] )
legend('Spectrum simulation','Filter simulation','Theoretical')


%% ACF
t = (-N_s*T_s:T_s:N_s*T_s);
figure(5)
for i = 1:3
    % calculate theoretical autocorr
    ACF_theory = K_c(i)^2 + besselj(0,2*pi*f_D*t);
    
    % calculate simulation autocorr
    [ACF_spectrum, lags_spectrum] = xcorr(spectrumMethod(f_D,T_s,N_s,K_c(i)),'unbiased');
    [ACF_filter, lags_filter] = xcorr(filterMethod(f_D,T_s,N_s,K_c(i)),'unbiased');
    subplot(3,1,i)
    plot(t,ACF_theory,lags_spectrum.*T_s,ACF_spectrum,lags_filter.*T_s,ACF_filter)
    grid on, legend('theoretical','spectrum method','filter method','Orientation','horizontal'), legend('boxoff')
    title(['Auto-correlation of the channels with K_c = ' num2str(K_c(i))])
    xlim([-0.2 1])
end
%% 2) Create Rician fading time and frequency varying channels
% given parameters
N = 300;
M = 64;
L = [1 2 3];
K_c = [0 1 10];
fDTs = [0.1 0.005];
% get fD from the normalized fDTs
Ts=T_s; % use the same T_s as before
fD = fDTs/Ts; % this is assuming fDTs = fD*Ts


% go through every 18 combinations to create the plots
for i = 1:3
    for j = 1:3
        for k = 1:2
            l = L(i);
            switch l % depending on the amount of channels L
                case 1
                    c_spectrum = spectrumMethod(fD(k),Ts,N,K_c(j));
                    c_filter = filterMethod(fD(k),Ts,N,K_c(j))';
                    c_paddedS = padWithZeros(c_spectrum,M,N);
                    c_paddedF= padWithZeros(c_filter,M,N);
                    C_S = fft(c_paddedS,M);
                    C_F = fft(c_paddedF,M);
                    figure()
                    subplot(2,2,1)
                    mesh(0:N-1,0:M-1, abs(C_S));
                    title(['mesh spectrum ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,2)
                    surf(1:N, 0:M-1,abs(C_S), 'MeshStyle', 'row'), view(0,90)
                    title(['surf spectrum ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,3)
                    mesh(0:N-1,0:M-1, abs(C_F));
                    title(['mesh fitler ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,4)
                    surf(1:N, 0:M-1,abs(C_F), 'MeshStyle', 'row'), view(0,90)
                    title(['surf filter ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                case 2
                    c_spectrum = spectrumMethod(fD(k),Ts,N,K_c(j));
                    c_spectrum(2,:) = spectrumMethod(fD(k),Ts,N,K_c(j));
                    c_filter = filterMethod(fD(k),Ts,N,K_c(j))';
                    c_filter(2,:) = filterMethod(fD(k),Ts,N,K_c(j))';
                    c_paddedS = padWithZeros(c_spectrum,M,N);
                    c_paddedF= padWithZeros(c_filter,M,N);
                    C_S = fft(c_paddedS,M);
                    C_F = fft(c_paddedF,M);
                    figure()
                    subplot(2,2,1)
                    mesh(0:N-1,0:M-1, abs(C_S));
                    title(['mesh spectrum ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,2)
                    surf(1:N, 0:M-1,abs(C_S), 'MeshStyle', 'row'), view(0,90)
                    title(['surf spectrum ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,3)
                    mesh(0:N-1,0:M-1, abs(C_F));
                    title(['mesh fitler ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,4)
                    surf(1:N, 0:M-1,abs(C_F), 'MeshStyle', 'row'), view(0,90)
                    title(['surf filter ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                case 3
                    c_spectrum = spectrumMethod(fD(k),Ts,N,K_c(j));
                    c_spectrum(2,:) = spectrumMethod(fD(k),Ts,N,K_c(j));
                    c_spectrum(3,:) = spectrumMethod(fD(k),Ts,N,K_c(j));
                    c_filter = filterMethod(fD(k),Ts,N,K_c(j))';
                    c_filter(2,:) = filterMethod(fD(k),Ts,N,K_c(j))';
                    c_filter(3,:) = filterMethod(fD(k),Ts,N,K_c(j))';
                    c_paddedS = padWithZeros(c_spectrum,M,N);
                    c_paddedF= padWithZeros(c_filter,M,N);
                    C_S = fft(c_paddedS,M);
                    C_F = fft(c_paddedF,M);
                    figure()
                    subplot(2,2,1)
                    mesh(0:N-1,0:M-1, abs(C_S))
                    title(['mesh spectrum ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,2)
                    surf(1:N, 0:M-1,abs(C_S), 'MeshStyle', 'row'), view(0,90)
                    title(['surf spectrum ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,3)
                    mesh(0:N-1,0:M-1, abs(C_F));
                    title(['mesh fitler ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
                    subplot(2,2,4)
                    surf(1:N, 0:M-1,abs(C_F), 'MeshStyle', 'row'), view(0,90)
                    title(['surf filter ' 'K_c=' num2str(K_c(j)) ' fDTs=' num2str(fDTs(k)) ' L=' num2str(L(i))])
            end
        end
    end
end



