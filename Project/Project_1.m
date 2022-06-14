clear all, close all 

% Parameter definition
par.fc = 2e9; % Carrier [GHz]
par.lambda = 3e8/par.fc; % Wavelength [m]
par.v = 30e3/3600; % Relative velocity [m/s]
par.Ts = 0.1e-3; % Sampling rate [s]
% par.BW = 1/par.Ts; % Bandwidth [Hz]
% par.Rs = par.BW; % Symbol rate [S/s]
par.fD = par.v/par.lambda; % Doppler spread [Hz], no aliasing; BW=10k > 2fD=110
par.L = 1; % Number of taps 
par.channel = 'Rician'; % Choose type {'Rayleigh', 'Rician'}
par.method = 'Spectrum'; % Choose method {'Filter', 'Spectrum'}
par.kc = 1; % Rician factor for dominating path
par.syms = 10; % Number of symbols transmitted
par.Ns = 3000; %2^12; % Number of samples for the signal
par.N = par.Ns*par.syms; % Number of samples for the filter
par.fs = 1/par.Ts; % Sampling frequency
par.DS = par.Ts*2; % Delay spread
par.M = 64; % Zero-padding (oversampling in delay domain?)

% Preallocation of c
c = NaN(par.L,par.Ns);

switch par.method
    
    case 'Filter'
        nTs = -par.N*par.Ts:par.Ts:par.N*par.Ts; % Time vector for sampling

        % Define filter g(t)
        gt = besselj(1/4,2*pi*par.fD.*abs(nTs))./nthroot(abs(nTs),4); % For t~=0
        gt(par.N+1) = nthroot((pi*par.fD),4)/gamma(5/4); % For t=0
        gt = gt/sqrt(mean(abs(gt).^2)); % Normalise to unit energy
        
        E_Eg = mean(abs(gt).^2); % Expected energy of g(t)

        % Generate Ns random complex distribution for the channel
        x = sqrt(1/2).*(randn(par.L,par.Ns)+1j*randn(par.L,par.Ns)); % Define complex random Gaussian distribution
        var_x = var(x,[],2); % Unit variance
        figure, histogram(abs(x)) % Plot the distribution
        title('Distribution of x')

        % Generate channel c(t)
        switch par.channel
            case 'Rayleigh'
                for l = 1:par.L
                    c(l,:) = conv(x(l,:),gt,'same');
                    c(l,:) = c(l,:)./sqrt(mean(abs(c(l,:)).^2)); % Normalize to unit energy
                end
            case 'Rician'
                for l = 1:par.L
                    c(l,:) = conv(x(l,:),gt,'same') + par.kc;
                    c(l,:) = c(l,:)./sqrt(mean(abs(c(l,:)).^2)); % Renormalize
                end
        end

        % Calculate expected values of c
        E_Ec = mean(abs(c).^2,2); % E(energy) = 1
        E_c = mean(real(c),2); % E(c) = 0
        
        % Plot g(t) in time
        figure, plot(nTs,gt);
        xlabel('Time [s]')
        ylabel('g(t)')
        title('Time domain filter')

        % Plot the distribution of c(t)
        figure, histogram(abs(c(1,:)),100) % Rayleigh distribution in linear (Log-normal) (Gaussian in dB)
        title('Distribution of c(t) 1st tap')

    case 'Spectrum'

        f = -par.fD:par.fs/par.Ns:par.fD; % Defining the frequency interval
        nTs = -par.Ns*par.Ts:par.Ts:par.Ns*par.Ts; % Time vector for sampling

        Sf = 1./(pi*par.fD.*sqrt(1-(f./par.fD).^2));
        Sf(1) = Sf(end); % Get rid of the 'Inf' at -fD
        
        Gf = sqrt(Sf);
        
        % The periodic extension of G for k = 0:1:par.Ns-1
        Gfp = [Gf(length(Gf)/2+1:end) zeros(1,par.Ns-length(Gf)) Gf(1:length(Gf)/2)]; 
        
        % PSD of Gf
        Sc_f = abs(Gf).^2; 
        figure(1),clf, plot(f,real(Gf))
        xlabel('Frequency [Hz]')
        ylabel('G(f)')
        title('Theoretical PSD of G(f)')

        % Random complex distribution
        X = sqrt(1/2)*(randn(par.L,par.Ns)+1j*randn(par.L,par.Ns)); % Random complex 
        var_X = var(X,[],2); % Unit variance

        C = Gfp.*X;
        c = ifft(C,par.Ns,2);
        
        switch par.channel
            case 'Rayleigh'
                for l = 1:par.L
                    c(l,:) = c(l,:)./sqrt(mean(abs(c(l,:)).^2,2)); % Normalize to unit energy
                end
            case 'Rician'
                for l = 1:par.L
                    c(l,:) = c(l,:) + par.kc;
                    c(l,:) = c(l,:)./sqrt(mean(abs(c(l,:)).^2,2)); % Renormalize
                end
        end
        
        figure, histogram(abs(c(1,:)));
        E_Ec = mean(abs(c).^2,2);
        
        x = 0:0.1:10;
        dist = fitdist(abs(c).','Rician');
        PDF = pdf(dist,x);
        CDF = cdf(dist,x);
        
        Ac = besselj(0,2*pi*par.fD.*nTs);
        Ac_est = autocorr(c.',500);
        
        if par.L > 1
            c = [c; zeros(par.M-par.L,par.Ns)]; % Zero padding
            Cf = fft(c,par.M,1);
        end
end

% Theoretical Doppler PSD
fD_vec = linspace(-par.fD,par.fD,100); % Doppler spectrum vector
Sc = 1./(pi*par.fD.*sqrt(1-(fD_vec./par.fD).^2));

% Plot theoretical and estimated Doppler spectrum
figure
subplot(1,2,1)
plot(fD_vec./par.fD,Sc)
title('Theoretical Doppler spectrum')
xlabel('f/f_D'), ylabel('S_c(f)')
subplot(1,2,2)
pwelch(c(1,:),[],[],[],par.fs) % Estimate of PSD
title('Estimated Doppler spectrum')

figure, plot(nTs,Ac)
title('Theoretical ACF')

figure
subplot(1,2,1)
plot(PDF)
title(['PDF k = ' num2str(par.kc)])
subplot(1,2,2)
plot(CDF)
title(['CDF k = ' num2str(par.kc)])


        % Plot the channel c(t)
%         figure(3),clf, plot3(par.L,nTs,real(c(par.L)))
%         xlabel('Time [s]')
%         ylabel('Re(c(t))')


