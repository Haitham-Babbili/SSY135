clear;
close all
% Section D Trials
M = 300;                   %given time samples;
N = 64;                    % given frequency samples
L = [1 2 3];               %given Number of L_taps
fDTs = [.1 .005];          %relative doppler shift
kc = [0 1 10];             % Power of LOS components
Ts = 0.01;                 % sampling interval =>fD = fDTs
Mf = 1000;                 %The filter length should be about 2*N+1
Ns = 100000;
method = 1;                %0 = filter method; 1 = spectrum method
%a trial
% c = rician_generator(Ts,fDTs(2)/Ts,kc,3,400,10000,1);
% C = fft(c,N);
% plot_rician(C(1:N,1:M),M,N)
%
for k = 1:3
    %kc loop
    for j = 1:2
        %fdTs loop
        for k = 1:3
            %L loop
            L_taps = L(k);
            % c(tc,t) (time domain)
            c = rician_gen(Ts,fDTs(j)/Ts,kc(k),L_taps,Mf,Ns,method);
            %FFT Calculation with respect to tc(tau)
            C = fft(c,N);
            figure
            fig = plot_rician(C(1:N,1:M),M,N);
            title(sprintf('k = %d, L_taps = %d, fDTs = %.3f',kc(k),L_taps,fDTs(j)))
            saveas(fig,sprintf('taskB_(%d,%d,%.3f)_(k,L_taps,fDTs).png',kc(k),L_taps,fDTs(j)),'png')
        end
    end
end