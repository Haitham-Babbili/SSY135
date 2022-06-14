function [ c ] = rician_channel(Ts,fD,kc,L_taps,Mf,Ns,method )
%
%    Ts: sampling time
%    kc: Line of Sight (LOS)dominant component 
%    tc: Time shift, tau
%    fD: Doppler frequency in Hertz
%   L_taps: number of taps
%   Mf: filter length
%   Ns: number of samples, related to frequency resolution??
%   Approach methods: 0 for filtering, 1 for spectrum method
c = zeros(L_taps *2 ,Ns);
switch method
    case 0
        for tc = 1:2:L_taps*2
            method = channelByFilter(Ns, Ts, fD, kc,Mf);
            c(tc,:) = method;
        end    
    case 1
        for tc = 1:2:L_taps*2
            method = channelBySpectrum(Ns, Ts, fD, kc);
            c(tc,:) = method;
        end        
    otherwise
        error('invalid method')
end

end

