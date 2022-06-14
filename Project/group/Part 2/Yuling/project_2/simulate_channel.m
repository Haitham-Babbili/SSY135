function Y = simulate_channel(rx, L, N0, Ncp)
% this function is used to generate the siganl at receiver side and without
% cyclic prefix，finally, y is in frequency domain


for i = 1: L
    awgn_samples = awgn(rx, N0); 
    y{i} = awgn_samples(Ncp+1: end); % remove cp
    Y = fft(y{i});
    
end


end