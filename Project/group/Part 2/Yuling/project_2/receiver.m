function [s_est_sym, s_est_bits] = receiver(r_cp, h, const, N, Ncp, i)

r_rx = r_cp(Ncp+1 : end); % remove cp
y = fft(r_rx,N); 

c = fft(h(i,:),N);
%c_l = c(i);

s_rx = y .* c';

d = abs(repmat(s_rx,1,1) - repmat(const, length(s_rx), 1)).^2;
[~, index] = min(d, [], 2); % pick out the bits with shortest distances to the symbols
s_est_sym = const(index');

index = index - 1;
s_binary = de2bi(index, 2, 'left-msb');  % 2: since QPSK
s_bi = s_binary';
s_est_bits = s_bi(:)';


end