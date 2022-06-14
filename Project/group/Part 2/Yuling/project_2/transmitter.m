function z_tx = transmitter(bits, N, Ncp)

% bits_de = bi2de(bits,'left-msb')+1;     % Make binary values into decimal values
% bits = reshape(bits(:), 2, N); % in QPSK, 2bits/symbol
% bits_de = bi2de(bits')+1; % matlab index starts from 1
% symbols = const(bits_de).';   % symbols on QPSK
ak = 1 - 2*bits(1:2:N);  % Real part
bk = 1 - 2*bits(2:2:N);  % Imaginary part
symbols = ak + 1j*bk;


z_tx = sqrt(N) * ifft(symbols);
z_cp = z_tx(length(z_tx)-Ncp+1 : end); % cyclic prefix
z_tx = [z_cp z_tx]; % add cp to z


end