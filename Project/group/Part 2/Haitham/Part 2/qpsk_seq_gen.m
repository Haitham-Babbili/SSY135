function [Z]= qpsk_seq_gen(N,E,T_s,Ncp) 
%N: lenght of the QPSK sequence
%E: energie per symbol
data = randi([0 1],2*N,1,'double');       %random bitstream of lengt 2*N
symbols = (sqrt(E))*nrSymbolModulate(data,'QPSK','OutputDataType','double');        %QPSK sequence of length N

z = sqrt(N/T_s) * ifft(symbols);
%add cp
Z = [z(end-Ncp+1:end);z];

end

