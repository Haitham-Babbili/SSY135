function [symbols, data]= qpsk_sequence(N,E) 
%N: lenght of the QPSK sequence
%E: energie per symbol
data = randi([0 1],2*N,1,'double');       %random bitstream of lengt 2*N
symbols = (sqrt(E))*nrSymbolModulate(data,'QPSK','OutputDataType','double');        %QPSK sequence of length N
end