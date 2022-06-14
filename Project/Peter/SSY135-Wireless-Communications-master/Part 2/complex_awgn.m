function output = complex_awgn(input, EbN0, rate, bits_per_symbol)
%This function file is to add the AWGN noise to the signal.
% input: The input signal assumed to have unit power.
% EbN0 - The target SNR in power per source bit. Measured in dB.
% rate - The overall code rate for the designed system
% QPSK: bits_per_symbol is 2
W_EbN0 = 10 ^ (EbN0 / 10); %Bit Energy conversion to Watt
W_EsN0 = W_EbN0 * bits_per_symbol  * rate;
N0 = sqrt(1 / (W_EsN0 * 2));
%Create the random noise component, modelled as complex Gaussian
awgn = N0 * (randn(length(input), 1) + 1i * randn(length(input), 1));
output = input + awgn;
end