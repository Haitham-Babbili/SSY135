function c_padded = padWithZeros(c,M,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[L, ~] = size(c);
padding = zeros(M-L,N);
c_padded = [c; padding];
end

