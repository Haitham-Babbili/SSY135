function C_padding = channel_response(c,L,M,N)

for m = 1:L
    for n = 1:N
        c_mtx(m,n) = c{m}(n);
    end
end
c_padding = [c_mtx; zeros(M-L,N)];
C_padding = fft(c_padding);


end