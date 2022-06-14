function g = g_filter(N,Ts,fd)

n=-N:N;

g=zeros(1,(2*N+1));  % create filter

 if n == 0
        g(1,N+1) = (pi.*fd)^(1/4)/gamma(5/4);
    else
        g(1,N+1) = besselj(1/4,2*pi*fd*abs(n.*Ts)/(abs(n.*Ts).^ (1/4)));
 end
 

  g_norm = g/sqrt(mean(abs(g).^2)); % normalise g to unit energy
  
  g= g_norm;
 
end
