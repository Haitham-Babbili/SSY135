function Gp = generate_filter(fd,fs)

f_Sc = (0 : fd);
Sc = (1/(pi*fd)).*(1./sqrt(1-(f_Sc./fd).^2)); 
Sc(1) = Sc(end);
G = sqrt(Sc);
Gp = zeros(1,fs);
Gp(1:length(G)) = G;
Gp(fs-length(G)+1:end) = sqrt(flip(Sc));

end