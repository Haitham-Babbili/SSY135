fc=2e9; %in Hz
Ts = 0.1e-3;
v = 30; %in km/h
Ns = 1e4;
N = 100;
fD = (v/(3.6*3e8))*fc;
max_Ts = 1/(2*fD);
c = channelByFilter(Ts,Ns,N,fD);
pwelch(c);
sc = zeros(1,2*ceil(fD));
for i = -ceil(fD):ceil(fD)
    sc(1,i+ceil(fD)+1) = Sc(i,fD);
end
figure(2);
plot(sc);