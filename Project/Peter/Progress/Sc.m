function Sc = Sc(f,fD)
if abs(f) <= fD
    Sc = (1/pi*fD)*(1/sqrt(1-(f/fD)^2));
else
    Sc = 0;
end
