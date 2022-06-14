function Gp = Gp(f,fD,fs)
if (0 <= f) && (f < fD)
    Gp = sqrt((1/(pi*fD))*(1./sqrt(1-(f/fD).^2)));
elseif  ((fs-fD) < f) && (f < fs)
    Gp = sqrt((1/(pi*fD))*(1./sqrt(1-((f-fs)/fD).^2)));
else
    Gp = 0;
end
