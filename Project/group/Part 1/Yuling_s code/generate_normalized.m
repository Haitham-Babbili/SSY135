function c_norm = generate_normalized(c)

%c_eng = mean(abs(c).^2);
c_norm = c/sqrt(mean(abs(c).^2));

end