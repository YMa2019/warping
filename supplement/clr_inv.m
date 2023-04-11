function gamma = clr_inv(f, T)
    temp = exp(f)/(trapz(T,exp(f)));
    gamma = cumsum(temp)/sum(temp); 
    gamma=(gamma-min(gamma))/(max(gamma)-min(gamma));
end