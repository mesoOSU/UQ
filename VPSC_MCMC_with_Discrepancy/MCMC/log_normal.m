function log_norm = log_normal(x,mu,variance)
%log normal up to a constant of proportionality
    log_norm = -((x-mu).^2)./(2.*variance);
end