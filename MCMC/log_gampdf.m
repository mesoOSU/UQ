function log_gampdf = log_gampdf(x,shape,rate)
%un-normalized gamma pdf
%a_tau = shape, b_tau = rate

log_gampdf = (shape - 1)*log(x) - (x*rate);

end