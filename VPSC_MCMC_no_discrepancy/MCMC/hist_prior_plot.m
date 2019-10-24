%plot priors with marginal historgram

subplot(2,2,1)
h = histogram(theta_chain(iter_start:k,1));
hold on
xgrid = min(theta_chain(iter_start:k,1))*.75:.01:max(theta_chain(iter_start:k,1))*1.25;
yval = exp(super_tau0_logprior(xgrid));
norm_y = yval/max(yval);
plot(xgrid,yval/norm_y * max(h.Values))
hold off
