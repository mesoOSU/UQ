function log_posterior = log_posterior(log_like,curr_theta,curr_Delta,curr_delta,a_tau,b_tau,V_not,v_not,D)


log_posterior = log_like + log_gampdf(curr_theta(25),.98,.014) + log_gampdf(curr_theta(26),.98,.014)...
        + log_normal(curr_theta(27),1800,1E6) + log_normal(curr_theta(28),600,1E6)...
        + log_wishart_pdf(V_not,v_not,curr_Delta,D) + log_gampdf(curr_delta(end),a_tau,b_tau);

end