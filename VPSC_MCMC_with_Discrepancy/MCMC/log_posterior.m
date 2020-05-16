%Copyright (c) 2020, Denielle Ricciardi
%All rights reserved.                  

function log_posterior = log_posterior(S,N,D,block_index,log_like,curr_theta,curr_t2,curr_R,curr_delta,curr_Delta,log_det_Gamma,inv_Gamma,a_delta,b_delta,R_not,r_not,a_t,b_t)

%Hyperparameters for theta (Gamma & Normal)
hyp1 = [.98,.98,1800,600];
hyp2 = [.014,.014,1E6,1E6];

overall_log_priors = {@log_gampdf,@log_gampdf,@log_normal,@log_normal};

overall_par = curr_theta(block_index{S+1});
overall_log_prior = 0;

for dd = 1:D
    overall_log_prior = overall_log_prior + overall_log_priors{dd}(overall_par(dd),hyp1(dd),hyp2(dd));
end


Z = -N*log(2*pi) - log_det_Gamma;

Delta_log_prior = (1/2)*Z - (1/2)*curr_Delta'*inv_Gamma*curr_Delta;

t2_log_prior = 0;
for d = 1:D
    t2_log_prior = t2_log_prior + log_gampdf(curr_t2(d),a_t,b_t);
end

R_log_prior = log_wishart_pdf(R_not,r_not,curr_R,D);

delta_log_prior = log_gampdf(curr_delta(end),a_delta,b_delta);

log_posterior = double(log_like + overall_log_prior + Delta_log_prior + t2_log_prior + R_log_prior + delta_log_prior);
end