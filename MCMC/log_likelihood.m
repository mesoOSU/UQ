%This is the log likelihood function for the random effects model.

%yture is (150x6) 150 data points in each of six data sets
% par is N x (4 x 7)

function [log_likelihood,prop_vpsc] = log_likelihood(curr_theta,curr_Delta,curr_delta,sigma_obs,strain_inc,sigma_vpsc,block)

S = size(sigma_obs,2);
D = 4;
N = length(sigma_obs);

curr_re_tau0 = curr_theta(1:D:S*D);
curr_re_tau1 = curr_theta(2:D:S*D);
curr_re_theta0 = curr_theta(3:D:S*D);
curr_re_theta1 = curr_theta(4:D:S*D);
curr_inv_delta = 1/curr_delta;
%par(25:28 hold the super parameters)

%covariance matrix for the observation error, including model discripancy
%for the initial 'fix' points

%err_cov_mat = diag(repmat(curr_inv_delta,1,T)); 
err_pr_mat = diag(repmat(curr_delta,1,N));

% Pre-allocate 
log_like_between = 0;
log_like_within = 0;

prop_vpsc = sigma_vpsc;
if block <= 6
prop_vpsc(:,block) = VPSC(curr_re_tau0(block),curr_re_tau1(block),curr_re_theta0(block),curr_re_theta1(block),strain_inc);
%replace one column of prop_vpsc with VPSC evaluated at proposed parameters
end

    Z1 = -D*log(2*pi) + log(det(curr_Delta));

    Z2 = -N*log(2*pi) + sum(log(diag(err_pr_mat)));

for i = 1:S
    
    log_like_between = log_like_between + (1/2)*Z1 - (1/2)*(curr_theta(D*i-(D-1):D*i) - curr_theta(S*D+1:(S+1)*D))...
        *curr_Delta*(curr_theta(D*i-(D-1):D*i) - curr_theta(S*D + 1: (S+1)*D))';
    log_like_within = log_like_within + (1/2)*Z2-(1/2)*(sigma_obs(:,i) - prop_vpsc(:,i))'...
        *err_pr_mat*(sigma_obs(:,i)-prop_vpsc(:,i)); 
end

log_likelihood = log_like_between + log_like_within;
end