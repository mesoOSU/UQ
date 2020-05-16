%Copyright (c) 2020, Denielle Ricciardi
%All rights reserved.                  

function [log_likelihood,prop_vpsc] = log_likelihood(curr_theta,curr_Lambda,curr_delta,curr_Delta,sigma_obs,strain_inc,sigma_vpsc,b,block_index)

S = size(sigma_obs,2);
D = 4;
N = length(sigma_obs);

%Error precision matrix
err_pr_mat = diag(repmat(curr_delta,1,N));

% Pre-allocate 
log_like_between = 0;
log_like_within = 0;

prop_vpsc = sigma_vpsc;
if b <= 6
	%Replace one column of prop_vpsc with VPSC evaluated at proposed parameters
    par = curr_theta(block_index{b});
    prop_vpsc(:,b) = VPSC(par(1),par(2),par(3),par(4),strain_inc);
end

Z1 = -D*log(2*pi) + log(det(curr_Lambda));
Z2 = -N*log(2*pi) + sum(log(diag(err_pr_mat)));

for i = 1:S    
    log_like_between = log_like_between + (1/2)*Z1 - (1/2)*(curr_theta(block_index{i}) - curr_theta(block_index{S+1}))...
        *curr_Lambda*(curr_theta(block_index{i}) - curr_theta(block_index{S+1}))';
    log_like_within = log_like_within + (1/2)*Z2-(1/2)*(sigma_obs(:,i) - prop_vpsc(:,i) - curr_Delta)'...
        *err_pr_mat*(sigma_obs(:,i)-prop_vpsc(:,i) - curr_Delta);
end

log_likelihood = log_like_between + log_like_within;

end