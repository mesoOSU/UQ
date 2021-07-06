function [f,gp_cov] = joint_pdf(y_proj,c_design,pars,theta,b_inv,c,sigma,gpsigma,tau,phi,method,v)

    % Log joint probability density of training and observed responses
    % y_proj : observation projected to spline coefficients
    % c_design : training 
    coeff=[c,c_design,y_proj];
    if nargin < 11
        method = "matern";
    elseif nargin > 11
        cov_thetas = cov_matern_v([theta;pars;theta],gpsigma,tau,phi,v);
    end
    
    if method == "matern"
        cov_thetas = cov_matern([theta;pars;theta],gpsigma,tau,phi);
    elseif method == "matern_phi"
        cov_thetas = cov_matern2([theta;pars;theta],gpsigma,tau,phi);
    elseif method == "exponential"
        cov_thetas = cov_exp([theta;pars;theta],gpsigma,tau,phi);
    elseif method == "squared_exponential"
        cov_thetas = cov_sqexp([theta;pars;theta],gpsigma,tau,phi);
    end
    
    gp_cov = zeros(size(c_design,1)*(size(c_design,2)+2));
    a = size(c_design,1);
    b = size(c_design,2);
    for i = 1:size(c_design,2)+2
        for j = 1:size(c_design,2)+2
            gp_cov((a*i-a+1):a*i,(a*j-a+1):a*j) = diag(repmat(cov_thetas(i,j),size(c_design,1),1));
        end
    end
    
    gp_cov(end-a+1:end,end-a+1:end) = gp_cov(end-a+1:end,end-a+1:end) + sigma^2*b_inv;
    
    f = logmvnpdf(reshape(coeff,(b+2)*a,1)',zeros(1,(b+2)*a),gp_cov);
end