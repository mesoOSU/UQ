function [rnd,loglik] = full_cond(c_design,y_proj,gp_cov)
    a = size(c_design,1);
    b = size(c_design,2);
    inv_cov = inv(gp_cov(a+1:end,a+1:end));
    mu = gp_cov(1:a,a+1:end)*inv_cov*[reshape(c_design,a*b,1);y_proj];
    full_cov = gp_cov(1:a,1:a)-gp_cov(1:a,a+1:end)*inv_cov*gp_cov(a+1:end,1:a);
    full_cov = (full_cov + full_cov')/2;
    rnd = mvnrnd(mu,full_cov)';
    loglik = logmvnpdf([rnd',reshape(c_design,a*b,1)',y_proj'],zeros(1,a*(b+2)),gp_cov);
end