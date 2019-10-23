function cov_mat_new = full_cov_mat_adj(par,D,cov_mat_orig)

    %adjust the ratio of the proposal variances based on the covariance of
    %the parameters
    def_by = 1;                    

    chaincorr = corr(par);
    chaincorr(isnan(chaincorr)) = 0;                    

    new_stepcor_chain = zeros(D,D);
    temp_stepvar_chain = zeros(D,D);

    for rr = 1:D
        for cc = 1:D
            if (rr==cc)
                new_stepcor_chain(rr,cc) = 1;
                elseif (abs(chaincorr(rr,cc))<0.6)
                    new_stepcor_chain(rr,cc)=0;
                else
                     new_stepcor_chain(rr,cc) = sign(chaincorr(rr,cc))*(abs(chaincorr(rr,cc))-0.4);
            end
        end
    end

    pvars = diag(cov_mat_orig)/def_by;

    for rr = 1:D
        for cc = 1:D
            temp_stepvar_chain(rr,cc) = new_stepcor_chain(rr,cc)*sqrt(pvars(rr)*pvars(cc));
        end
    end

    tempvar = diag(temp_stepvar_chain);

    for dd = 1:D
        if tempvar(dd)<=0
            temp_stepvar_chain(dd,dd) = cov_mat_orig(dd,dd);
        end
    end

    cov_mat_new = temp_stepvar_chain;
    
%**% In your case, I don't forsee needing these lines
%     [~,posdef] = chol(temp_stepvar_chain);
% 
%     if posdef ~= 0 || ~isnan(sum(sum(temp_stepvar_chain(logical(eye(D))))))
%         cov_mat_adj(logical(eye(D))) = temp_stepvar_chain(logical(eye(D)));
%     end 
end

                
     
