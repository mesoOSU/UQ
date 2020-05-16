function log_wishpdf = log_wishart_pdf(V_not,v_not,Delta,D)
%Delta is a DxD inverse variance-covariance matrix
%This is the log of the un-normalized wishard pdf

log_wishpdf = (-v_not/2)*log(det(V_not))+((v_not-D-1)/2)*log(det(Delta))-(1/2)*trace(V_not\eye(D)*Delta);
end