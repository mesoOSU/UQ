function C = cov_sqexp(z,sigma2,tau2,phi)
    C = zeros(size(z,1));
for i = 1:size(z,2)
    D=pdist(z(:,i));
    t=squareform(D);
    C = C-(t.^2)./phi(i);
end
C = sigma2*exp(C);
C(C==sigma2) = sigma2+tau2;
end
