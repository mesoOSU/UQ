function qqt = sq_exp(u1,v1,w,sigma)

if (size(u1,1)<size(u1,2))
    u1 = u1';
end

if (size(v1,1)<size(v1,2))
    v1 = v1';
end

u = repmat(u1,1,length(v1));
v = repmat(v1',length(u1),1);


qqt = sigma.*exp(-(u-v).^2./(2*w^2));

end