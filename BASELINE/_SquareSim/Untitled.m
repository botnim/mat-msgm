mu = mean(mean(P,1),2);
P_ = bsxfun(@minus,P,mu);
P_ = reshape(P_,size(P,1) * size(P,2), size(P,3));
sigma = std(P_,0,2);
MU = mean(sigma)


mu = mean(U,2);
U_ = bsxfun(@minus,U,mu);
sigma = std(U_,0,2);
MU = mean(sigma)