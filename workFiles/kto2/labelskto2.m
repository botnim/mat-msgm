function L_ = labelskto2(L,N,K)

L_ = repmat(1:K,N,1);
L_ = bsxfun(@ge,L_,L);
L_ = L_(:,1:end-1)';
L_ = L_(:) + 1;