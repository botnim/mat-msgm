function eng = Energy_gen(D,E,W,vL)

%
% computes energy
%

vL = vL(:);

%
% pairwise term
[ii, jj, ww] = find(triu(E));
eng = sum(W(sub2ind(size(W),ww,vL(ii),vL(jj))));

%
% unary term
eng = eng + sum(D(sub2ind(size(D),1:size(D,1),vL')));
