function [U_, E_, P_, const] = kto2(U,E,P,gP)
%K2to transform a k-labeled pairwise MRF to boolean (pairwise) MRF
%
% algorithm from "Transforming an arbitrary minsum problem into a binary
% one", Dmitrij Schlesinger and Boris Flach, see pg. 11
%

CONST = gP.CONST;        % replaces INF

N = size(U,1);  % number of vars
M = size(E,1);  % number of edges
K = size(U,2);  % number of labels

N_ = N*(K-1);   % new number of vars
% M_ = N*(K-2) + ((K-1)^2)*M; % new number of edges


%
% (*) construct unary
% unary term of a variable which is indexed by 'i' is kept at U(i,:); the
% i-th variable is mapped to (K-1) variables indexed by
% (i-1)*(K-1) + 1 : i*(K-1)
% 

% from unary term
vdUfU = -1 * diff(U,1,2)';
vdUfU = vdUfU(:);

% from pairwise
mPreg = squeeze(P(:,K,:));
mPrev = squeeze(P(K,:,:));
mPreg = -1 * diff(mPreg);
mPrev = -1 * diff(mPrev);
vidxReg = bsxfun(@plus,(E(:,1)-1)*(K-1),1:(K-1))'; 	% see (*)
vidxRev = bsxfun(@plus,(E(:,2)-1)*(K-1),1:(K-1))';  % see (*)
vdPfP = accumarray([vidxReg(:); N_], [mPreg(:); 0]) + ...
    accumarray([vidxRev(:); N_], [mPrev(:); 0]);

% total
U_ = zeros(N*(K-1),2);
U_(:,2) = vdUfU + vdPfP;




%
% (**) construct pairwise

% "self loops" (sl)
E_sl_src = bsxfun(@plus,((1:N)'-1)*(K-1),1:(K-2));
E_sl_trgt = bsxfun(@plus,((1:N)'-1)*(K-1),2:(K-1));
P_sl = zeros(2,2,N*(K-2));
P_sl(2,1,:) = CONST;       	% heavy penalty on illegal labeling

% "real" edges (re)
mAlpha = P(1:(K-1),1:(K-1),:) + P(2:K,2:K,:) - ...
    P(2:K,1:(K-1),:) - P(1:(K-1),2:K,:);
v_re_src_idx(1,1,:) = (E(:,1)-1) * (K-1);
v_re_trgt_idx(1,1,:) = (E(:,2)-1) * (K-1);
E_re_src = bsxfun(@plus,v_re_src_idx,repmat((1:K-1)',[1,(K-1)]));
E_re_trgt = bsxfun(@plus,v_re_trgt_idx,repmat((1:K-1),[(K-1), 1]));
P_re = zeros(2,2,((K-1)^2)*M);
P_re(2,2,:) = mAlpha(:);

% new pairwise
E_ = [E_sl_src(:), E_sl_trgt(:); ...
    E_re_src(:), E_re_trgt(:)];
P_ = cat(3,P_sl,P_re);



%
% calculate const (for debugging)
const = sum(U(:,K)) + sum(P(K,K,:));


end

