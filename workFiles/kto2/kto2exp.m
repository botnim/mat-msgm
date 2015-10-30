function [U_, E_, P_, const] = kto2exp(U,E,P,gP)
%Kto2exp transform a k-labeled pairwise MRF to binary (pairwise) MRF
%
% we call this binary representation 'explicit' because in contrary to
% Schlesinger's reprarameterization, in this representation each K-labeled
% node is represented by (K-1) nodes such that (i_k = 1) iff (i's lbl is k)
%
% the price to pay for explicit representation is that each node transforms
% to (K-1)-sized clique, instead of K-1 chain
%

CONST = gP.CONST;        % replaces INF

N = size(U,1);  % number of vars
M = size(E,1);  % number of edges
K = size(U,2);  % number of labels


% reparameterize the graph such that U(i,1) = 0 for all i,
% and P(li,1,k) = 0, P(1,lj,k) = 0 for all li, lj, k
ui = P(:,1,:);
P = bsxfun(@minus,P,ui);
uj = P(1,:,:);
P = bsxfun(@minus,P,uj);
ui = squeeze(ui)';
uj = squeeze(uj)';
for i = 1 : K
    
    U(:,i) = U(:,i) + accumarray(E(:,1),ui(:,i),[N,1]) + ...
        accumarray(E(:,2),uj(:,i),[N,1]); 
end
const = sum(U(:,1));
U = bsxfun(@minus,U,U(:,1));
U = U(:,2:end)';

%
% build proxy graph

% unary term
U_ = zeros((K-1)*N,2);
U_(:,2) = U(:);

% pairwise (cliques)
E_clique = nchoosek(1:(K-1),2);     % all connections within a clique
E_clq_src = bsxfun(@plus,(0:(N-1))' * (K-1),E_clique(:,1)');
E_clq_trgt = bsxfun(@plus,(0:(N-1))' * (K-1),E_clique(:,2)');
P_clq = zeros(2,2,numel(E_clq_src));
P_clq(2,2,:) = CONST;

% "real" edges (re)
mAlpha = P(2:K,2:K,:);
v_re_src_idx(1,1,:) = (E(:,1)-1) * (K-1);
v_re_trgt_idx(1,1,:) = (E(:,2)-1) * (K-1);
E_re_src = bsxfun(@plus,v_re_src_idx,repmat((1:K-1)',[1,(K-1)]));
E_re_trgt = bsxfun(@plus,v_re_trgt_idx,repmat((1:K-1),[(K-1), 1]));
P_re = zeros(2,2,((K-1)^2)*M);
P_re(2,2,:) = mAlpha(:);

P_ = cat(3,P_clq,P_re);
E_ = [E_clq_src(:), E_clq_trgt(:); E_re_src(:), E_re_trgt(:)];





