function [vEdgeList, Hji] = ScoreEdges_short(U,E,P,prior,L,gP)

%
% SelectEdges(U,E,P,L,gP) - create an ordered list of edges to be processed
% by CompCoarseGraph(..)
%
% OUTPUT:
%
% vEdgeList -   an ordered list of edges, Mx1 (or 2Mx1) array
%               vEdgeList(k) holds the index (wrt E,P) of the k-th best
%               edge.
%                   Denote: i := E(k,1)
%                           j := E(k,2)
%                   if vEdgeList(k)>0 then edge interpolates i -> j
%                   if vEdgeList(k)<0 then edge interpolates j -> i
%

if gP.wDegEnt
    % the unary term of each vertex is normalized by
    % the vertex' degree, compute the degree for all the vertices
      
    N = size(U,1);                      % number of vertices
    vDeg = accumarray(E(:,1),1,[N,1]);  % number of edges for which the
                                        % vertex is on the left side
                                        
    vDeg = vDeg + accumarray(E(:,2),1,[N,1]);   % degree of each vertex
end

%
% sum all energy terms into P,
% such that P(li,lj,k) = phi_i(li) + phi_j(lj) + phi_ij(li,lj)
if gP.wDegEnt
    % normalize the unary term by the degree
    % s.t. P(li,lj,k) = phi_i(li)/deg(i) + phi_j(lj)/deg(j) + phi_ij(li,lj)
    U = bsxfun(@rdivide,U,vDeg);
end


%
% calculate conditional entropies
[Hji, Hij] = CalcCondEnt(U,E,P,prior);

indRev = [zeros(size(Hji)); ones(size(Hij))];
Hji = [Hji; Hij];
[Hji, ind] = sort(Hji,'ascend');

M = size(E,1);
vEdgeList = ind;
vEdgeList(vEdgeList > M) = - (vEdgeList(vEdgeList > M) - M);



function [Hji, Hij] = CalcCondEnt(U,E,P,prior)


M = size(P,3);
K = size(P,1);

Hji = zeros(M,1);
Hij = zeros(M,1);

%
% loop over edges, calculate entropy in both directions
for iM = 1 : M
    
    ii = E(iM,1);   % left vertex of edge
    jj = E(iM,2);   % right vertex of edge
    
    %
    % double-loop over the pairwise term,
    % to add the unary terms and to exponentiate
    % ..online update of conditional entropy
    % ..and online updates of the marginals
    Hji_ = zeros(1,K);      % cond. ent. of jj given ii
    Hij_ = zeros(1,K);      % cond. ent. of ii given jj
    pi = zeros(1,K);        % unnormalized marginal of ii
    pj = zeros(1,K);        % unnormalized marginal of jj
    mi = zeros(1,K);        % (normalized) marginal of ii
    mj = zeros(1,K);        % (normalized) marginal of jj
    for i = 1 : K
        
        if any(prior)
            mi(i) = prior(ii,i);
            mj(i) = prior(jj,i);
        end
        
        for j = 1 : K
            
            % add unary terms
            P(i,j,iM) = P(i,j,iM) + U(ii,i) + U(jj,j);
            
            % exp.
            P(i,j,iM) = exp(-P(i,j,iM));
            
            % update marginals
            pi(i) = pi(i) + P(i,j,iM);
            pj(j) = pj(j) + P(i,j,iM);
            
            % update entropy
            Hji_(i) = Hji_(i) + (P(i,j,iM) * log(P(i,j,iM)));
            Hij_(j) = Hij_(j) + (P(i,j,iM) * log(P(i,j,iM)));
            
        end
    end
    
    
    %
    % prepare the marginal of ii & jj,
    % when there's no prior
    if ~any(prior)
        % need to normalize pi,pj to proper distribution
        si = sum(pi);
        sj = sum(pj);
        mi = pi / si;
        mj = pj / sj;
    end
        
    
    %
    % calculate the conditional entropy
    Hjitmp = 0;
    Hijtmp = 0;
    for i = 1 : K
        
        Hjitmp = Hjitmp + mi(i) * ( - Hji_(i) / pi(i) + log(pi(i)));
        Hijtmp = Hijtmp + mj(i) * ( - Hij_(i) / pj(i) + log(pj(i)));
        
    end
    
    Hji(iM) = Hjitmp;
    Hij(iM) = Hijtmp;
    
end