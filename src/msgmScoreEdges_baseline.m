function vEdgeList = msgmScoreEdges_baseline(U,E,P,prior,L,gP, bBin)

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

%%% TMP, TWEAK TO AVOID PASSING PARAMS
if isempty(gP)
    gP = msgmParams;
end
L = [];
%%%%%%%%%%%%%%%%

rng(1);

    % the unary term of each vertex is normalized by
    % the vertex' degree, compute the degree for all the vertices
      
    N = size(U,1);                      % number of vertices
    vDeg = accumarray(E(:,1),1,[N,1]);  % number of edges for which the
                                        % vertex is on the left side
                                        
    vDeg = vDeg + accumarray(E(:,2),1,[N,1]);   % degree of each vertex


%
% sum all energy terms into P,
% such that P(li,lj,k) = phi_i(li) + phi_j(lj) + phi_ij(li,lj)

    % normalize the unary term by the degree
    % s.t. P(li,lj,k) = phi_i(li)/deg(i) + phi_j(lj)/deg(j) + phi_ij(li,lj)
    U = bsxfun(@rdivide,U,vDeg);

ui(:,1,:) = U(E(:,1),:)';           % unary term of vertex i
uj(1,:,:) = U(E(:,2),:)';           % unary term of vertex j
P = bsxfun(@plus,P,ui);
P = bsxfun(@plus,P,uj);

%
% shift to probability space
P = bsxfun(@minus,P,mean(mean(P,1),2));
P = exp(-P);
P = bsxfun(@rdivide,P,sum(sum(P,1),2));
assert(~any(isnan(P(:))));


%
% calculate conditional entropies
Hji = CalcCondEnt(P,E,1,prior);     % H(j|i), 'regular' direction
Hij = CalcCondEnt(P,E,2,prior);     % H(i|j), 'reverse' direction


%
% assertion, remove later
assert(~any(isnan(Hji)));
assert(~any(isnan(Hij)));


%
% sort the edges
M = size(E,1);      % number of edges

    % make big aggregators, keep edges in both directions
       
    Hji = [Hji; Hij];
    
    if (gP.numEntropyBins > 0)
        % 'binning' edge-scores, so that vEdgeList
        % has a random nature; this results in a somewhat
        % random coarse-graph, which helps getting out
        % of local minima
        % binning should be avoided when there's no labeling
        % because in that case the 'optimal' ordering is
        % preferable over a random orderding
        
        Hji = Hji / log(size(U,2));         % normalizing by maximal entropy
        Hji = round(Hji * gP.numEntropyBins);     % bin the entropy score
        Hji = Hji + 0.5 * rand(size(Hji));  % add randomness, avoid bin-mixing
               
    end
      
    [~, idx] = sort(Hji,'ascend');
    
  	idx(idx > M) = -1 * (idx(idx > M) - M);	% remove the offset of revers-
                                % ed edges (whose indices are > M), and
                                % apply a (-) sign to mark them.

                                                             
vEdgeList = idx;




function H = CalcCondEnt(P,E,dim,prior)

%
% calculate the conditional entropy of pairwise terms
%
% INPUT:
%
% P,E   -   pairwise and adjacency matrix of a graphical model
% dim   -   the dimension for which the cond. entropy is calculated
%           if (dim = 1) we calculate H(j|i)
%           if (dim = 2) we calculate H(i|j)
% prior -   estimate of the marginal distribution
%

dim1 = dim;                     % the given dimension
dim2 = not(dim - 1) + 1;        % the "other" dimension

Pr = sum(P,dim2);
H = bsxfun(@rdivide,P,Pr);      % normalize Pr(j|i=li)
H = -1 * H .* log(H);     
H(isnan(H)) = 0;                % NaNs here result from 0*log(0)
H = sum(H,dim2);
% H now holds H(j|i=li) for all 'li' and for all pairwise terms,
% to get H(j|i) we need to take the weighted average:
% H(j|i) = sum_i { p(i = li) * H(j|i=li) }
% ... we do this next:
if isempty(prior)
    % prior is the marginal distribution
    H = squeeze(Pr .* H)';
else
    % use the estimated prior
    if (size(prior,2) > size(P,1))
        % working with histogram, need to normalize
       
        prior = bsxfun(@rdivide,prior(:,1:end-1),prior(:,end));
    end
    H = prior(E(:,dim1),:) .* squeeze(H)';
end
H = sum(H,2);