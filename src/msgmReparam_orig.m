function G = msgmReparam_orig(G)

%
% ReparamGraph(U,E,P) - reparameterize the unary and pairwise terms of the
% given graphical model (U,E,P).
%
% Reparameterization of the unary and pairwise is an equivalent
% representation of these terms, such that for all labeling assignments the
% energy of the assignment will be equal regardless of which
% parameterization was used to calculate it.
%
% This step turned out to be paramount for the algorithm, as it impacts the
% order by which the edges are selected in SelectEdges(..), and the
% interpolation rule which is set in CompCoarseGraph()
%


U = G.u;
E = G.adj;
P = G.p;

% for new kto2 framework, mark "invalid" edges
idxInf = zeros(size(P,3),1);



%
% for each edge, calculate the mean value of the pairwise,
% both for the right and for the left vertices, this term is
% going to be represented on the unary term instead
pi = mean(P(:,:,~idxInf),2);                         % left vertices
pj = mean(P(:,:,~idxInf),1);                         % right vertices
        % can this be implemented in K^2 instead of 2*K^2 ?


%
% subtract the mean from the pairwise, so that
% the energy of the graph is preserved
P(:,:,~idxInf) = bsxfun(@minus,P(:,:,~idxInf),pi);	% left vertices
P(:,:,~idxInf) = bsxfun(@minus,P(:,:,~idxInf),pj);  % right vertices


%
% ...and add the mean to the unary term
N = size(U,1);              % number of variables
K = size(U,2);              % number of labels
pi = squeeze(pi)';
pj = squeeze(pj)';
if (size(pj,2) == 1)
    % fix dimensions error when there's only a single edge
    pj = pj';
end
    
for i = 1 : K
    % ..for each label separately
    U(:,i) = U(:,i) + accumarray(E(~idxInf,1),pi(:,i),[N,1]);
    U(:,i) = U(:,i) + accumarray(E(~idxInf,2),pj(:,i),[N,1]);
end


G.u = U;
G.p = P;
