function vEdgeList = msgmScoreEdges(G, param, bInitialized)
% msgmScpreEdges(G, param, bInitialized) score the edges according to the
% local conditional entropy criterion
%
% output
%
%   - vEdgeList     :   list of edges, ordered according to the local conditional
%                       entropy scores
%   - vbReverseEdge :   flags that the conditional entropy corresponds to a "reverse"
%                       (v2, v1) edge, i.e. to H(v1|v2), and not the "default" direction
%                       (v1, v2) by which the edge's pairwise potential is stored                        
%


    % normalize the unary potential by the variable's connectivity degree
    vDeg = accumarray(G.adj(:), 1, [size(G.u, 1), 1]);
    u = bsxfun(@rdivide, G.u, vDeg);

    % the local energy
    % TODO: reference to paper
    u1(:,1,:) = u(G.adj(:,1),:)';
    u2(1,:,:) = u(G.adj(:,2),:)';
    localEnergy = G.p;
    localEnergy = bsxfun(@plus, localEnergy, u1);
    localEnergy = bsxfun(@plus, localEnergy, u2);
    
    % from local energy to approximate marginal probabilities
    % TODO: reference paper
    localEnergy = bsxfun(@minus, localEnergy, mean(mean(localEnergy,1),2));
    localMarginals = exp(-localEnergy);
    localMarginals = bsxfun(@rdivide, localMarginals, sum(sum(localMarginals,1),2));
    assert(~any(isnan(localMarginals(:))));

    % compute local conditional entropies
    H21 = msgmCondEntropy(localMarginals, 1);	% H(v2|v1), 'regular' direction
    H12 = msgmCondEntropy(localMarginals, 2);	% H(v1|v2), 'reverse' direction   
    Hcond = [H21; H12];
    
    % bin entropy scores
    % TODO: reference to paper
    if ((param.numEntropyBins > 0) && bInitialized)
        
        Hcond = Hcond / log(size(U,2));             % normalizing by maximal entropy
        Hcond = round(Hcond * gP.numEntropyBins);   % bin the entropy score
        Hcond = Hcond + 0.5 * rand(size(Hcond));    % add randomness, avoid bin-mixing
    end
      
    M = size(G.p, 3);
    [~, idx] = sort(Hcond, 'ascend');
  	idx(idx > M) = -1 * (idx(idx > M) - M);                                                         
    vEdgeList = idx;

end


%% Helper

function Hcond = msgmCondEntropy(localMarginals, dim)
% msgmCondEntropy(localMarginals, dim) compute the local conditional
% entropy for to the direction given in 'dim'
%
% input:
%
%   localMarginals  :   the approximate marginal probabilities for pairs of variables
%   dim             :   the direction for which the conditional entropy is computed
%                       H21 is computed when (dim = 1)
%                       H12 is computed when (dim = 2)
%

    
    dim1 = dim;                     % the given dimension
    dim2 = mod(dim1, 2) + 1;        % the "other" dimension
    
    % the conditional probability Pr(x_dim2 | x_dim1)
    p1_marginal = sum(localMarginals, dim2);
    p21_cond = bsxfun(@rdivide, localMarginals, p1_marginal);

    % the conditional entropy
    Hcond = p21_cond .* log(p21_cond);     
    Hcond(isnan(Hcond)) = 0;	% NaNs here result from 0*log(0)
    Hcond = -1 * sum(Hcond, dim2);
    Hcond = p1_marginal .* Hcond;
    Hcond = sum(Hcond, dim1);
    Hcond = Hcond(:);
    
    % NaNs indicate bugs
    assert(~any(isnan(Hcond)));
end