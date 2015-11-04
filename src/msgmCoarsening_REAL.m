function [Gc, xc, mapFineToCoarse, mapInterpolation] = msgmCoarsening_REAL(G, x)
%msgmCoarsening apply coarsening procedure for given variable-grouping
%
% output:

    
    % select a variable-grouping
    [vg, mapFineToCoarse] = msgmVariableGrouping(G);

    % set an interpolation rule
    % TODO: reference to equations
    mapInterpolation = msgmSetInterpolationRule(G, x, vg);
    
    % set the coarse potentials
    Gc = msgmSetCoarsePotentials(G, vg, mapFineToCoarse, mapInterpolation);
    
    % intialize coarse scale's labeling
    xc = msgmInherit(x);

end


%% Coarsening subroutines

function map = msgmSetInterpolationRule(G, x, vg)
% msgmSetInterpolationRule set an interpolation-rule, from labeling of a
% coarse scale to labeling of a fine scale

    % intialize the interpolation-rule table
    map = zeros(size(G.u,1), G.numLabels);
    
    % itearte seed variables
    for i = 1 : numel(vg)
               
        % iterate variables in the seed's group
        for j = 1 : numel(vg(i).vars)
           
            if (vg(i).vars(j) == vg(i).seed)
                % seed variable maps to itself
                % TODO: reference to paper
                
                map(vg(i).vars(j),:) = 1 : G.numLabels;
                
            else
                % TODO: look for bugs in symmetry here!!
                % find the minimizer (Eq. (3))
                pairwise = squeeze(G.p(:,:,vg(i).edges(j-1)));
                if (~vg(i).binv(j-1))
                    % transpose the pairwise s.t. seed is on 2nd dim

                    pairwise = pairwise';
                end
                pairwise = bsxfun(@plus, pairwise, G.u(vg(i).vars(j),:)');
                [~, map_] = min(pairwise,[],1);            

                if any(x)
                    % labels are intiailized,
                    % reset the interpolatin rule (Eq. (6))

                    map_(x(vg(i).seed)) = x(vg(i).vars(j));
                end
                
                map(vg(i).vars(j),:) = map_;
            end
        end
    end
end

function Gc = msgmSetCoarsePotentials(G, vg, mapFineToCoarse, mapInterpolation)

    % keep track which edges have been accounted for
    vbTouch = false(size(G.p,3), 1);

    % set the coarse unary terms (Eq. (4))
    uc = zeros(numel(vg), G.numLabels);
    for i = 1 : numel(vg)
        
        % TODO: consider arrayfun here
        % TODO: consider general case when numLabels is not fixed
        uGroup = zeros(1, G.numLabels);
        
        % sum j's unary term, considering the interpolation
        % this is the first term in Eq. (4)        
        for j = 1 : numel(vg(i).vars)
            
            vj = vg(i).vars(j);
            uGroup = uGroup + G.u(vj, mapInterpolation(vj,:));            
        end
        
        % sum the pairwise energy, considering the interpolation
        % this is the second term in Eq. (4)
        for j = 1 : numel(vg(i).edges)
            
            vbTouch(vg(i).edges(j)) = true;
            pairwise = G.p(:,:,vg(i).edges(j));
            
            % interpolation rule applied to v1, v2
            v1 = G.adj(vg(i).edges(j), 1);
            v2 = G.adj(vg(i).edges(j), 2);
            map1 = mapInterpolation(v1,:);
            map2 = mapInterpolation(v2,:);
            
            uGroup = uGroup + diag(pairwise(map1, map2))';
        end
        
        uc(i,:) = uGroup;
    end
    
    % set the coarse pairwise terms (Eq. (5))
    adjMat = zeros(numel(vg));    % coarse scale's adjacency matrix
    pc = [];
    adjc = [];
    vNoTouch = find(~vbTouch);
    for i = 1 : numel(vNoTouch)

        iEdge = vNoTouch(i);
              
        % coarse representatives of edge's variables
        v1 = G.adj(iEdge,1);
        v2 = G.adj(iEdge,2);
        v1c = mapFineToCoarse(v1);
        v2c = mapFineToCoarse(v2);
        map1 = mapInterpolation(v1,:);
        map2 = mapInterpolation(v2,:);
        
        % the pairwise term, considering the interpolation
        pairwise = G.p(:,:,iEdge);
        pairwise = pairwise(map1, map2);
        
        % check if the edge defines a self-loop
        if (v1c ~= v2c)

            % check if (v1c,v2c) is already defined on the coarse scale
            if (~adjMat(v1c,v2c))
                % append as a new edge
                
                % track changes to adjacency
                adjMat(v1c,v2c) = size(adjc,1) + 1;     % pairwise stored as (v1c,v2c)
                adjMat(v2c,v1c) = -(size(adjc,1) + 1);  % flag that pairwise is transposed
                
                % update adjacency and pairwise
                adjc = cat(1, adjc, [v1c, v2c]);
                pc = cat(3, pc, pairwise);

            else
                % add to an existing edge
                
                idx = adjMat(v1c,v2c);
                if (idx < 0)
                    pairwise = pairwise';
                end
                pc(:,:,abs(idx)) = pc(:,:,abs(idx)) + pairwise;
            end
            
        else
            % self-loop, add the edge to respective coarse unary term
            
            uc(v1c,:) = uc(v1c,:) + diag(pairwise)';
        end
    end
    
    Gc.u = uc;
    Gc.p = pc;
    Gc.adj = adjc;
    Gc.numLabels = G.numLabels;
    Gc.bProcessed = false;
end
