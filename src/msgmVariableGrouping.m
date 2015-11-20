function [vg, mapFineToCoarse] = msgmVariableGrouping(G, bInitialized)
% msgmVariableGrouping(G, bInitialized) select a variable grouping for the coarsening stage
%
% input:
%
%       - G : graphical model
%       - bBin : boolean flag for whether labels are initialized
%
% output:
%       - vg, struct array
%           vg(i).seed  : the i-th seed
%           vg(i).vars  : variables in the i-th group (including seed)
%           vg(i).edges : (some of the...) edges in the i-th group
%           vg(i).binv  : boolean flag specifying if the edge is 'inversed'
%       
%       - mapFineToCoarse : mapping of variables from a fine scale to a
%                           coarse scale, i.e. mapFineToCoarse(v1) = vc
%                           means that fine-variable indexed by v1 is
%                           mapped to a coarse-variable indexed by vc

% TODO: for fine-scale, pre-processing (scoring) is not necessary!  

    % output data structures
    vg = [];
    mapFineToCoarse = zeros(size(G.u,1),1);

    % VARS and SEEDS
    vbVars = true(size(G.u,1),1);
    vbSeeds = false(size(G.u,1),1);
    
    % local-conditional-entropy scores
    % vbInvEdge defines whether the pairwise is stored
    % in G.p as (v1,v2) or as (v2,v1)
    % TODO: scoreEdges only if G.bProcessed = false
    orderedEdgeList = msgmScoreEdges(G.u, G.adj, G.p, [], [], [], bInitialized);
    vbInvEdge = (orderedEdgeList < 0);
    orderedEdgeList = abs(orderedEdgeList);
    
    % assign SEED variables and their respective group
    iEdge = 0;
    while (any(vbVars) &&  (iEdge < numel(orderedEdgeList)))
       
        % get the next edge
        iEdge = iEdge + 1;
        v1 = G.adj(orderedEdgeList(iEdge), 1);
        v2 = G.adj(orderedEdgeList(iEdge), 2);
        if (vbInvEdge(iEdge))
            % the relevant direction is (v2,v1)
            % TODO: consider transposing the pairwise and getting rid of vg.binv
           
            v_ = v1;
            v1 = v2;
            v2 = v_;
        end
        
        % verify that v2 has not been assigned
        % ..and that v1 is either seed or not assigned
        if (vbVars(v2) && (vbVars(v1) || vbSeeds(v1)))
            
            % set v1 to be v2's seed
            if (vbSeeds(v1))
                % v1 is already a seed variable
                
                v1c = mapFineToCoarse(v1);
                
                % update v1c's group
                vg(v1c).vars = cat(1, vg(v1c).vars, v2);
                vg(v1c).edges = cat(1, vg(v1c).edges, orderedEdgeList(iEdge));
                vg(v1c).binv = cat(1, vg(v1c).binv, vbInvEdge(iEdge));
                
                % update v2's coarse representative
                mapFineToCoarse(v2) = v1c;
                
            else
                % v1 is a new seed variable
                
                % construct v1c's group              
                v1group.seed = v1;
                v1group.vars = [v1; v2];
                v1group.edges = orderedEdgeList(iEdge);
                v1group.binv = vbInvEdge(iEdge);
                vg = cat(1, vg, v1group);
                
                % update v1,v2 coarse representative
                mapFineToCoarse(v1) = numel(vg);
                mapFineToCoarse(v2) = numel(vg);
            end
            
            % update seeds and vars
            vbVars([v1, v2]) = false;
            vbSeeds(v1) = true;
        end
    end
    
    % collect "leftover" variables
    leftoverVars = find(vbVars);
    for i = 1 : numel(leftoverVars)
       
        % the group is singleton
        v1group.seed = leftoverVars(i);
        v1group.vars = leftoverVars(i);
        v1group.edges = [];
        v1group.binv = [];
        vg = cat(1, vg, v1group);
        mapFineToCoarse(leftoverVars(i)) = numel(vg);
    end
end

