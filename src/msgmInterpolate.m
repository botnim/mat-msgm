function x = msgmInterpolate(G, vg, xc, mapFineToCoarse, mapInterpolation, param)
% msgmInterpolate(xc, mapFineToCoarse, mapInterpoation)
% given a labeling 'xc' of the coarse scale, interpolate the assignment
% according to the interpolation rule 'mapInterpolation'
%
   
    % all variables in a group get the label of their seed
    x = xc(mapFineToCoarse);

    if (~param.bSoftInterpolation)
        % apply the interpolation rule for all variables
    
        inds = sub2ind(size(mapInterpolation), ...
                (1 : numel(mapFineToCoarse))', ...
                x);
        x = mapInterpolation(inds);
    
    else
        % 'soft' interpolation:
        % fix labels of seed variables and optimize labels of all the rest

        % condition upon the seed variables
        vb = false(size(G.u, 1), 1);
        vb([vg.seed]) = true;
        [Gcond, ~] = msgmConditionalDist(G, x, vb);
        
        % optimize labels of the conditional graph
        % TODO: numRelax was set to 5 here
        % TODO: initialize xcond 'winner takes all' MOVE TO RELAX GRAPH! 
        [~, xcond] = min(Gcond.u, [], 2);
		xASSERT = x;
		xASSERT(~vb) = xcond;
        xcond = msgmOptimizeScale(Gcond, xcond, param);
        
        % map the optimized labels to the original variables
        x(~vb) = xcond;
		
		assert(msgmEnergy(G, x) <= msgmEnergy(G, xASSERT));
    
    end
end