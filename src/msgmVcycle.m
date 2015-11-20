function x = msgmVcycle(G, x, param)
% msgmVcycle(G, x, param) the recursive construction and optimization of
% the multiscale hierarchy
%
% input:
%
%   G   -   graphical model (see msgm() for full description)
%
%   x   -   an initial labeling assignment to the variables (column vector,
%           empty if an initial guess is not available)
%
%   param  -   set of parameters, see msgmParams().
%
%
% output:
%
%   x   -   a labeling assignment of the variables
%

    % check stopping condition
    % TODO: make sure LSA forground is not messing up
    if (size(G.u,1) <= param.numMinVars)
        
        if (isempty(x))
            % labels not initialized,
            % current move-making methods require initialization
           
            % initialize according to unary term
            [~, x] = min(G.u, [], 2);
        end
        x = msgmOptimizeScale(G, x, param);
        
        return;
    end


    % reparameterize the graph
    % TODO: take out of here!
    G1 = msgmReparam_orig(G);
    G2 = msgmReparam(G);

    % run inference on the current scale
    % TODO: it may work better if the graph is first processed (read:
    % reparamterized)
    x = msgmOptimizeScale(G, x, param);
  

    % coarsen the graph
    % TODO: remove vg if compatible relaxations is unecessary,
    %       otherwise, come up with more elegant solution
    [Gc, xc, mapFineToCoarse, mapInterpolation, vg] = msgmCoarsening(G, x);
    msgmEnergyAssert(G, x, Gc, xc);


    % recursive call
    xc = msgmVcycle(Gc, xc, param);


    % interpolate solution
    x = msgmInterpolate(G, vg, xc, mapFineToCoarse, mapInterpolation, param);

    
    % run inference on the current scale
    x = msgmOptimizeScale(G, x, param);


end