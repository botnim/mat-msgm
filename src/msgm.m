function [x, e, t] = msgm(U,E,P,x,params)

%
% TODO:
%   -   consider writing u,adj,p as struct arrays
%   -   test >1 V-cycles, specifically, where labels are initialized
% 
% msgm - main entrypoint to multiscale optimization of graphical models.
%
% Vcycle(U,E,P,L,gP) -  apply the recursive multiscale framework
% iteratively, as in classical MG methods.
%
% INPUT:
%
%   U,E,P,L -   unary term, adjacency matrix, pairwise term and initial
%               labeling, correspondingly. See msmrf(..) for
%               further details.
%
%   gP      -   set of parameters:
%
%
%
% NAMING CONVENTION:
% While the edges aren't directed in general, they are directed in
% practice. This shows up in two cases:
%   (1) The pairwise terms stored in 'P' are generally not symmetric, i.e.
%       P(x,y,:) != P(y,x,:). The vertex that corrsponds to the first
%       dimension of P will usually be called the "left" vertex, and will
%       ususally be denoted using i's. The vertex that corresponds to the
%       second dimension of P will be called the "right" vertex and will
%       usually be denoted with j's.
%   (2) When contracting pairs of vertices, we will have an 'interpolator'
%       vertex, usually denoted using i's, and an 'interpolated' vertex,
%       usually denoted with j's.
%
%

    tic;

    if isempty(params)
        params = setParams;
    end


    %
    % initial reparameteriation of the fine graph
    % can be done outside msmrf() for computational efficiency,
    % otherwise it will be done at each Vcycle iteration
    %[U, P] = msgmReparamGraph(U,E,P,params);

    %
    % Vcycles
    ebest = Inf;
    G.u = U;
    G.adj = E;
    G.p = P;
    G.numLabels = size(G.u,2);
    for i = 1 : params.numV

        x = msgmVcycle(G, x, params);
        e = Energy(G.u, G.adj, G.p, x); 
        e = round(e * 1e6) / 1e6;


    %     if (e < ebest)
    %         % found better labeling assignment
    %         
    %         ebest = e;
    %         Lbest = L;
    %         
    %     end

        assert(e <= ebest);

    end


    %
    % final relaxation of graph
    params.bFineGraph = true;
    if isempty(x)
        % provide initial guess
        [~, x] = min(G.u, [], 2);
    end
    x = msgmRelaxGraph(G.u, G.adj, G.p, x, params);



    %
    % energy & time
    t = toc;
    e = Energy(G.u, G.adj, G.p, x);

end