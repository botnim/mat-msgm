function [L, e, t, Q] = msmrf(U,E,P,L,gP)

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

if isempty(gP)
    gP = setParams;
end


%
% initial reparameteriation of the fine graph
% can be done outside msmrf() for computational efficiency,
% otherwise it will be done at each Vcycle iteration
[U, P] = ReparamGraph(U,E,P,gP);


%
% set prior
prior = SetPrior(size(U,1),size(U,2),gP);


% variables to keep track of best labeling/energy
Lbest = [];
ebest = Inf;


%
% Vcycles
for i = 1 : gP.numV
    
    [L, prior, Q(i+1)] = Vcycle(U,E,P,L,prior,gP,true,0);
    e = Energy(U,E,P,L); 
    e = round(e * 1e6) / 1e6;


    if (e < ebest)
        % found better labeling assignment
        
        ebest = e;
        Lbest = L;
        
    end
    
%     assert(e <= ebest);

end
L = Lbest;



%
% final relaxation of graph
gP.bFineGraph = true;
if isempty(L)
    % provide initial guess
    [~, L] = min(U,[],2);
end
L = RelaxGraph(U,E,P,L,gP);



%
% energy & time
t = toc;
e = Energy(U,E,P,L);

end